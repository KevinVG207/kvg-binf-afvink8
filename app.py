import datetime
import math

import matplotlib.pyplot as plt
import mysql.connector
import numpy as np
from Bio import Entrez
from Bio import Medline
from flask import Flask, render_template, request

app = Flask(__name__)


@app.route('/')
def homepage():
    return render_template("index.html")


@app.route('/ensembl', methods=["POST", "GET"])
def ensembl():
    if request.method == "POST" and request.form.get("search", "") != "":
        # POST request
        rows = []
        search = request.form.get("search")
        try:
            connection = mysql.connector.connect(host='ensembldb.ensembl.org',
                                                 port=3306,
                                                 user='anonymous',
                                                 db='homo_sapiens_core_95_38')
            cursor = connection.cursor()
            query = """select * from gene where description like '%""" + search + """%' limit 100;"""
            cursor.execute(query)
            out_rows = cursor.fetchall()
            cursor.close()
            connection.close()

            for row in out_rows:
                rows.append(row[9].replace(search, "<b>" + search + "</b>"))

        except mysql.connector.Error as error:
            rows = [error]
        return render_template("ensembl.html", rows=rows)
    else:
        return render_template("ensembl.html", rows=[])


@app.route('/pubmed1', methods=["POST", "GET"])
def pubmed1():
    if request.method == "POST" and request.form.get("search", "") != "":
        terms = request.form.get("search").split(",")
        all_years = []
        now = datetime.datetime.now()
        min_year = now.year
        for term in terms:

            Entrez.email = "kevinvg207@gmail.com"
            term = term.strip()
            count = 0
            handle = Entrez.egquery(term=term)
            record = Entrez.read(handle)
            for row in record["eGQueryResult"]:
                if row["DbName"] == "pubmed":
                    count = row["Count"]
            handle.close()

            print(term, count)

            search_handle = Entrez.esearch(
                usehistory="y",
                db="pubmed",
                term=term,
                retmax=count
            )
            search_results = Entrez.read(search_handle)
            search_handle.close()
            id_list = search_results["IdList"]
            webenv = search_results["WebEnv"]
            query_key = search_results["QueryKey"]

            batch_size = 500
            count = len(id_list)
            records = []
            for start in range(0, count, batch_size):
                end = min(count, start + batch_size)
                print("Going to download record %i to %i" % (start + 1, end))
                medline_handle = Entrez.efetch(
                    db="pubmed",
                    rettype="medline",
                    retmode="text",
                    retstart=start,
                    retmax=batch_size,
                    webenv=webenv,
                    query_key=query_key
                )
                medline_records = Medline.parse(medline_handle)
                records = records + list(medline_records)
                medline_handle.close()

            years = []
            for record in records:
                year = now.year
                try:
                    year = int(record.get("DP", "?").split(" ")[0])
                except ValueError:
                    print("incorrect date format, skipping")
                    print(record.get("DP", "?"))
                    pass
                years.append(year)
            all_years.append(years)

        # Generating the plot. This became a big mess help
        # I started with a plot for 1 term and tried to make it work for
        # multiple terms but that created a big mess.
        # Should probably redo this but can't be bothered.
        for group in all_years:
            if group and min(group) < min_year:
                min_year = min(group)
        min_year = min_year - (min_year % 5)
        year_diff = now.year - min_year
        N = math.ceil(year_diff / 5)  # -1?
        labels = []
        for i in range(N):
            cur_year = min_year + 5 * i
            labels.append(f"{cur_year} - {cur_year + 4}")
        fig, ax = plt.subplots()
        ind = np.arange(N)
        width = 0.5
        bars = []
        groups = len(all_years)
        for group in all_years:
            grouped_years = []
            for i in range(N):
                grouped_years.append(0)
            for year in group:
                grouped = math.floor((year - min_year) / 5)
                grouped_years[grouped] += 1
            index = all_years.index(group)
            bars.append(ax.bar(ind - width/groups + width/groups * index, grouped_years, width/groups, label=terms[index].strip()))
        ax.axhline(0, color='grey', linewidth=0.8)
        ax.set_ylabel('Counts')
        ax.set_title('Counts per 5 years')
        ax.set_xticks(ind)
        ax.set_xticklabels(labels)
        ax.legend()

        for tick in ax.xaxis.get_major_ticks():
            tick.label.set_fontsize(5)
            tick.label.set_rotation("vertical")

        # Dynamic file name because browsers keep previous image cached
        base_url = "static/images/plot/"
        datestamp = str(now.year) + str(now.month) + str(now.day) + str(now.hour) + str(now.minute) + str(now.second)
        full_url = base_url + datestamp + ".png"
        plt.savefig(full_url)

        return render_template("pubmed1.html", plot=True, plot_url=full_url)
    else:
        return render_template("pubmed1.html", plot=False)


if __name__ == '__main__':
    app.run()
