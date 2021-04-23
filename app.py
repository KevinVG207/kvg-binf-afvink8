from flask import Flask, render_template, request
from flask_mysqldb import MySQL
import mysql.connector

app = Flask(__name__)


@app.route('/', methods=["POST", "GET"])
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
                rows.append(row[9].replace(search, "<b>"+search+"</b>"))

        except mysql.connector.Error as error:
            rows = [error]
        return render_template("ensembl.html", rows=rows)
    else:
        return render_template("ensembl.html", rows=[])


if __name__ == '__main__':
    app.run()
