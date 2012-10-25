import sqlite3

connection = sqlite3.connect("experiments.db")

cursor = connection.cursor()

cursor.execute("SELECT FirstName, LastName FROM Person;")

results = cursor.fetchall()

for r in results:
    print r[0], r[1]

cursor.close()
connection.close()
