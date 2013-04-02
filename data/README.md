# Data management

Mike Jackson, The Software Sustainability Institute. 

**Based on materials by Greg Wilson.**

## What is a database?

A *database* is a system designed to manage complex data in a systematic, structured way. It allows you to store, update, retrieve and delete data held within.

You can retrieve data from a database by issuing commands, called *queries*. A query specifies what data you are interested in. The database takes care of how to get it. Queries allow data to be extracted, filtered, combined, summarised and analysed in various ways. This makes databases a good choice when you are exploring data. 

If you're used to using spreadsheets to manage your data, running queries can be viewed as akin to inserting formulae or new sheets for your analyses. Unlike spreadsheets, databases do not typically have built in charting and visualising tools. But, it's always possible to export the results of a query to be used in other tools.

However, one advantage of databases is that thousands of programmer-years have gone into their design and implementation and some have evolved to satisfy the real-time demands of commercial environments. As a result, databases are designed so that they can handle very large datasets - terabytes or more - quickly and reliably, far beyond the capabilities of a spreadsheet.

## What is a relational database?

The most common type of database is a [relational database](http://en.wikipedia.org/wiki/Relational_database) (or RDBMS). A relational database stores data as tables which represent relations e.g. Experiments, Observations, Projects. The columns of each table represent fields or properties of the data e.g. name, address, post-code. The rows, or records, contain the data itself. For example,

    --------------------------------------------------
    | ID           | Project | Experiment ID | Hours |
    --------------------------------------------------
    | fred_boggs   | ATLAS   | 1             | 3.5   |
    | fred_boggs   | CMS     | 4             | 8.2   |
    | julie_hickey | ATLAS   | 2             | 7.7   |
    | julie_hickey | CMS     | 1             | 18.2  |

Multiple tables may share the same column names, and this allows data from multiple tables to be combined.

Queries are written in [SQL](http://en.wikipedia.org/wiki/SQL) (Structured Query Language), an expressive language to express complex data extraction, filtering, aggregation, sorting and combination tasks,

    SELECT ID, SUM(Hours) FROM Experiment GROUP BY ID;

Many popular databases including MySQL, Postgres, IBM DB2 and Microsoft SQL Server are relational databases.

For a detailed introduction to relational databases, see the [Software Carpentry online lecture](http://software-carpentry.org/4_0/databases/intro.html).

## What is a NoSQL database?

In recent years, [NoSQL](http://en.wikipedia.org/wiki/NoSQL) databases have been proposed as an alternative to relational databases (NoSQL means "not only SQL"). Instead of tables, NoSQL databases store collections of 'things'.

NoSQL databases are designed for the storage and retrieval of large volumes of data, data which may be unstructured, messy or unpredictable. where the relationships between individual data items may not be of so much concern, or interest, as properties across the data as a whole.

Qualities of NoSQL databases include,

* Simplicity
* Lightweight
* Scalablility, at low cost
* Availability

Applications can NoSQL databases include,

* Logging Twitter posts or the internet server logs from a large numbers of users.
* Logging sensor data at periodic intervals.
* Recording observations from researchers in the field.
* Recording particle collision data in the [Muon Ion Cooling Experiment](http://mice.iit.edu/)
* Recording data with incomplete or unclear meta-data in the [Compact Muon Solenoid](http://cms.web.cern.ch) experiment.

Are NoSQL databases better than relational databases? Yes! No! It depends upon your data and what you want to do with it. As always, in software development, there are no one-size-fits-all solution.

## Types of NoSQL database

There are 3 main classes of NoSQL database:

* [Graph databases](http://en.wikipedia.org/wiki/Graph_database) in which each node is a data item, and edges represent relationships between the data. Both data and relationships are expressed via key-value properties. This is the essence of [linked data](http://en.wikipedia.org/wiki/Linked_data).
* [Key-value store](http://en.wikipedia.org/wiki/Key/value_store) in which opaque data objects (the values) are each associated with unique keys which allows the data to be very efficiently inserted, retrieved and deleted. If you're familiar with programming concepts such as association lists, hashtables or Python dictionaries, then key-value stores can be viewed as analogous but with the data persisted in a database.
* [Document-oriented databases](http://en.wikipedia.org/wiki/Document-oriented_database) which store structured documents, which may be in XML or JSON (which we'll look at soon), sometimes flat, sometimes as a heirarchical set of collections, and which allow queries based on the content of these documents, and also for them to be updated. The [CMS experiment uses a document-oriented database](http://readwrite.com/2010/08/26/lhc-couchdb).

## What is MongoDB?

[MongoDB](http://www.mongodb.org/) is a document-oriented database.

* A MongoDB server holds 0 or more databases. 
* Each database holds 1 or more collections.
* Each collection 0 or more documents. 

MongoDB documents are written in [JSON](http://json.org/) is JavaScript Object Notation, a simple notation for sstructured data. A document is defined as a collection of key-value pairs. Each key is a unique string and each value can be either,

* A string,
* A number,
* 'true' or  'false',
* An ordered list of values,
* A collection of key-value pairs.

In fact, MongoDB documents are not JSON but [BSON](http://bsonspec.org/) or Binary-encoded JSON, as they can also store dates and binary data. This allows MongoDB to store data as BSON documents or embedded within BSON documents, for example if you had binary data, the BSON document can record this data with additional meta-data expressed using the standard JSON notation.

MongoDB is *schema free* - the documents can be all of the same structure or of different structures. Mongo 

### pymongo

Many databases support command-line and graphical user interfaces. However, they also support, or have developed for them, drivers, so they can be accessed and used from within various programming languages.

MongoDB has [drivers](http://api.mongodb.org/) in many languages including C, C++, Java, Perl and Python. We'll be using [pymongo](http://api.mongodb.org/python/current/tutorial.html).

## Let's dive in...

Next: [Connecting to our database, browsing and selecting documents](Select.md)
