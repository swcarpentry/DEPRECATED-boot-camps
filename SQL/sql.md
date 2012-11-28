Databases using SQL
===================

**If you don't have Firefox installed already please install it now**

Relational databases
--------------------

* Relational databases store data in tables and with fields (columns) and records (rows)
* Data in tables has types, just like in Python, and all values in a field have the same type
* Queries let us look up data or make calculations based on columns
* The queries are distinct from the data, so if we change the data we can just rerun the query

Getting setup
-------------

1. Install Firefox
2. Install the SQLite Manager add on **Tools -> Add-ons -> Search -> SQLite Manager -> Install -> Restart**
3. Download the Portal Database **TODO**
4. Open SQLite Manage **Tools -> SQLite Manager**

The data
--------

This is data on a small mammal community in southern Arizona over the last 35 years.

This is a real dataset that has been used in over 100 publications.
I've simplified it just a little bit for the workshop,
but you can download the [full dataset](

Import
------

1. Start a New Database **Database -> New Database**
2. Start the import **Database -> Import**
3. Select the file to import
4. Give the table a name (or use the default)
5. If the first row has column headings, check the appropriate box
6. Make sure the delimiter and quotation options are correct
7. Press **OK**
8. When asked if you want to modify the table, click **OK**
9. Set the data types for each field

***Exercise: Import the plots and species tables***

Basic queries
-------------

Let’s write an SQL query that selects only the scientific names from the species table.

    SELECT scientific_name FROM species;

We have capitalized the words SELECT and FROM because they are SQL keywords.
SQL is case insensitive, but it helps for readability – good style.

If we want more information, we can just add a new column to the list of fields, right after SELECT:

    SELECT scientific_name, species_id FROM species;

Or we can select all of the columns in a table using the wildcard *

    SELECT * FROM species;

***Exercise: Write a query that returns only the taxa column***

### Unique values

So, we've all written a query that pulls out the taxa column from the database,
but what if we want only the unique values so that we can quickly see what
taxa have been sampled.

    SELECT DISTINCT taxa FROM species;

If we select more than one column, then the distinct pairs of values are returned

    SELECT DISTINCT taxa, scientific_name FROM species;

Calculated values
-----------------
Now let's switch to the **surveys** table.
Here we have data on every individual that was captured at the site,
including when they were captured, what plot they were captured on,
their sex and their weight in grams.

Now that we're using a table with numbers in it we can see that in addition
to selecting columns we can also do calculations with their values.
For example, if we wanted the mass in kg instead of g we would use

    SELECT plot, species, sex, wgt, wgt / 1000.0 from surveys

When we run the query, the expression ``wgt / 1000.0`` is evaluated for each row
and appended to that row, in a new column. 
Expressions can use any fields, any arithmetic operators (+ - * /)
and a variety of built-in functions (). For example, we could round the values to
make them easier to read.

    SELECT plot, species, sex, wgt, ROUND(wgt / 1000.0, 2) FROM surveys;

Filtering
---------
One of the most powerful features of a database is the ability to filter data –
selecting only the data meeting certain criteria.
For example, let’s say we only want data for the species Dipodomys merriami,
which has a species code of DM.
We need to add a WHERE clause to our query:

    SELECT * FROM surveys WHERE species="DM";

We can do the same thing with numbers.
Here, we only want the data since 2000:

    SELECT * FROM surveys WHERE year >= 2000;

We can user more sophisticated conditions by combining tests with AND and OR.
For example, suppose we want to data on Dipodomys merriami startinging in the year 2000:

    SELECT * FROM surveys WHERE (year >= 2000) AND (species = "DM");

Note that the parentheses aren’t needed, but again, they help with readability.
They also ensure that the computer combines AND and OR in the way that we intend.

If we wanted to get data for any of the Dipodomys species,
which have species codes DM, DO, and DS we could combine the tests using OR:

    SELECT * FROM surveys WHERE (species = "DM") OR (species = "DO") OR (species = "DS");

***Exercise: Write a query that returns the day, month, year, species ID, and weight
for individuals caught on plot 1 that weigh more than 75 grams***

Exporting results of queries
----------------------------
Getting the result of your query out to work with elsewhere is as easy
as clicking the **Actions** button and choosing **Save Result to File**.

Building more complex queries
-----------------------------

Now, lets combine the above queries to get data for the 3 Dipodomys species
from the year 2000 on.
This time, let’s use IN as one way to make the query easier to understand.
It is equivalent to saying ``WHERE (species = "DM") OR (species = "DO") OR (species = "DS")``,
but reads more neatly:

    SELECT * FROM surveys WHERE (year >= 2000) AND (species IN ("DM", "DO", "DS"));

We started with something simple, then added more clauses one by one,
testing their effects as we went along.
For complex queries, this is a good strategy, to make sure you are getting what you want.
Sometimes it might help to take a subset of the data that you can easily see in a temporary
database to practice your queries on before working on a larger or more complicated database.

Sorting
-------

We can also sort the results of our queries by using ORDER BY.
For simplicity, let’s go back to our species table and alphabetize it by Genus.

    SELECT * FROM species ORDER BY Genus ASC;

The keyword ASC tells us to order it in Ascending order.
We could alternately use DESC to get descending order.

    SELECT * FROM species ORDER BY Genus DESC;

ASC is the default.

We can also sort on several fields at once.
To truly be alphabetical, we might want to order by genus then species.

    SELECT * FROM species ORDER BY genus ASC, species ASC;

***Exercise: Write a query that returns all of the data in the plots table,
sorted alphabetically by plot type and then (within each plot type),
in descending order of the plot ID.***

Order of execution
------------------

Another note for ordering:
We don’t actually have to display a column to sort by it.
For example, let’s say we want to order by the species ID, but we only want to see genus and species. 

    SELECT genus, species FROM species ORDER BY species_id ASC;

We can do this because sorting occurs earlier in the computational pipeline than field.

The computer is doing this:

1. Filtering rows according to WHERE
2. Sorting results according to ORDER BY
3. Displaying requested columns or expressions.


***Exercise: Let’s try to combine what we’ve learned so far in a single query.  Let’s go back to the surveys table and lets say that we want to display
the three date fields, species ID, and weight in kilograms (rounded to two 
decimal places),  for rodents captured in 1999, ordered alphabetically by 
the species ID.***

The order of the clauses is dictated by SQL: SELECT, FROM, WHERE, ORDER BY
and we often write each of them on their own line for readability.

Aggregation
-----------
Aggregation allows us to combine results group records based on value
and calculate combined values in groups.

Let’s go to the surveys table and find out how many individuals there are.
Using the wildcard simply counts the number of records (rows)

    SELECT COUNT(*) FROM individuals

We can also find out how all of those individuals weigh.

    SELECT SUM(wgt) FROM individuals

***Do you think you could output this value in kilograms,
rounded to 3 decimal places?***

    SELECT ROUND(SUM(wgt)/1000.0, 3) FROM surveys

There are many other aggregate functions included in SQL including
MAX, MIN, and AVG.
 
***From the surveys table, can you use one query to output the total weight, average weight, and the min and max weights?***

Now, let's try to see how many individuals were counted in each species?

    SELECT species, COUNT(*) FROM individuals

Why does this NOT work?

The database counted all the individuals,
but we didn’t specify anything about the species,
so it just picked an arbitrary value and returned it.

If we select fields and aggregate at the same time, values from 
unaggregated fields can be any value – we need to tell the data how to 
aggregate, and we can do that using GROUP BY clause

    SELECT species, COUNT(*)
    FROM surveys
    GROUP BY species

***Exercise: How many individuals were counted in each year?***
SELECT month, COUNT(DISTINCT sp_code) FROM individuals GROUP BY month

***Exercise: How many individuals were counted in each species in each year?***

We can order the results of our aggregation by a specific column, 
including the aggregated column.
Let’s count the number of individuals of each species captured,
ordered by the count

SELECT species, COUNT(*)
FROM surveys
GROUP BY species
ORDER BY COUNT(sp_code)

***Exercise: Write a query that let's us look at which years contained the most individuals and which had the least?***

***Exercise: Write a query that shows us which species had the largest individuals on average?***

Database Design
---------------
Each field in a database should store a single value.
Information should not be duplicated in a database.
Each table should be about a single subject (avoids unnecessary replication).
When naming fields, you should think about meaning, not presentation.

When we divide our data between several tables,
we need a way to bring it back together again.
The key is to have an identifier in common between tables – shared columns.
This will allow us to join tables, which we will discuss now.

For example, the species ID is included in the surveys table,
but we don’t know what the species ID stands for.
That information is stored in the Species table and can be
linked to if we need it.
This means that we don't have to record the full genus, species,
and taxa information for the several thousand individuals of each species. 

Joins
-----
To combine data from two tables we use the SQL JOIN command,
which comes after the FROM command.

    SELECT * FROM surveys JOIN species

But this didn’t do what we wanted because we didn’t tell the database 
manager how the tables are related to each other.
Look at the number of rows it returned!
Simply adding the JOIN combines every row of 
one table with every row of the other -
it creates a cross-product of the sets of rows.
We need to tell SQL how the tables are related.

To do this we indicated which columns provide the link between
the two tables using the word ON.
What we want is to join the data with the same species codes.

    SELECT *
    FROM surveys
    JOIN species ON surveys.species = species.species_id

ON is like WHERE, it filters things out according to a test condition.
We use the table.colname to tell the manager what column in which table
we are referring to.

We often won't want all of the fields from both tables,
so anywhere we would have used a field name in a non-join query,
we can use *table.colname*

For example, what if we wanted information on when individuals of each
species were captured, but instead of their species ID we wanted their
actual species names.

    SELECT surveys.year, surveys.month, surveys.day, species.genus, species.species
    FROM surveys
    JOIN species ON surveys.species = species.species_id

***Exercise: Write a query that the genus, the species, and the weight of every individual captured at the site***

Joins can be combined with sorting, filtering, and aggregation.
So, if we wanted average mass of the individuals on each different
type of treatment, we could do something like

    SELECT plots.plot_type, AVG(surveys.wgt)
    FROM surveys
    JOIN plots
    ON surveys.plot = plots.plot_id
    GROUP BY plots.plot_type



Creating tables
---------------

Adding data to existing tables
------------------------------
