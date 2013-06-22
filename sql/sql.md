Databases using SQL
===================

This material is provided with a [Creative Commons Attribution licence CC-BY 3.0](http://creativecommons.org/licenses/by/3.0/us/).
It was originally developed for a Software Carpentry Bootcamp by Ethan White and has been modified by Sarah Supp.

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
3. Download the [Portal Database](https://github.com/swcarpentry/2012-11-UNC/raw/master/SQL/portal_mammals.sqlite)
4. Open SQLite Manager **Firefox Button -> Web Developer -> SQLite Manager** 
    or **Tools -> SQLite Manager**

The data
--------

This is data on a small mammal community in southern Arizona over the last 35 years.
This is part of a larger project studying the effects of rodents and ants
on the plant community.
The rodents are sampled on a series of 24 plots, with different 
experimental manipulations of which rodents are allowed to access the plots.

This is a real dataset that has been used in over 100 publications.
I've simplified it just a little bit for the workshop,
but you can download the [full dataset](http://esapubs.org/archive/ecol/E090/118/)
and work with it using exactly the same tools we'll learn about today.

Import
------

Data can be added that is already in a sqlite databae, or by entering CSV or TXT files manually.

1. Start a New Database **Database -> New Database**
2. Start the import **Database -> Import**
3. Select the file to import
4. Give the table a name (or use the default)
5. If the first row has column headings, check the appropriate box
6. Make sure the delimiter and quotation options are correct
7. Press **OK**
8. When asked if you want to modify the table, click **OK**
9. Set the data types for each field

***EXERCISE: Import the plots and species tables***

***  Explore the different options for displaying the data in the SQL Manager**

You can also use this same approach to append new data to an existing table.

Basic queries
-------------
Let's start by using the **surveys** table.
Here we have data on every individual that was captured at the site,
including when they were captured, what plot they were captured on,
their species ID, sex and weight in grams.

Let’s write an SQL query that selects only the year column from the surveys table.

    SELECT year FROM surveys;

We have capitalized the words SELECT and FROM because they are SQL keywords.
SQL is case insensitive, but it helps for readability – good style.

If we want more information, we can just add a new column to the list of fields, right after SELECT:

    SELECT day, month, year FROM surveys;

Rows and columns aren't stored in any particular order, they are just displayed that way. 
We can rearrange columns by changing the order in our select statement:

    SELECT year, month, day FROM surveys;
    
Or we can even ask it to repeat columns (although that may not be very useful):

    SELECT year year FROM surveys;

Or we can select all of the columns in a table using the wildcard *

    SELECT * FROM surveys;
    
***EXERCISE: write a query that will return the species column from the database***

### Unique values

**How can we remove duplicate values from our query?**

So, we've all written a query that pulls out the species column from the database,
but what if we want only the unique values so that we can quickly see what
species have been sampled.

    SELECT DISTINCT species FROM surveys;

If we select more than one column, then the distinct pairs of values are returned

    SELECT DISTINCT year, species FROM surveys;
    
***EXERCISE: write a query that returns the months that were sampled in each year***

### Calculated values

We can also do calculations with the values in a query.
For example, if we wanted to look at the mass of each individual
on different dates, but we needed it in kg instead of g we would use

    SELECT month, day, year, wgt/1000.0 from surveys

When we run the query, the expression ``wgt / 1000.0`` is evaluated for each row
and appended to that row, in a new column. 
Expressions can use any fields, any arithmetic operators (+ - * /)
and a variety of built-in functions (). For example, we could round the values to
make them easier to read.

    SELECT plot, species, sex, wgt, ROUND(wgt / 1000.0, 2) FROM surveys;

***EXERCISE: Write a query that returns
             The day, month, year, speciesID and weight in mg***

Filtering
---------
One of the most powerful features of a database is the abilinmmty to filter data –
selecting only the data meeting certain criteria.
For example, let’s say we only want data for the species Dipodomys merriami,
which has a species code of DM.
We need to add a WHERE clause to our query:

    SELECT * FROM surveys WHERE species="DM";

We can do the same thing with numbers.
Here, we only want the data since 2000:

    SELECT * FROM surveys WHERE year >= 2000;

We can user more sophisticated conditions by combining tests with AND and OR.
For example, suppose we want to data on Dipodomys merriami starting in the year 2000:

    SELECT * FROM surveys WHERE (year >= 2000) AND (species = "DM");

Note that the parentheses aren’t needed, but again, they help with readability.
They also ensure that the computer combines AND and OR in the way that we intend.

If we wanted to get data for any of the Dipodomys species,
which have species codes DM, DO, and DS we could combine the tests using OR:

    SELECT * FROM surveys WHERE (species = "DM") OR (species = "DO") OR (species = "DS");

***EXERCISE: Write a query that returns
   The day, month, year, species ID, and weight (in kg) for
   individuals caught on plot 1 that weigh more than 0.075 kg***

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

    SELECT *
	FROM surveys
	WHERE (year >= 2000) AND (species IN ("DM", "DO", "DS"));

We started with something simple, then added more clauses one by one,
testing their effects as we went along.
For complex queries, this is a good strategy, to make sure you are getting what you want.
Sometimes it might help to take a subset of the data that you can easily see in a temporary
database to practice your queries on before working on a larger or more complicated database.

Sorting
-------

We can also sort the results of our queries by using ORDER BY.
For simplicity, let’s go back to our species table and alphabetize it by species.

    SELECT * FROM species ORDER BY species ASC;

The keyword ASC tells us to order it in Ascending order.
We could alternately use DESC to get descending order.

    SELECT * FROM species ORDER BY species DESC;

ASC is the default.

We can also sort on several fields at once.
To truly be alphabetical, we might want to order by genus then species.

    SELECT * FROM species ORDER BY genus ASC, species ASC;

***Exercise: Write a query that returns
             The genus, species, and taxon, sorted alphabetically by taxon.***
             
***Exercise: Write a query that returns the genus, species of rodents, sorted alphabetically by scientific name***

Order of execution
------------------

Another note for ordering:
We don’t actually have to display a column to sort by it.
For example, let’s say we want to order by the species ID, but we only want to see genus and species. 

    SELECT genus, species FROM species ORDER BY taxon ASC;

We can do this because sorting occurs earlier in the computational pipeline than field.

The computer is doing this:

1. Filtering rows according to WHERE
2. Sorting results according to ORDER BY
3. Displaying requested columns or expressions.

We can even sort based on the value of an expression. 
For example, sometimes ecologists want to look at samples or sites in a random order to avoid bias 
(i.e., sampling our favorite sites, or the ones that were most successful last time).
We can use the RANDOM function to return a pseudorandom integer for each row in a table. Watch it change:

    SELECT RANDOM() FROM plots
    
We can then generate a randomly ordered list of our plots for sampling.

    SELECT plot_id FROM plots ORDER BY RANDOM()


***Exercise: Let's try to combine what we've learned so far in a single query. 
Let’s go back to the surveys table and lets say that we want to display
the three date fields, species ID, and weight in kilograms (rounded to two 
decimal places),  for rodents captured in 1999, ordered alphabetically by 
the species ID.***

The order of the clauses is dictated by SQL: SELECT, FROM, WHERE, ORDER BY
and we often write each of them on their own line for readability.

**BREAK**

Aggregation
-----------
Aggregation allows us to combine results group records based on value
and calculate combined values in groups.

Let’s go to the surveys table and find out how many individuals there are.
Using the wildcard simply counts the number of records (rows)

    SELECT COUNT(*) FROM individuals

We can also find out how much all of those individuals weigh.

    SELECT SUM(wgt) FROM individuals

***Do you think you could output this value in kilograms or pounds,
rounded to 3 decimal places? Choose your favorite units.***

    SELECT ROUND(SUM(wgt)/1000.0, 3) FROM surveys

There are many other aggregate functions included in SQL including
MAX, MIN, and AVG.
 
***From the surveys table, can you use one query to output the total weight, average weight, and the min and max weights?***

Now, let's try to see how many individuals were counted in each species?

    SELECT species, COUNT(*)
    FROM surveys
    GROUP BY species

Why doesn't this work?

    SELECT COUNT(species)
    FROM surveys

***EXERCISE: Write queries that return:
1. How many individuals were counted in each year
2. Average weight of each species in each year***

We can order the results of our aggregation by a specific column, 
including the aggregated column.
Let’s count the number of individuals of each species captured,
ordered by the count

    SELECT species, COUNT(*)
    FROM surveys
    GROUP BY species
    ORDER BY COUNT(sp_code)

***Exercise: Write a query that lets us look at which years contained the most individuals and which had the least?***

***Exercise: Write a query that shows us which species had the largest individuals on average?
Which plots have the smallest individual on average?***

***Short break***

Database Design
----------------
Each field in a database should store a single value.
Information should not be duplicated in a database.
Each table should be about a single subject (avoids uneccesary replication).
When naming fields, think about meaning, not presentation.

When we divide our data between several tables, we need a way to bring it back together. 
The key is to have an identifier in common between tables - shared columns. 
This will allow us to JOIN tables.
This is what we will discuss now

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
Remember that you are smarter than your computer!

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

    SELECT plots.plot_type, ROUND(AVG(surveys.wgt),2)
    FROM surveys
    JOIN plots
    ON surveys.plot = plots.plot_id
    GROUP BY plots.plot_type

If you query starts to look messy, you can save yourself some typing by using aliases for the table names.

    SELECT p.plot_type, ROUND(AVG(s.wgt),2)
    FROM surveys s
    JOIN plots p
    ON s.plot = p.plot_id
    GROUP BY p.plot_type

Its also getting difficult to read some of our column names after running our queries. 
We can fix this too using AS:

    SELECT p.plot_type, ROUND(AVG(s.wgt),2) AS avg_wgt
    FROM surveys s
    JOIN plots p
    ON s.plot = p.plot_id
    GROUP BY p.plot_type
    
***Exercise: Find your query from earlier where you the total weight, average weight, and the min and max weights
and clean up your column names using concise descriptive names***


***Some's good, more's better***
We can join as many tables as we want using more JOIN clauses. 
Let's link to the species table so we can make sure that we are only including rodent species in our output.

    SELECT plots.plot_type, ROUND(AVG(surveys.wgt),2)
    FROM surveys
    JOIN plots
    JOIN species
    ON (surveys.plot = plots.plot_id) AND (surveys.species = species.species_id)
    WHERE taxa = Rodent
    GROUP BY plots.plot_type
    
***Exercise: At our site Onychomys (grasshopper mouse) is considered to be functionally different from our other species.
Can you make a query that joins the three tables and returns the data where the genus is not Onychomys?***
    
Primary Key
------------
The primary key is a value or combination of values that uniquely identifies each row in a table.
What is the primary key in each of our tables?

NULL
-----
Databases should use a special value for holes in the data.
Depending on your field, you may have seen many different things used -- What have you seen?
We recommend using a blank or NULL as best, because it can be read into many programs, and is unlikely to be mistaken as a 'real' value.
(e.g., it is common to use -999 in some fields, but this can be problematic!)

NULL is a special value. 
It does not equal 0, FALSE, or an empty string.
It means "no data here".
How might you return NULL data? Try:

    SELECT * FROM surveys WHERE wgt = NULL


    SELECT * FROM surveys WHERE wgt != NULL
    
This won't work because NULL cannot be compared to anything else. It is special. 
These statements will always be false and will never return any records.    
To find, we need to use a special operator, IS NULL or IS NOT NULL. 

    SELECT * FROM surveys WHERE wgt IS NULL


    SELECT * FROM surveys WHERE wgt IS NOT NULL
    
If there are NULL fields, we can sometimes run into problems with data analysis.
For example, there may be cases when we need to find missing data in our database, or remove rows with missing data.

***Exercise, return all the rows where species was not recorded, ordered by date***

Keep in mind that most aggregate functions skip NULL. 
Be aware of it and decide how you would like to proceed.
For example, if we are taking averages, the function will take the average of the non-null data and
divide it by the number of rows that are not null, rather than the actual number of rows.

***Break***
    
    
Adding data to existing tables
------------------------------
Sometimes we need to add information to an existing table, or add a new table.
We can also do these things right in SQL.

Let's say that for our study, we would like to start tracking weather conditions.
When creating a new table we need to name it, name the columns, and state the variable type for each column.

    CREATE TABLE weather (rowid integer primary key, day integer, month integer, year integer, 
    precip float, temp float, description text)
    
We can then update this new table by adding new information each time we collect more data.
Note: we don't need to update rowid because it is the primary key and it should autoincrement.

    INSERT INTO weather (day, month, year, precip, temp, description 
    VALUES 22, 06, 2013, 3.1, 25, 'sunny')


    INSERT INTO weather (day, month, year, precip, temp, description
    VALUES 23, 06, 2013, 0, 21, 'cloudy')
    
We can also change data in existing cells if we make a mistake, using UPDATE.

    UPDATE weather SET description = 'sunny' WHERE (day = 23) AND (month = 06) AND (year = 2013)
    
Finally, we can delete rows in tables:

    DELETE FROM weather WHERE (day = 23) AND (month = 06) AND (year = 2013)
    
***Exercise. We captured a new species at our site and have to update our species list.
It was a bird -- a canyon towhee (Melozone fusca). Give it a new species code and add it to the species table.
Also, the non-census rodent, Ammospermophilus harrisi should be spelled "harrisii" (Wilson and Reeder 2005).
Edit the table to correct it.





    



