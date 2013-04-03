## Filtering and sorting data

One of the most powerful features of a database is the ability to filter through our data to find that which matches certain criteria. We've done that above, by simple matching on fields.

Let's now do something different, and look for all stars that are 1 parsec or closer,

    >>> docs = stars.find({'Distance':{ '$lt': 1}})
    >>> for doc in docs: print doc

Just one, our Sun! How about 2 parsecs,

    >>> docs = stars.find({'Distance':{ '$lt': 2}})
    >>> for doc in docs: print doc

Let's just return the `ProperName` and `Distance` again,

    >>> docs = stars.find({'Distance':{ '$lt': 4}}, fields)
    >>> for doc in docs: print doc

To find the names and distances of stars between 4 and 5 parsecs we can add a conditional operator,

    >>> docs = stars.find({'$and': [{'Distance':{'$gt': 4}}, {'Distance':{ '$lt': 5}}]}, fields)
    >>> for doc in docs: print doc

If we want to just see the how many distinct names are in our results, we can do,

    >>> docs.distinct('ProperName')

We can use this across the collection as a whole,

    >>> stars.distinct('ProperName')

As another example of filtering, let's return all the documents whose `ProperName` is in a list of those we're interested in,

    >>> names = ['Sol', 'Aldebaran', 'Betelgeuse']
    >>> query = {'ProperName': { '$in': names}}
    >>> docs = stars.find(query, fields)
    >>> for doc in docs: print doc

What MongoDB does here is checks that that the value of the `ProperName` field is in the given list (the `'$in': names` bit), and then just returns us the `ProperName` and `Distance`.

MongoDB has a rich set of [operators](http://docs.mongodb.org/manual/reference/operators/) that allow searches based on comparisons, logical operators, list inclusion and exclusion and geo-spatial operations based on distances, proximity or bounding boxes.

## Sorting

We can also request the data be sorted. So let's query for all the documents, pull out the name and distance, then sort by name. But we'll ignore stars that have no name, so,

    >>> from pymongo import ASCENDING
    >>> docs = stars.find({'ProperName':{'$ne':''}}, fields, sort=ASCENDING)

And to sort by name in ascending order, we do,

    >>> docs = stars.find({'ProperName':{'$ne':''}}, fields, sort=[('ProperName', ASCENDING)])
    >>> for doc in docs: print doc

And for descending,

    >>> from pymongo import DESCENDING
    >>> docs = stars.find({'ProperName':{'$ne':''}}, fields, sort=[('ProperName', DESCENDING)])
    >>> for doc in docs: print doc

Note how `sort` specifies a list of pairs. This allows us to determine how to sort by a secondary field should some of the values of a primary field all be equal,

    >>> docs = stars.find({'ProperName':{'$ne':''}}, fields, sort=[('Distance', ASCENDING), ('ProperName', ASCENDING)])

Previous: [Connecting to our database, browsing and selecting documents](Select.md) Next: [Counting and aggregating data](CountAggregate.md)
