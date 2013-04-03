## Counting and aggregating data

Let's suppose we're interested in the number of stars within 10 parsecs. We'd run

    >>> docs = stars.find({'Distance':{ '$lt': 10}})

Rather than iterate through our cursor and count the number of documents, we can ask MongoDB to do that for us,

    >>> print docs.count()

This is an example of an *aggregation* operator, since it summarises data from multiple documents.

MongoDB has an [aggregation framework](http://docs.mongodb.org/manual/applications/aggregation/) to support such operations across documents. So, for example, we may want to sum the distances of all the stars. To do this we first write,

    >>> result = stars.aggregate([{'$group':{'_id':None, 'total':{'$sum':'$Distance'}}}])
    >>> print result

What this does is,

* Group the documents according to the value of the `_id` field above. In this case the value of the field is not defined (`None`) so all the documents are bundled into a single group.
* For each group...
 * Create a new document, with an `_id` equal to the group ID and a `total` field with value 0.
 * For each document in the group, get the value of its `Distance` field and sum this to the `total` field.
* Return a list of the documents for each group. Again, for us there is only 1 group so the list has one element.

We can then extract our total by...

    >>> total = result['result'][0]['total']
    >>> print total

We can also use this as another way of calculating the number of documents in the collection, by replacing `'$Distance'` with `1`, so that when each document is processed, 1 is added to the `count` field.

    >>> result = stars.aggregate([{'$group':{'_id':None, 'count':{'$sum':1}}}])
    >>> print result
    >>> count = result['result'][0]['count']
    >>> print count

And just to check,

    >>> print stars.count()

And now we can calculate our average,

    >>> print total / count

We can sum and count in one go,

    >>> result = stars.aggregate([{'$group':{'_id':None, 'total':{'$sum':'$Distance'}, 'count':{'$sum':1}}}])
    >>> print result
    >>> total = result['result'][0]['total']
    >>> count = result['result'][0]['count']

We saw when selecting documents that we could add conditions to select documents according to various criteria e.g. the number of stars less than a parsec away. We can do this in the aggregation framework too, for example, to count the number of stars less than a parsec away and the total distance of these,

    >>> stars.aggregate([{'$match': {'Distance':{ '$lt': 1}}}, {'$group':{'_id':None, 'total':{'$sum':'$Distance'}, 'count':{'$sum':1}}}])

This is just our Sun, which we can confirm by doing,

    >>> print stars.find_one({'Distance':{ '$lt': 1}}, ['ProperName', 'Distance'])

And for those less than 2 parsecs,
    
    >>> stars.aggregate([{'$match': {'Distance':{ '$lt': 2}}}, {'$group':{'_id':None, 'total':{'$sum':'$Distance'}, 'count':{'$sum':1}}}])

And, we can confirm the count,

    >>> docs = stars.find({'Distance':{ '$lt': 2}}, ['ProperName', 'Distance']).count()
    >>> print docs.count()

This is just a dip in the water of the capabilities of MongoDB's aggregation framework.

Rather than pulling all the data out of the database and writing a program to process it and calculate these values, instead we are driving our computation to the database, telling it what we want and letting it worry about how to do it.

Previous: [Filtering and sorting data](FilterSort.md) Next: [Missing data](MissingData.md)
