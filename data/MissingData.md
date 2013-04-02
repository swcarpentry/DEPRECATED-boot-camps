## Missing data

Remember, when we added our stars named after ourselves and we only provided a `ProperName` value but no values for the other data in our database, for example a `Distance`. 

    >>> print stars.find_one({'ProperName':'YOURNAME'})

These documents can be viewed as having 'holes' or missing data. We can query for such documents,

    >>> print stars.find_one({'Distance':{'$exists':False}})

This returns the first document that does not have a `Distance` field.

MongoDB defines a special value, `null`, which, in pymongo, is represented as `None` (the conversion happens behind the scenes). This can be used to represent undefined data. So, let's remove our document,

    >>> stars.remove({'ProperName':'YOURNAME'})

This behaves in the same way as `find` but then removes all matching documents. Now, let's create a new version of our document,

    >>> doc = {}
    >>> doc['ProperName'] = 'YOURNAME'
    >>> doc['Group'] = 'EGI'
    >>> doc['Distance'] = None

We can look for documents with `null` values,

    >>> print stars.find_one({'Distance':None})

Note that while both versions of our document don't define a distance, in the original one we didn't define a distance at all, whereas in our new one, we have a distance whose value is undefined. We can see the difference as follows,

    >>> print stars.find_one({'Distance':{'$exists':False}})

This returns no document as it only succeeds if there is a document that has no `Distance` field.

The presence of `null` also affects the behaviour of operators. If we count the number of stars we have we get,

    >>> print docs.count()

If we now query the number of stars which have a distance greater than or equal to 0, we get,

    >>> print stars.find({'Distance':{'$gte':0}}).count()

Our documents with a `null` `Distance` are not counted.

It is up to us when we create our databases, or use them, to decide what it means for a field to be absent, or for a field to be present but with a `null` value. Likewise, it is up to us how we'll handle them. For example, when calculating the average distance, do we want to treat `null` or missing distances as 0, or not include documents containing such values in our total counts.

As we've seen MongoDB does allow both missing and `null` values to be checked for and handled explicitly.

Previous: [Counting and aggregating data](CountAggregate.md) Next: [Conclusions and further information](Conclusion.md)
