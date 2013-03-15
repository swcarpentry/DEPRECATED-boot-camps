## Connect to our database server

First we will connect to our database server which is running on our local server on MongoDB's default port,

    $ python
    >>> from pymongo import MongoClient
    >>> connection = MongoClient()

We can now see the databases on the server,

    >>> print connection.database_names()

## Browse databases and collections

To use the `astronomy` database we do,

    >>> astro = connection['astronomy']

This returns us an object, or handle, that represents our database. If the database of that name does not exist, MongoDB would create a new, empty, database for us.

We can see the collections in `astronomy` by doing,

    >>> print astro.collection_names()

The first is a *system collection*, which are for MongoDB's internal use. We're interested in `stars` so, to use this, we do,

    >>> stars = astro['stars']

This returns us an object, or handle, that represents our collection. If one did not exist, MongoDB would create a new, empty, collection for us.

## Looking at a document

Let's get a document from our collection,

    >>> doc = stars.find_one()
    >>> print doc

This gets a document from the collection then prints it. This is a BSON document and can be used in the same way as a Python dictionary,

    >>> print doc['ProperName']

Note the `_id` field,

    >>> print doc['_id']

This is the unique ID of the document in the MongoDB database. When we insert documents into the database, this is automatically added.

## Selecting documents

`find_one` returns us a single document, but we want more control over what we return. We can issue a *query* to MongoDB telling us what we want to give us. We leave MongoDB to worry about how to get it.

In running `find_one` above we have issued a query, albeit the most general, indifferent one we can, which says "find me a document".

Let's say we want to find a star by name. We can do this as follows,

    >>> print stars.find_one({'ProperName':'Sol'})

This says "find me a document which has a `ProperName` field whose value is `Sol`". The query is expressed as a Python dictionary. Let's try some more,

    >>> print stars.find_one({'ProperName':'Betelgeuse'})
    >>> print stars.find_one({'ProperName':'Aldebaran'})
    >>> print stars.find_one({'ProperName':'Proxima Centauri'})

These queries return the whole document. We may just be interested in parts of each document, so we can provide a list specifying the fields we're interested in,

    >>> fields = ['ProperName', 'Distance']
    >>> print stars.find_one({'ProperName':'Betelgeuse'}, fields)

This is termed a *projection*. Since BSON documents can be nested - a field's value is a document itself - MongDB also allows the sub-documents to be projected.

## Adding a document

Let's add our own star named after ourselves. First we can create a BSON document,

    >>> doc = {}
    >>> doc['ProperName'] = 'YOURNAME'
    >>> doc['Group'] = 'EGI'

Though our `stars` documents have many fields, because MongoDB is schema-less, we don't have to worry about providing values for them all. We can also add additional fields, like `Group` if we want.

Now, we can add it to our collection, and we'll get it's ID,

    >>> id = stars.insert(doc)
    >>> print id

Now we can use the ID to retrieve our document,

    >>> print stars.find_one({'_id':id})

And we should see our document, which is now in the database. And, let's run another query,

    >>> print stars.find_one({'ProperName':'YOURNAME'})

Now, let's try,

    >>> docs = stars.find({'Group':'EGI'})

Would anyone like to guess what this will do?

    >>> print docs

`docs` holds a reference to a *cursor* which holds references to the matching documents. We can then iterate through this cursor,

    >>> for doc in docs: print doc

We can only iterate through the cursor once.

    >>> for doc in docs: print doc

To do so again, we need to *rewind*,

    >>> docs.rewind()
    >>> for doc in docs: print doc

## Back to selecting...

Let's use `find` to get the names and distances of all stars,

    >>> print fields
    >>> docs = stars.find({}, fields, limit = 10)
    >>> for doc in docs: print doc

Our query is `{}` which means get all the documents, and we've added a `limit` argument which limits us to the first 10 documents matching our query.

Suppose we don't want to see the `_id`. We can exclude it as follows,

    >>> fields = {'ProperName':True, 'Distance':True, '_id':False}
    >>> docs = stars.find({}, fields, limit = 10)
    >>> for doc in docs: print doc

We now use a dictionary rather than a list to specify the fields we want, and we explicitly state we want `_id` not to be included.

Note how many the `ProperName`s are `''`. Our data is incomplete, not all stars have been named.

