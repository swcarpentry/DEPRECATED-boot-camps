## If there is time...

### Create a database and collection

    >>> database = connection['YOURNAME-db']
    >>> examples = database['examples']

### Insert a document

    >>> doc = {}
    >>> doc['name'] = 'YOURNAME'
    >>> doc['age'] = YOURAGE
    >>> doc['organisation'] = None
    >>> examples.insert(doc)
    >>> print examples.find_one()

### Update a document

Update an existing field,

    >>> examples.update({'name':'YOURNAME'}, {'$set': {'organisation':'YOURORGANISATION'}})
    >>> print examples.find_one()

The first argument is a query. If more than one document matches the query then all the matching documents are updated.

Add a new field,

    >>> examples.update({'name':'YOURNAME'}, {'$set': {'position':'YOURPOSITION'}})
    >>> print examples.find_one()

### Remove a document

    >>> examples.remove({'name':'YOURNAME'})
    >>> print examples.find_one()

The first argument is a query. If more than one document matches the query then all the matching documents are removed. When they're gone, they're gone for good, so beware!

### Delete a collection

    >>> database.drop_collection('examples')
    >>> print database.collection_names()

When it's gone, it's gone for good, so beware!

### Delete a database

    >>> connection.drop_database('YOURNAME-db')
    >>> print connection.database_names()

When it's gone, it's gone for good, so beware!

### Map-Reduce

http://docs.mongodb.org/manual/applications/map-reduce/

* Distance - parsecs
* Mag - apparent visual magnitude
* AbsMag - absolute visual magnitude (its apparent magnitude from a distance of 10 parsecs)

* Sum: 42501640892.5
* Count: 119617
* Average distance: 355314.385853

Example,

    >>> from bson.code import Code
    >>> map = Code("function () { emit('totals', {count: 1, distance: this.Distance})}")
    >>> reduce = Code("function (key, values) { var totals = {count: 0, distance: 0}; values.forEach(function(value) { totals.count += value.count; totals.distance += value.distance;}); return totals;}")
    >>> finalize = Code("function(key, values) { return {count: values.count, distance: values.distance, avgdistance: values.distance / values.count}}")
    >>> stars.map_reduce(map, reduce, {'inline':1}, query={'ProperName':'Sol'}, finalize=finalize)
    >>> stars.map_reduce(map, reduce, {'inline':1}, finalize=finalize)

Replace `inline` with `result` to create new a `result` collection.

    >>> map = Code("function () { emit('totals', {count: 1, distance: this.Distance, mag: this.Mag, absmag: this.AbsMag})}")
    >>> reduce = Code("function (key, values) { var totals = {count: 0, distance: 0, mag: 0, absmag: 0}; values.forEach(function(value) { totals.count += value.count; totals.distance += value.distance; totals.mag += value.mag; totals.absmag += value.absmag;}); return totals;}")
    >>> finalize = Code("function(key, values) { return {count: values.count, distance: values.distance, mag: values.mag, absmag: values.absmag, avgdistance: values.distance / values.count, avgmag: values.mag / values.count, avgabsmag: values.absmag / values.count}}")
    >>> stars.map_reduce(map, reduce, {'inline':1}, query={'ProperName':'Sol'}, finalize=finalize)
    >>> stars.map_reduce(map, reduce, {'inline':1}, finalize=finalize)

Previous: [Conclusions and further information](Conclusion.md)
