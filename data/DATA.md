# Sample data

The example data used in this session can be created using 
the HYG star database archive,

https://github.com/astronexus

For more background, see http://astronexus.com/node/34.

Assuming you have a MongoDB database listening on its 
standard/default port, you can populate it as follows,

Download the AstroNexus HYG database,

    $ git clone https://github.com/astronexus/HYG-Database.git

Run the Python script, specifying the path to `hygxyz.csv`,
a database name and a collection name,

    $ python csv_to_json.py HYG-Database/hygxyz.csv astronomy stars

See `csv_to_json.py` for more information.
