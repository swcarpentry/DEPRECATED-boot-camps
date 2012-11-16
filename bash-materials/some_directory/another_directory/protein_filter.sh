#!/bin/bash

DATE=`date`
ATOM="H"
PDB_EXT="pdb"

echo "Processing date: $DATE"
for PDB in `find . -name "*.$PDB_EXT"`; do
    COUNT=`grep -w $ATOM $PDB | wc -l`
    echo "$PDB $COUNT"
done
echo "Processing completed!"
