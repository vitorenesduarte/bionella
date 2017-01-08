#!/bin/bash
if [ -z "$QUERY_DIR" ]; then
  echo ">>> QUERY_DIR is not configured; please export it."
  exit 1
fi

cd $QUERY_DIR

for QUERY in $(ls)
do
  echo "Running blastp for query $QUERY"
  $BLAST_BIN/blastp -query $QUERY -db swissprot -outfmt 5 -out $QUERY.xml
done

