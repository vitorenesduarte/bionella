#!/bin/bash

ENV_VARS=(
  QUERY_DIR
  DB
)

# Verificar que todas as variáveis de ambiente necessárias estão configuradas.
for ENV_VAR in "${ENV_VARS[@]}"
do
  if [ -z "${!ENV_VAR}" ]; then
    echo ">>> ${ENV_VAR} is not configured; please export it."
    exit 1
  fi
done

# Correr o blast.
cd $QUERY_DIR
for IN_FILE in $(ls)
do
  OUT_FILE=$IN_FILE".xml"
  $BLAST_BIN/blastp -query $IN_FILE -db $DB -outfmt 5 -out $OUT_FILE
done
