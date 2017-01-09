from Bio import Entrez
from Bio import SeqIO
import util.www as www
import util.rw as rw

def extract_features(record):
    """
    Das várias features encontradas no ficheiro genbank,
    decidimos extrair:
      - gene
      - CDS
      - tRNA
      - rRNA

    Esta função retorna uma lista apenas com as features
    dos tipos indicados acima.
    """

    result = []
    features = record.features
    feature_count = len(features)
    types = ["gene", "CDS", "tRNA", "rRNA"]

    for i in range(feature_count):
        feature = features[i]
        should_extract = feature.type in types

        if should_extract:
            result.append(feature)

    return result

def features_to_dictionary(features):
    """
    A cada uma das features, extraimos a seguintes propriedades:
      - db_xref
      - EC_number
      - function
      - gene
      - gene_synonym
      - locus_tag
      - note
      - product
      - protein_id
      - translation
    """
    dictionary = {}
    properties = ["db_xref",
                  "EC_number",
                  "function",
                  "gene",
                  "gene_synonym",
                  "note",
                  "product",
                  "protein_id",
                  "translation"]

    for feature in features:
        tag = feature.qualifiers["locus_tag"][0]

        if not tag in dictionary:
            # se esta tag ainda não existe,
            # criar um dicionário vazio para ela
            dictionary[tag] = {}

        if feature.type in ["tRNA", "rRNA"]:
            # Taggar as features com estes dois tipos
            dictionary[tag]["type"] = feature.type
        else:
            # se não, taggar com "mRNA"
            dictionary[tag]["type"] = "mRNA"

        for prop in properties:
            if prop in feature.qualifiers:
                value = feature.qualifiers[prop][0]

                if prop in dictionary[tag]:
                    # se já encontramos esta propriedade
                    # para este locus_tag
                    # verificar que é a mesma
                    current_value = dictionary[tag][prop]
                    assert value == current_value
                else:
                    # se não, adicionar
                    dictionary[tag][prop] = value

    return dictionary

def main():
    start = 270001
    end = 505535
    ncbi_gb_path = ".ncbi.gb"
    ncbi_json_path = ".ncbi.json"
    table_json_path = ".table.json"

    ## 1.1
    record = www.fetch_genbank(start, end)
    rw.write_genbank(record, ncbi_gb_path)

    record = rw.read_genbank(ncbi_gb_path)
    features = extract_features(record)

    dictionary = features_to_dictionary(features)
    rw.write_json(dictionary, ncbi_json_path)

    ## 1.2
    table = www.fetch_table()
    rw.write_json(table, table_json_path)

    dictionary = rw.read_json(ncbi_json_path)
    table = rw.read_json(table_json_path)
    ## TODO
    ## verificar que o que está em "table"
    ## também esta em dictionary
    ## se não estiver, adicionar
    ## se estiver e for diferente,
    ## guardar ambos

if __name__ == "__main__":
    main()
