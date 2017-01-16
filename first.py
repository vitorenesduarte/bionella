from Bio import Entrez
from Bio import SeqIO
import util.www as www
import util.rw as rw

def extract_features(record):
    """
    Das várias features encontradas no ficheiro genbank, extraímos:
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

def features_to_dictionary(start, features):
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

    Também extraímos a localização:
      - start
      - end
      - strand
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

        location = feature.location
        dictionary[tag]["location"] = {}
        dictionary[tag]["location"]["start"] = location._start + start
        dictionary[tag]["location"]["end"] = location._end + start - 1
        dictionary[tag]["location"]["strand"] = location._strand

        if feature.type in ["tRNA", "rRNA"]:
            # Taggar as features com estes dois tipos
            dictionary[tag]["type"] = feature.type
        else:
            # se não, taggar com "mRNA"
            dictionary[tag]["type"] = "mRNA"

        for prop in properties:
            if prop in feature.qualifiers:
                value = feature.qualifiers[prop][0]

                if prop == "db_xref":
                    # renomear esta propriedade e remover "GeneID:"
                    prop = "gene_id"
                    value = value[7:]

                if prop == "translation":
                    # se translation, também guardar o tamanho
                    dictionary[tag]["length"] = len(value)

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

def verify(dictionary, table):
    """
    Verifica se os valores presentes na tabela correspondem
    aos valores presentes no dicionário.
    Se encontrar alguma diferença retorna um dicionário
    para análise manual dos resultados.
    """
    diff = {}

    for tag in dictionary:
        # para cada uma das locus_tag da nossa zona
        if tag in table:
            # se estiver na tabela
            values = table[tag]
            for prop in values:

                dict_value = dictionary[tag][prop]
                table_value = values[prop]

                if not table_value == dict_value:
                    diff[tag] = {}
                    diff[tag]["prop"] = prop
                    diff[tag]["dictionary"] = dict_value
                    diff[tag]["table"] = table_value

    return diff

def retrieve_uniprot_ids(dictionary):
    """
    Adicionar ao dictionário uma nova propriedade:
    - uniprot_id
    """
    gene_id_to_tag = {}

    for tag in dictionary:
        if "gene_id" in dictionary[tag]:
            gene_id = dictionary[tag]["gene_id"]
            gene_id_to_tag[gene_id] = tag

    gene_ids = list(gene_id_to_tag.keys())
    gene_id_to_uniprot_id = www.gene_ids_to_uniprot_ids(gene_ids)

    for gene_id in gene_id_to_uniprot_id:
        uniprot_id = gene_id_to_uniprot_id[gene_id]
        tag = gene_id_to_tag[gene_id]

        # adicionary uma nova propriedade ao dicionário: "uniprot_id"
        dictionary[tag]["uniprot_id"] = uniprot_id

    return dictionary

def main():
    start = 270001
    end = 505535
    ncbi_gb_path = ".ncbi.gb"
    ncbi_json_path = ".ncbi.json"
    table_json_path = ".table.json"

    ## 1.1
    #record = www.fetch_genbank(start, end)
    #rw.write_genbank(record, ncbi_gb_path)

    record = rw.read_genbank(ncbi_gb_path)
    features = extract_features(record)

    dictionary = features_to_dictionary(start, features)
    rw.write_json(dictionary, ncbi_json_path)

    ## 1.2
    #table = www.fetch_table()
    #rw.write_json(table, table_json_path)

    dictionary = rw.read_json(ncbi_json_path)
    #table = rw.read_json(table_json_path)
    #diff = verify(dictionary, table)
    ### análise manual de diff

    dictionary = retrieve_uniprot_ids(dictionary)
    rw.write_json(dictionary, ncbi_json_path)

if __name__ == "__main__":
    main()
