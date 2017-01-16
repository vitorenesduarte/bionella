import util.www as www
import util.rw as rw
import util.blast as blast

def get_tags_and_proteins(dictionary):
    """
    Esta função retorna uma lista de pares
    em que a primeira componente é a tag
    e a segunda a sequência de aminoácidos.
    """

    tags_and_proteins = []

    for tag in dictionary:
        if dictionary[tag]["type"] == "mRNA":
            # se for do tipo "mRNA" estamos na presença de uma proteína
            protein = dictionary[tag]["translation"]
            both = (tag, protein)
            tags_and_proteins.append(both)

    return tags_and_proteins

def add_info_to_blast_results(blast_results, uniprots):
    """
    Esta função adiciona aos blast results informação
    retirada do site da uniprot.
    """

    for tag in blast_results:
        for i in range(len(blast_results[tag])):
            uniprot_id = blast_results[tag][i]["uniprot_id"]

            properties = [
                "molecular_functions",
                "comment_functions",
                "sequence",
                "status"
            ]

            for p in properties:
                value = uniprots[uniprot_id][p]
                blast_results[tag][i][p] = value

    return blast_results

def main():
    ncbi_json_path = ".ncbi.json"
    blast_results_json_path = ".blast_results.json"
    uniprots_json_path = ".uniprots.json"

    dictionary = rw.read_json(ncbi_json_path)
    tags_and_proteins = get_tags_and_proteins(dictionary)

    blast_results = blast.blastp(
        tags_and_proteins,
        "swissprot",
        "expasy"
    )
    rw.write_json(blast_results, blast_results_json_path)

    blast_results = rw.read_json(blast_results_json_path)

    uniprot_ids = set()
    for tag in blast_results:
        for result in blast_results[tag]:
            uniprot_ids.add(result["uniprot_id"])

    uniprot_ids = list(uniprot_ids)

    uniprots = www.fetch_uniprots(uniprot_ids)
    rw.write_json(uniprots, uniprots_json_path)

    uniprots = rw.read_json(uniprots_json_path)

    blast_results = add_info_to_blast_results(
        blast_results,
        uniprots
    )
    rw.write_json(blast_results, blast_results_json_path)
    

if __name__ == "__main__":
    main()
