import util.www as www
import util.rw as rw
import util.blast as blast

def add_info_to_dictionary(dictionary, uniprots):
    """
    Esta função adiciona ao dicionário com informação do NCBI
    a informação retirada do site da uniprot.
      - short_name
      - product
      - EC_number
      - comment_functions
      - cofactors
      - pathologies
      - pdbs
      - molecular_functions
      - biological_processes
      - locations
      - translation
      - length
      - mass
    """

    # propriedades que já tinhamos:
    old_properties = [
        "EC_number",
        "translation",
        "length"
    ]

    new_properties = [
        "short_name",
        "product",
        "comment_functions",
        "cofactors",
        "pdbs",
        "pathologies",
        "molecular_functions",
        "biological_processes",
        "locations",
        "mass"
    ]

    diff = {}

    for tag in sorted(dictionary.keys()):
        if tag in uniprots:
            for p in old_properties:
                new = uniprots[tag][p]
                if p in dictionary[tag]:
                    # se já existe, verifica que é o mesmo
                    current = dictionary[tag][p]
                    if not new == current:
                        # guarda informação sobre o que diferiu

                        if tag in diff:
                            diff[tag].append(p)
                        else:
                            diff[tag] = [p]
                
                # substitui sempre pelo valor da uniprot
                dictionary[tag][p] = new

            for p in new_properties:
                dictionary[tag][p] = uniprots[tag][p]

    return (diff, dictionary)

def tag_as_key(dictionary, uniprots):
    """
    Retorna um dicionário com os resultados da uniprot,
    em que as chaves são locus tag em vez de uniprot_id.
    """
    d = {}

    for tag in dictionary:
        if "uniprot_id" in dictionary[tag]:
            uniprot_id = dictionary[tag]["uniprot_id"]

            info = uniprots[uniprot_id]
            d[tag] = info

    return d

def main():
    ncbi_json_path = ".ncbi.json"
    uniprots_json_path = ".uniprots.json"
    ncbi_uniprot_json_path = ".ncbi_uniprot.json"
    ncbi_uniprot_diff_json_path = ".ncbi_uniprot_diff.json"

    dictionary = rw.read_json(ncbi_json_path)

    uniprot_ids = set()
    for tag in dictionary:
        if "uniprot_id" in dictionary[tag]:
            uniprot_id = dictionary[tag]["uniprot_id"]
            uniprot_ids.add(uniprot_id)

    uniprot_ids = list(uniprot_ids)

    uniprots = www.fetch_uniprots(uniprot_ids)
    uniprots = tag_as_key(dictionary, uniprots)
    rw.write_json(uniprots, uniprots_json_path)

    uniprots = rw.read_json(uniprots_json_path)

    (diff, dictionary) = add_info_to_dictionary(
        dictionary,
        uniprots
    )
    rw.write_json(diff, ncbi_uniprot_diff_json_path)
    rw.write_json(dictionary, ncbi_uniprot_json_path)
    

if __name__ == "__main__":
    main()
