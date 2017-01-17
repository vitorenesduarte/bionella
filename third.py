import util.www as www
import util.rw as rw
import util.blast as blast

def add_info_to_dictionary(dictionary, uniprots):
    """
    Esta função adiciona ao dicionário com informação do NCBI
    a informação retirada do site da uniprot.
    """

    return dictionary

def main():
    ncbi_json_path = ".ncbi.json"
    ncbi_uniprot_json_path = ".ncbi_uniprot.json"

    dictionary = rw.read_json(ncbi_json_path)

    uniprot_ids = set()
    for tag in dictionary:
        uniprot_id = dictionary["tag"]["uniprot_id"]
        uniprot_ids.add(uniprot_id)

    uniprot_ids = list(uniprot_ids)

    uniprots = www.fetch_uniprots(uniprot_ids)
    rw.write_json(uniprots, uniprots_json_path)

    uniprots = rw.read_json(uniprots_json_path)

    dictionary = add_info_to_dictionary(
        dictionary,
        uniprots
    )
    rw.write_json(dictionary, ncbi_uniprot_json_path)
    

if __name__ == "__main__":
    main()
