import util.www as www
import util.rw as rw
import util.blast as blast

def extract_uniprot_id(hit_def):
    """
    Extrair o uniprot id do 'hit_def' do alinhamento.
    """
    parts = hit_def.split("|")
    return parts[1]

def main():
    ncbi_json_path = ".ncbi.json"
    dictionary = rw.read_json(ncbi_json_path)

    tags = []
    proteins = []

    for lt in dictionary:
        if dictionary[lt]["type"] == "mRNA":
            # se for do tipo "mRNA" estamos na presença de uma proteína
            protein = dictionary[lt]["translation"]

            tags.append(lt)
            proteins.append(protein)


    #proteins = proteins[40:50]

    out_files = blast.blastp(proteins, "swissprot", "local")

    for i in range(len(proteins)):
        lt = tags[i]
        out_file = out_files[i]
        print("TAG:", lt)

        ids = []

        handle = rw.read_blast(out_file)
        for a in handle.alignments:
            uniprot_id = extract_uniprot_id(a.hit_def)
            ids.append(uniprot_id)

        m = www.fetch_uniprots(ids)
        rw.write_json(m, lt + ".json")

if __name__ == "__main__":
    main()
