import util.helper as helper
import util.blast as blast

def extract_uniprot_id(hit_def):
    parts = hit_def.split("|")
    return parts[1]

def main():
    ncbi_json_path = ".ncbi.json"

    dictionary = helper.read_json(ncbi_json_path)
    for lt in dictionary:
        print("LT:", lt, " | type:", dictionary[lt]["type"])

        if dictionary[lt]["type"] == "mRNA":
            # se for do tipo "mRNA" estamos na presença de uma proteína
            protein = dictionary[lt]["translation"]
            [out_file] = blast.blastp([protein], "swissprot", "docker")

            handle = helper.read_blast(out_file)
            for a in handle.alignments:
                uniprot_id = extract_uniprot_id(a.hit_def)
                print("UI:", uniprot_id)

if __name__ == "__main__":
    main()
