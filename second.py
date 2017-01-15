import util.www as www
import util.rw as rw
import util.blast as blast

def main():
    ncbi_json_path = ".ncbi.json"
    blast_results_json_path = ".blast_results.json"

    dictionary = rw.read_json(ncbi_json_path)
    tags_and_proteins = []

    for tag in dictionary:
        if dictionary[tag]["type"] == "mRNA":
            # se for do tipo "mRNA" estamos na presença de uma proteína
            protein = dictionary[tag]["translation"]
            both = (tag, protein)
            tags_and_proteins.append(both)

    blast_results = blast.blastp(
        tags_and_proteins,
        "swissprot",
        "local"
    )
    rw.write_json(blast_results, blast_results_json_path)

if __name__ == "__main__":
    main()
