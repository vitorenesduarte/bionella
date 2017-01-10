import util.www as www
import util.rw as rw
import util.blast as blast

def extract_uniprot_id(hit_def):
    """
    Extrair o uniprot id do 'hit_def' do alinhamento.
    NOTA: Hit def's que vêm do ncbi blast usando o biopython
    não contêm esta informação.
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

    #out_files = blast.blastp(proteins, "swissprot", "local")
    out_files = blast.blastp(["MRYTNIELLKRIPQHLKGVMEYYPDVLLFQLNQIQYTHLWHELSSAKYLYTNGFVIKPTHWLMYAFQTVKGWLGFDNHCQPEKISYVLDKLSYYGYTKQFLQPDFSSITNYSISPEISALVVKSRDDLTTAQLQNKLVNNYFKVEPHLGIDYSYQKLNPNHRFGESWMHAHEWGLIPQLDPQDDSLIAEVISTLDKKGISAKDVFFLQHSKYAKAAAQYYCDKAKNTIVPSFFFRLLWTDPRPRYLALALAYDPQIAKRDAQKFIEYHLAQKEYDNAFNLLGILNDSNLVLKFLLAIPETERHALIQKDTAIAAIMAKYYIEKKQYLLATQFYTNIEEISPNAAFAIEIQEQNYEKAYDIFKKYNSPNLFSMPERKLLAKIFYSDAETAYVAGKTYRGNKNWEKAKQYYLQSLEQKKAAHHLDPADEYLEDLYAHKRLYALLLIDADIDLNKAEDSDIASIQKAITLLRECRSSNKEEQEHHTRALAAGLMRRVDTLREKIAFNYYSSDFESIRKHKIEHQQDIAILIKTLEELIALLEGTQDKALRLQLGKAHYLLADVQCFFDINAPDINQHYKMAMKAVPGNPFYVLRVAELFEEEKSKLQKIGITQLKNMGYQVFDFLHWSEERWCKRDDIIHNIKDIHMPPSEPINSSSWNFTF"], "swissprot", "docker")

    for i in range(len(proteins)):
        lt = tags[i]
        out_file = out_files[i]
        print("TAG:", lt)

        ids = []

        handle = rw.read_blast(out_file)
        for a in handle.alignments:
            uniprot_id = extract_uniprot_id(a.hit_def)
            ids.append(uniprot_id)

        print(ids)
        print(www.fetch_uniprots(ids))

if __name__ == "__main__":
    main()
