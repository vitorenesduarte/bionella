from Bio.Blast import NCBIWWW
from Bio.Blast.Applications import NcbiblastpCommandline
import helper

def blastp(protein, file_path):
    """
    Corre o blast para a proteína passada como argumento
    e grava o resultado em file_path
    """
    handle = NCBIWWW.qblast("blastp", "nr", protein)
    helper.write_blast(handle, file_path)

def local_blastp(protein, file_path):
    """
    Corre o blast localmente para a proteína passada como argumento
    e grava o resultado em file_path.
    Assume a existência de uma base de dados local chamada
    'swissprot'.
    """

    # Gravar a proteína num ficheiro
    query_path = ".query"
    helper.write_file(protein, query_path)

    # Correr o blast
    blastp_cline = NcbiblastpCommandline(
        query=query_path,
        db="swissprot",
        evalue=10,
        outfmt=5,
        out=file_path)
    blastp_cline()

def extract_uniprot_id(hit_def):
    parts = hit_def.split("|")
    return parts[1]

def main():
    ncbi_json_path = ".ncbi.json"

    dictionary = helper.read_json(ncbi_json_path)
    for lt in dictionary:
        if dictionary[lt]["type"] == "mRNA":
            # se for do tipo "mRNA" estamos na presença de uma proteína
            protein = dictionary[lt]["translation"]
            # ...

    p1 = "MIGCCLIIFPCNRDYDKIVRTNYSLVRRLMKHNLLDKAYKHCVNHGYRFTEPRERVLKILVDERKPLGAYDILQRLSMEVDNPKPPTVYRAIQFWHQEGFIHCIDSLKSYVACLHGHHVGQAQFLICNQCDFVKELECVVDFTAVTEAANSIQFSIINCTVEIKGLCSDCNLTNLKN"
    p2 = "MKKAFRIMATYYGQCACGKIQFMCTGEPVFTQYCHCNKCREIASLSQKKQDKSGYSLTAAYLTRDFRILSEDNDFEEIIKEKAKLFLCSYCHSLVYGIALDPAQQDSIGINVNNFSFDTSIPDSFKPVRHIWYASRIMDCDDQLPKFKDAPKEQFGSGELFELPDAN"

    p1_f = ".p1.xml"
    p2_f = ".p2.xml"

    #print("LOCAL")
    #local_blastp(p1, p1_f)
    #r = helper.read_blast(p1_f)
    #for a in r.alignments:
    #    uniprot_id = extract_uniprot_id(a.hit_def)
    #    print(uniprot_id)

    #print("\n\n\nWEB")
    #blastp(p1, p1_f)
    r = helper.read_blast(p1_f)
    for a in r.alignments:
        #print(a.hit_def)
        #uniprot_id = extract_uniprot_id(a.hit_def)
        #print(uniprot_id)


        hsp = a.hsps[0]
        print("- e:", hsp.expect)
        print("- score:", hsp.score)
        print("- align len:", hsp.align_length)
        print("- matches:", hsp.identities)
        print("- accession:", a.accession)
        print("*", a.title)
        print("*", a.hit_id)
        print("*", a.length)
        print("*", a.hit_def)

if __name__ == "__main__":
    main()
