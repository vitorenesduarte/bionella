from Bio.Blast import NCBIWWW
import helper

def blastp(protein, file_path):
    """
    Corre o blast para a proteína passada como argumento
    e grava o resultado em file_path
    """
    handle = NCBIWWW.qblast("blastp", "nr", protein)
    helper.write_blast(handle, file_path)

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

    blastp("CCC", p1_f)
    #blastp(p2, p2_f)
    p1_r = helper.read_blast(p1_f)
    print(p1_r)

if __name__ == "__main__":
    main()
