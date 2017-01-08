from Bio.Blast.Applications import NcbiblastpCommandline
import util.helper as helper
import util.blast as blast

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

    out_files = blast.blastp([p1, p2], "swissprot")

    for out_file in out_files:
        r = helper.read_blast(out_file)
        for a in r.alignments:
            print(a.hit_def)
            uniprot_id = extract_uniprot_id(a.hit_def)
            print(uniprot_id)

if __name__ == "__main__":
    main()
