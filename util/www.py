from Bio import Entrez, SeqIO
from urllib.request import urlretrieve
from xml.etree import ElementTree
import util.rw as rw

def fetch_xml_tree(url, add_root=False):
    """
    Faz download de um ficheiro xml, faz parsing dele,
    e retorna um 'ElementTree'.
    """
    local_file, _ = urlretrieve(link)

    if add_root:
        rw.wrap_file("<root>", "</root>", local_file)

    fd = open(local_file, "r")
    tree = ElementTree.parse(fd)
    fd.close()

    return tree

def fetch_genbank(start, end):
    """
    Procura no NCBI pelo Accession NC_002942.5 e 
    faz download de toda a informação em formato genbank,
    retornando o respectivo record com as features
    entre a posição start e end.
    """

    Entrez.email = "vitorenesduarte@gmail.com"

    handle = Entrez.esearch(db="nucleotide", term="NC_002942.5")
    record = Entrez.read(handle)
    genbank_id = record["IdList"][0]
    handle.close()

    handle = Entrez.efetch(db="nucleotide",
                           rettype="gbwithparts",
                           retmod="text",
                           id=genbank_id,
                           seq_start=start,
                           seq_stop=end)
    record = SeqIO.read(handle, "gb")
    handle.close()

    return record

def fetch_table():
    """
    A página https://www.ncbi.nlm.nih.gov/genome/proteins/416?genome_assembly_id=166758
    faz um pedido http ao url_prefix para preencher a tabela.
    Este método implementa a paginação que seria feita manualmente
    no website, e guarda a informação que queremos num dicionário.

      - A coluna 5  corresponde à propriedade db_xref    do genbank
      - A coluna 6  corresponde à propriedade gene       do genbank
      - A coluna 7  corresponde à propriedade locus_tag  do genbank
      - A coluna 8  corresponde à propriedade protein_id do genbank
      - A coluna 11 corresponde à propriedade product    do genbank
    """

    url_prefix = "https://www.ncbi.nlm.nih.gov/genomes/Genome2BE/genome2srv.cgi?action=GetFeatures4Grid&type=Proteins&genome_id=416&genome_assembly_id=166758&gi=&mode=2&is_locus=1&is_locus_tag=1&optcols=1,1,1,0,0&replicons=52840256,NC_002942.5,chr"

    end = False
    page = 1
    page_size = 100

    locus_tag_index = 7
    column_mapping = {
        5: "db_xref",
        6: "gene",
        8: "protein_id",
        11: "product"
    }
    nested_values = [5, 8]

    dictionary = {}

    while not end:
        # o crawling acaba quando a página que pedi for vazia

        try:
            url = url_prefix + "&page=" + str(page) + "&pageSize=" + str(page_size)
            tree = fecth_xml_tree(url, add_root=True)
            rows = tree.findall(".TR")

            for row in rows:
                cols = row.findall(".TD")
                locus_tag = cols[locus_tag_index].text

                # criar um dicionário vazio para esta locus tag
                dictionary[locus_tag] = {}

                # para cada uma das colunas que nos interessa,
                # guardar o valor no dicionário se for
                # diferente de -
                for index in column_mapping:
                    prop = column_mapping[index]
                    if index in nested_values:
                        value = cols[index].find(".a").text
                    else:
                        value = cols[index].text

                    if value != "-":
                        dictionary[locus_tag][prop] = value

            page += 1
            end = len(rows) == 0
        except lxml.etree.ParserError:
            end = True

    return dictionary

def fetch_uniprot(uniprot_id):
    """
    Faz download da informação presente no site uniprot
    em formato xml dado um 'uniprot_id'.
    É retornado um object com essa informação.
    """

    url_prefix = "http://www.uniprot.org/uniprot/"
    url = url_prefix + uniprot_id + ".xml"
    tree = fecth_xml_tree(url)
    print(tree)
