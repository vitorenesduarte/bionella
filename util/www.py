from Bio import Entrez, SeqIO
from urllib.request import urlopen
import lxml.html, lxml.etree

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
    faz um pedido http ao link_prefix para preencher a tabela.
    Este método implementa a paginação que seria feita manualmente
    no website, e guarda a informação que queremos num dicionário.

      - A coluna 5  corresponde à propriedade db_xref    do genbank
      - A coluna 6  corresponde à propriedade gene       do genbank
      - A coluna 7  corresponde à propriedade locus_tag  do genbank
      - A coluna 8  corresponde à propriedade protein_id do genbank
      - A coluna 11 corresponde à propriedade product    do genbank
    """

    link_prefix = "https://www.ncbi.nlm.nih.gov/genomes/Genome2BE/genome2srv.cgi?action=GetFeatures4Grid&type=Proteins&genome_id=416&genome_assembly_id=166758&gi=&mode=2&is_locus=1&is_locus_tag=1&optcols=1,1,1,0,0&replicons=52840256,NC_002942.5,chr"

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

    dictionary = {}

    while not end:
        # o crawling acaba quando a página que pedi for vazia

        try:
            link = link_prefix + "&page=" + str(page) + "&pageSize=" + str(page_size)
            link_content = urlopen(link).read()
            doc = lxml.html.fromstring(link_content)
            rows = doc.getchildren()

            for row in rows:
                cols = row.getchildren()
                locus_tag = cols[locus_tag_index].text_content()

                # criar um dicionário vazio para esta locus tag
                dictionary[locus_tag] = {}

                # para cada uma das colunas que nos interessa,
                # guardar o valor no dicionário se for
                # diferente de -
                for index in column_mapping:
                    prop = column_mapping[index]
                    value = cols[index].text_content()

                    if value != "-":
                        dictionary[locus_tag][prop] = value

            page += 1
            end = len(rows) == 0
        except lxml.etree.ParserError:
            end = True

    return dictionary
