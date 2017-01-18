from Bio import Entrez, SeqIO
from urllib.request import urlretrieve
from xml.etree import ElementTree
import requests
import util.rw as rw
import util.util as util

def parse_xml(file_path, add_root=False, start=0, end=0):
    """
    Faz parsing de um ficheiro xml e retorna um 'ElementTree'.
    """
    if start > 0 or end > 0:
        rw.truncate_file(start, end, file_path)

    if add_root:
        rw.wrap_file("<root>", "</root>", file_path)

    fd = open(file_path, "r")
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

      - A coluna 5  corresponde à propriedade db_xref    do genbank (gene_id para nós)
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
        5: "gene_id",
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
            file_path, _ = urlretrieve(url)
            tree = parse_xml(file_path, add_root=True)
            rows = tree.findall("TR")

            for row in rows:
                cols = row.findall("TD")
                locus_tag = cols[locus_tag_index].text

                # criar um dicionário vazio para esta locus tag
                dictionary[locus_tag] = {}

                # para cada uma das colunas que nos interessa,
                # guardar o valor no dicionário se for
                # diferente de -
                for index in column_mapping:
                    prop = column_mapping[index]
                    if index in nested_values:
                        value = cols[index].find("a").text
                    else:
                        value = cols[index].text

                    if value != "-":
                        dictionary[locus_tag][prop] = value

            page += 1
            end = len(rows) == 0
        except Exception as ex:
            print(ex)
            end = True

    return dictionary

def gene_ids_to_uniprot_ids(ids):
    """
    Mapeia a lista de gene ids passada como argumento para uniprot ids.
    Retorna um dictionário em que as chaves são gene_id e os valores
    são uniprot_id.
    """

    url = "http://www.uniprot.org/mapping/"
    query = " ".join(ids)

    data = {
        "format": "tab",
        "from": "P_ENTREZGENEID",
        "to": "ACC",
        "query": query
    }

    # HTTP POST
    response = requests.post(
        url,
        data=data
    )

    result = {}

    # ignorar a primeira linha
    for row in response.text.splitlines()[1:]:
        gene_id, uniprot_id = row.split("\t")

        result[gene_id] = uniprot_id

    return result

def expasy_blast(query):
    """
    Corre o blast na expasy.
    """

    url = "http://web.expasy.org/cgi-bin/blast/blast.pl"
    data = {
        "format": "xml",
        "prot_db1": "UniProtKB",
        "matrix": "BLOSUM62",
        "showsc": 50, # best scoring sequences
        "showal": 50, # best alignments
        "seq": query
    }

    response = requests.post(
        url,
        data=data
    )

    return response.text

def fetch_uniprots(ids):
    """
    Faz download da informação presente no site uniprot
    em formato xml dada uma lista de ids uniprot.
    """

    url = "http://www.uniprot.org/uploadlists/"
    max_retries = 100

    # Cada pedido à uniprot vai no máximo com 1000 ids.
    queries  = [ids[i:i+1000] for i in range(0, len(ids), 1000)]

    files = []

    for query in queries:
        query_all = " ".join(query)
        md5 = util.md5(query_all)

        data = {
            "format": "xml",
            "from": "ID",
            "to": "ACC",
            "uploadQuery": query_all,
            "jobId": md5
        }

        done = False
        attempt = 0

        # enquanto não acabar ou não esgotar as tentativas todas
        while not done and attempt < max_retries:

            # HTTP POST
            response = requests.post(
                url,
                data=data
            )

            if response.status_code == 200:
                # Gravar o xml
                file_path = "/tmp/" + md5 + ".xml"
                rw.write_file(response.text, file_path)
                files.append(file_path)
                done = True

            attempt += 1

        if not done:
            print("Não consegui fazer download de " + str(query))

    uniprots = {}

    for file_path in files:
        tree = parse_xml(file_path, add_root=True, start=2, end=1)
        entries = tree.findall(".//entry")

        for entry in entries:
            (uniprot_id, info) = extract_uniprot_info(entry)
            uniprots[uniprot_id] = info

    return uniprots

def extract_uniprot_info(entry):
    """
    Extrai a informação que necessitamos da uniprot.
    """
    # dataset
    ds = entry.get("dataset")

    if ds == "Swiss-Prot":
        status = "reviewed"
    elif ds == "TrEMBL":
        status = "unreviewed"
    else:
        print("Não conheço o dataset " + ds)

    # accessions
    accessions = [a.text for a in entry.findall(".//accession")]
    accession = accessions[0]

    # short name
    short_name = entry.find("name").text

    # full name
    full_name = entry.find(".//fullName").text

    # ec number
    ec_number = entry.find(".//ecNumber")
    if ec_number != None:
        ec_number = ec_number.text

    # organism
    organisms = entry.findall(".//name[@type='scientific']")
    if len(organisms) == 1:
        organism = organisms[0].text
    else:
        print("Encontrei um número de organismos diferente de 1: ", accession)
        organism = None

    # encontrar o texto função que costuma estar no início da
    # página da uniprot
    comment_functions = [c.find("text").text for c in entry.findall(".//comment[@type='function']")]

    # cofator
    cofactors = [c.find("cofactor").find("name").text for c in entry.findall(".//comment[@type='cofactor']")]

    # patologias
    pathologies = [c.find("disease").find("name").text for c in entry.findall(".//comment[@type='disease']")]

    # PDB
    pdbs = []
    PDB = entry.findall(".//dbReference[@type='PDB']")
    for pdb in PDB:
        id = pdb.get("id")
        method = pdb.find("property[@type='method']")
        chains = pdb.find("property[@type='chains']")

        if (not method == None) and (not chains == None):
            method = method.get("value")
            chains = chains.get("value")

            d = {}
            d["id"] = id
            d["method"] = method
            d["chains"] = chains
            pdbs.append(d)

    # encontrar GO - Molecular Function
    # e GO - Biological process
    molecular_functions = []
    biological_processes = []
    locations = []
    GO = entry.findall(".//dbReference[@type='GO']")
    for go in GO:
        value = go.find("property[@type='term']").get("value")
        is_molecular_function = value.startswith("F:")
        is_biological_process = value.startswith("P:")
        is_location = value.startswith("C:")

        value = value[2:]

        if is_molecular_function:
            molecular_functions.append(value)
        elif is_biological_process:
            biological_processes.append(value)
        elif is_location:
            locations.append(value)

    # sequência
    sequences = entry.findall(".//sequence[@length]")
    assert len(sequences) == 1
    sequence = sequences[0].text.replace("\n", "")
    length = int(sequences[0].get("length"))
    mass = int(sequences[0].get("mass"))

    # modified residue
    mrs = []
    res = entry.findall(".//feature[@type='modified residue']")
    for r in res:
        desc = r.get("description")
        position = r.find("location").find("position").get("position")
        mrs.append({"desc": desc, "position": position})

    dictionary = {}
    dictionary["status"] = status
    dictionary["accessions"] = accessions
    dictionary["short_name"] = short_name
    dictionary["product"] = full_name
    dictionary["EC_number"] = ec_number
    dictionary["organism"] = organism
    dictionary["comment_functions"] = comment_functions
    dictionary["cofactors"] = cofactors
    dictionary["pathologies"] = pathologies
    dictionary["pdbs"] = pdbs
    dictionary["molecular_functions"] = molecular_functions
    dictionary["biological_processes"] = biological_processes
    dictionary["locations"] = locations
    dictionary["translation"] = sequence
    dictionary["length"] = length
    dictionary["mass"] = mass
    dictionary["modified_residues"] = mrs

    return (accession, dictionary)
