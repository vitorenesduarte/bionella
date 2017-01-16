from Bio.Blast.Applications import NcbiblastpCommandline
import shutil, os, subprocess
import util.rw as rw
import util.www as www

def fasta_it(tag):
    """
    Retorna o nome do in_file dada a locus tag da proteina.
    """
    return tag + ".fasta"

def xml_it(in_file):
    """
    Retorna o nome do out_file dado o in_file.
    """
    return in_file + ".xml"

def write_queries_to_dir(tags_and_proteins, directory):
    """
    Grava as proteínas passadas como argumento na
    diretoria também passada como argumento.
    """

    tag_to_files = {}

    # Apagar a diretoria e criar uma nova
    if os.path.exists(directory):
        shutil.rmtree(directory)
    os.makedirs(directory)

    for (tag, protein) in tags_and_proteins:
        # Gravar a proteína num ficheiro
        in_file = directory + "/" + fasta_it(tag)
        rw.write_file(protein, in_file)
        
        # registar esta informação num dicionário
        out_file = directory + "/" + xml_it(fasta_it(tag))
        tag_to_files[tag] = (in_file, out_file)

    return tag_to_files

def local_blastp(tag_to_files, db):
    """
    Corre o blast localmente.
    """

    for tag in tag_to_files:
        (in_file, out_file) = tag_to_files[tag]
    
        blastp_cline = NcbiblastpCommandline(
            query=in_file,
            db=db,
            evalue=10,
            outfmt=5,
            out=out_file
        )
        blastp_cline()

def docker_blastp(directory, db):
    """
    Corre o blast numa instância docker.
    """
    cmd = "docker run -e QUERY_DIR=" + directory \
                  + " -e DB=" + db \
                  + " -v $PWD/" + directory + ":/" + directory \
                  + " -ti vitorenesduarte/swissprot_blast"

    p = subprocess.Popen(cmd, shell=True, stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
    p.wait()

def expasy_blastp(tag_to_files):
    """
    Corre o blast na expasy.
    """
    for tag in tag_to_files:
        (in_file, out_file) = tag_to_files[tag]
        query = rw.read_file(in_file)
        blast_result = www.expasy_blast(query)
        rw.write_file(blast_result, out_file)

def blastp(tags_and_proteins, db, type="local"):
    """
    Corre o blast para a proteínas passadas como argumento,
    contra a base de dados também passada como argumento.

    O argumento type é opcional e pode ter dois valores:
        - local
        - docker
        - expasy

    A cada hit do blast extraímos:
        - uniprot_id
        - evalue
        - score
        - identity

    É retornado um dicionário com os resultados.
    """

    directory = ".query_dir"
    tag_to_files = write_queries_to_dir(tags_and_proteins, directory)

    # correr o blast
    if type == "local":
        local_blastp(tag_to_files, db)
    elif type == "docker":
        docker_blastp(directory, db)
    elif type == "expasy":
        expasy_blastp(tag_to_files)
    else:
        raise Exception("Unsupported type: " + type)

    blast_results = extract_blast_info(tag_to_files, type)
    return blast_results

def extract_blast_info(tag_to_files, type):
    """
    Extrai a informação que necessitamos dos resultados do blast.
      - uniprot_id
      - evalue
      - score
      - identity
    """

    blast_results = {}

    for tag in tag_to_files:
        (_, out_file) = tag_to_files[tag]
        handle = rw.read_blast(out_file)

        result = []

        for a in handle.alignments:
            # extrair o uniprot id
            if type in ["local", "docker"]:
                # nos blast locais, o uniprot_id está no hit_def
                uniprot_id = a.hit_def.split("|")[1]
            elif type == "expasy":
                # nos blast expasy, o uniprot_id está no hit_id
                uniprot_id = a.hit_id.split("|")[1]

            # escolher sempre o primeiro hsp
            hsp = a.hsps[0]
            evalue = hsp.expect
            score = hsp.score
            identities = hsp.identities
            align_length = hsp.align_length
            identity = (identities * 100) / align_length

            hit = {}
            hit["uniprot_id"] = uniprot_id
            hit["evalue"] = evalue
            hit["score"] = score
            hit["identity"] = identity
            result.append(hit)

        blast_results[tag] = result

    return blast_results
