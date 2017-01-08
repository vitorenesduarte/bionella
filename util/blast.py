from Bio.Blast.Applications import NcbiblastpCommandline
import hashlib, shutil, os, subprocess
import util.helper as helper

def md5(protein):
    """
    Retorna a hash da proteina passada como argumento.
    """
    return hashlib.md5(protein.encode("utf-8")).hexdigest()

def fasta_it(md5):
    """
    Retorna o nome do in_file dada a hash da proteina.
    """
    return md5 + ".fasta"

def xml_it(in_file):
    """
    Retorna o nome do out_file dado o in_file.
    """
    return in_file + ".xml"

def get_out_files(directory):
    """
    Retorna a lista de out_files dada a directoria
    com os in_files.
    """
    in_files = os.listdir(directory)
    return [xml_it(in_file) for in_file in in_files]

def write_queries_to_dir(proteins, directory):
    """
    Grava as proteínas passadas como argumento na
    diretoria também passada como argumento.
    """

    # Apagar a diretoria e criar uma nova
    if os.path.exists(directory):
        shutil.rmtree(directory)
    os.makedirs(directory)

    for protein in proteins:
        # Gravar a proteína num ficheiro
        h = md5(protein)
        in_file = directory + "/" + fasta_it(h)
        helper.write_file(protein, in_file)

def local_blastp(directory, db):
    """
    Corre o blast localmente.
    """

    os.chdir(directory)
    in_files = os.listdir()

    for in_file in in_files:
        out_file = xml_it(in_file)
    
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

    print(p.stdout.read())
    print(p.stderr.read())

def blastp(proteins, db, type="local"):
    """
    Corre o blast para a proteínas passadas como argumento,
    contra a base de dados também passada como argumento.
    O argumento type é opcional e pode ter dois valores:
        - local
        - docker
    Cada um dos blast é gravado num ficheiro xml.
    É retornada uma lista com os ficheiros xml.
    """

    directory = ".query_dir"
    write_queries_to_dir(proteins, directory)
    out_files = get_out_files(directory)

    if type == "local":
        local_blastp(directory, db)
    elif type == "docker":
        docker_blastp(directory, db)
    else:
        raise Exception("Unsupported type: " + type)

    return out_files
