from Bio.Blast.Applications import NcbiblastpCommandline
import hashlib, shutil, os, subprocess
import util.helper as helper

def md5(protein):
    """
    Retorna a hash da proteina passada como argumento
    """
    return hashlib.md5(protein.encode("utf-8")).hexdigest()

def write_queries_to_dir(proteins, directory):
    """
    Grava as proteínas passadas como argumento na
    diretoria também passada como argumento.
    Retorna uma lista de pares onde a primeira
    componente é o ficheiro onde foi guardada
    a proteína e a segunda componente é onde
    deve ser guardado o resultado do blast.
    """

    result = []

    # Apagar a diretoria e criar uma nova
    if os.path.exists(directory):
        shutil.rmtree(directory)
    os.makedirs(directory)

    for protein in proteins:
        # Gravar a proteína num ficheiro
        h = md5(protein)
        in_file = directory + "/" + h + ".fasta"
        out_file = directory + "/" + h + ".xml"
        helper.write_file(protein, in_file)

        # Adicionar o par à lista de resultados
        pair = (in_file, out_file)
        result.append(pair)

    return result

def local_blastp(in_file, out_file, db):
    """
    Corre o blast localmente.
    """
    blastp_cline = NcbiblastpCommandline(
            query=in_file,
            db=db,
            evalue=10,
            outfmt=5,
            out=out_file
    )
    blastp_cline()

def docker_blastp(directory, in_file, out_file, db):
    """
    Corre o blast numa instância docker.
    """
    cmd = "docker run -e IN_FILE=" + in_file \
                  + " -e OUT_FILE=" + out_file \
                  + " -e DB=" + db \
                  + " -v $PWD/" + directory + ":/" + directory \
                  + " -ti vitorenesduarte/swissprot_blast"

    p = subprocess.Popen(cmd, shell=True, stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
    p.wait()

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
    pair_list = write_queries_to_dir(proteins, directory)

    for (in_file, out_file) in pair_list:
        if type == "local":
            local_blastp(in_file, out_file, db)
        elif type == "docker":
            docker_blastp(directory, in_file, out_file, db)
        else:
            raise Exception("Unsupported type: " + type)

    # Retornar os out_files
    return [out_file for (in_file, out_file) in pair_list]
