from Bio import SeqIO
from Bio.Blast import NCBIXML
import json

def read_genbank(file_path):
    """
    Lê file_path e retorna um record genbank.
    """
    return SeqIO.read(file_path, "gb")

def write_genbank(record, file_path):
    """
    Grava o record em file_path no formato genbank.
    """
    SeqIO.write(record, file_path, "gb")

def read_blast(file_path):
    """
    Lê file_path e retorna uma lista de blast records.
    """
    fd = open(file_path, "r")
    record = NCBIXML.parse(fd)
    return record

def write_blast(handle, file_path):
    """
    Grava o blast record em file_path.
    """
    print(handle)
    records = [r for r in handle]
    print(records)
    print(len(records))
    assert len(records) == 1
    record = records[0]
    fd = open(file_path, "w")
    fd.write(record.read())
    fd.close()

def read_json(file_path):
    """
    Lê file_path e retorna dicionário.
    """
    fd = open(file_path, "r")
    dictionary = json.load(fd)
    fd.close()

    return dictionary

def write_json(dictionary, file_path):
    """
    Grava o dicionário em file_path no formato json.
    """
    fd = open(file_path, "w")
    json.dump(dictionary, fd, sort_keys=True, indent=2)
    fd.close()
