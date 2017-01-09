from Bio import SeqIO
from Bio.Blast import NCBIXML
import json, os

def read_file(file_path):
    """
    Lê file_path e retorna uma string com o conteúdo do ficheiro.
    """
    fd = open(file_path, "r")
    string = fd.read()
    fd.close()
    return string

def write_file(string, file_path):
    """
    Grava a string em file_path
    """
    fd = open(file_path, "w")
    fd.write(string)
    os.fsync(fd)
    fd.close()

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
    record = NCBIXML.read(fd)
    return record

def write_blast(handle, file_path):
    """
    Grava o blast record em file_path.
    """
    fd = open(file_path, "w")
    fd.write(handle.read())
    os.fsync(fd)
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
    os.fsync(fd)
    fd.close()

def wrap_file(start_line, end_line, file_path):
    """
    Adiciona uma linha no início do ficheiro e outra no fim.
    NOTA: not memory-friendly
    """
    fd = open(file_path, "r")
    lines = [start_line] + fd.readlines() + [end_line]
    fd.close()
    fd = open(file_path, "w")
    fd.writelines(lines)
    os.fsync(fd)
    fd.close()

def truncate_file(start, end, file_path):
    """
    Remove as 'start' primeiras linhas
    e as 'end' últimas linahs de um ficheiro.
    """
    fd = open(file_path, "r")
    lines = fd.readlines()
    fd.close()
    fd = open(file_path, "w")
    fd.writelines(lines[start:-end])
    os.fsync(fd)
    fd.close()
