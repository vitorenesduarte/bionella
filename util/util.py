import hashlib

def md5(string):
    """
    Retorna a hash da string passada como argumento.
    """
    return hashlib.md5(string.encode("utf-8")).hexdigest()
