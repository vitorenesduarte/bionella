import hashlib

def md5(string):
    """
    Retorna a hash da string passada como argumento.
    """
    return hashlib.md5(string.encode("utf-8")).hexdigest()

def get_tags_and_proteins(dictionary):
    """
    Esta função retorna uma lista de pares
    em que a primeira componente é a tag
    e a segunda a sequência de aminoácidos.
    """

    tags_and_proteins = []

    for tag in dictionary:
        if dictionary[tag]["type"] == "mRNA":
            # se for do tipo "mRNA" estamos na presença de uma proteína
            protein = dictionary[tag]["translation"]
            both = (tag, protein)
            tags_and_proteins.append(both)

    return sorted(tags_and_proteins, key=lambda tp : tp[0])
