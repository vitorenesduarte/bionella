from collections import defaultdict
import util.www as www
import util.rw as rw
import util.blast as blast
import util.util as util

def add_info_to_blast_results(blast_results, uniprots):
    """
    Esta função adiciona aos blast results informação
    retirada do site da uniprot.
    """

    for tag in blast_results:
        for i in range(len(blast_results[tag])):
            uniprot_id = blast_results[tag][i]["uniprot_id"]

            # esta proteína aparece nos resultados do blast como
            # Q8TC84-2 mas a uniprot retorna a sua informação
            # como Q8TC84
            if uniprot_id == "Q8TC84-2":
                uniprot_id = "Q8TC84"

            properties = [
                "status",
                "organism",
                "sequence",
                "molecular_functions"
            ]

            for p in properties:
                value = uniprots[uniprot_id][p]
                blast_results[tag][i][p] = value

    return blast_results

def infer_function(blast_results):
    """
    Esta função tenta inferir a função das proteínas baseando-se
    nos resultados do blast.
    """

    all = {}

    for tag in blast_results:
        inferred = []
        leaderboard = defaultdict(int)

        # Para os primeiros 10 resultados contar as ocorrências
        # de cada função
        top_ten = blast_results[tag][0:10]
        for result in top_ten:
            for function in result["molecular_function"]:
                leaderboard[function] += 1

        # As funções que aparecem em pelo menos 50% dos resultados
        # são consideradas potenciais funções
        min = len(leaderboard) / 2

        for function in leaderboard:
            if leaderboard[function] >= min:
                inferred.append(function)

        all[tag] = inferred

    return all

def main():
    ncbi_json_path = ".ncbi.json"
    blast_results_json_path = ".blast_results.json"
    uniprots_json_path = ".uniprots.json"

    #dictionary = rw.read_json(ncbi_json_path)
    #tags_and_proteins = util.get_tags_and_proteins(dictionary)

    #blast_results = blast.blastp(
    #    tags_and_proteins,
    #    "expasy"
    #)
    #rw.write_json(blast_results, blast_results_json_path)

    blast_results = rw.read_json(blast_results_json_path)

    #uniprot_ids = set()
    #for tag in blast_results:
    #    for result in blast_results[tag]:
    #        uniprot_ids.add(result["uniprot_id"])

    #uniprot_ids = list(uniprot_ids)

    #uniprots = www.fetch_uniprots(uniprot_ids)
    #rw.write_json(uniprots, uniprots_json_path)

    uniprots = rw.read_json(uniprots_json_path)

    blast_results = add_info_to_blast_results(
        blast_results,
        uniprots
    )
    rw.write_json(blast_results, blast_results_json_path)
    

if __name__ == "__main__":
    main()
