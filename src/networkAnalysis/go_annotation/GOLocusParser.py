__author__ = 'guorongxu'

import re
import logging

## Parsing the go gene file and return GO ID list with all gene lists.
def parse_gene_info_file(gene_info_file):
    entrez_id_list = {}

    with open(gene_info_file) as fp:
        lines = fp.readlines()

        for line in lines:
            if line.startswith("#Format"):
                continue
            fields = re.split(r'\t', line)
            entrez_id_list.update({fields[1]:fields[2]})

    fp.closed

    return entrez_id_list

## Parsing the go gene file and return GO ID list with all gene lists.
def parse_go_gene_file(go_gene_file, entrez_id_list):
    go_gene_list = {}
    gene_names = []

    with open(go_gene_file) as fp:
        lines = fp.readlines()

        for line in lines:
            if line.startswith("GOterm"):
                gene_names = re.split(r'\t', line)
            else:
                gene_list = []
                go_id = ""
                fields = re.split(r'\t', line)
                for index, elem in enumerate(fields):
                    if index == 0:
                        go_id = elem
                    if index > 0 and elem == "1":
                        if gene_names[index] in entrez_id_list:
                            gene_list.append(entrez_id_list.get(gene_names[index]))

                go_gene_list.update({go_id: gene_list})

    #print len(go_gene_list.get("GO:0008150"))

    fp.closed

    return go_gene_list, len(gene_names) - 1

## Parsing the go term file and return GO ID list with descriptions.
def parse_go_term_file(go_term_file):
    GO_Term_list = {}

    with open(go_term_file) as fp:
        lines = fp.readlines()

        accepted = False
        go_term_id = ""
        go_term_def = ""
        for line in lines:
            if line.startswith("[Term]"):
                accepted = True
                continue

            if accepted and line.startswith("id: GO:"):
                go_term_id = line[len("id: "):].rstrip()
            if accepted and line.startswith("name: "):
                go_term_name = line[len("name: "):].rstrip()
            if accepted and line.startswith("def: "):
                go_term_def = line[len("def: \""):line.index("[") - 2].rstrip()
            if accepted and line.startswith("\n"):
                GO_Term_list.update({go_term_id: [go_term_name, go_term_def]})
                accepted = False
                go_term_id = ""
                go_term_def = ""

    fp.closed

    return GO_Term_list, len(GO_Term_list)

## Parsing go gene file and go term file.
def parse(go_gene_file, gene_info_file, go_term_file):

    logging.info("parsing GO term file...")
    GO_Term_list, total_unique_term = parse_go_term_file(go_term_file)

    logging.info("parsing GO gene info files...")
    entrez_id_list = parse_gene_info_file(gene_info_file)

    logging.info("parsing GO gene and term file...")
    go_gene_list, total_unique_gene = parse_go_gene_file(go_gene_file, entrez_id_list)

    print str(len(go_gene_list)) + ":" + str(total_unique_gene) + ":" + str(len(GO_Term_list))
    return go_gene_list, total_unique_gene, GO_Term_list

## Main entry
if __name__ == "__main__":
    #gene_info_file = sys.argv[1]
    #gene_term_file = sys.argv[2]
    go_gene_file = "/Users/guorongxu/Desktop/SearchEngine/GO/GO2all_locus.txt"
    gene_info_file = "/Users/guorongxu/Desktop/SearchEngine/GO/Homo_sapiens.gene_info"
    go_term_file = "/Users/guorongxu/Desktop/SearchEngine/GO/go.obo"

    parse(go_gene_file, gene_info_file, go_term_file)