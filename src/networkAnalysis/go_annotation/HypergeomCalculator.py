__author__ = 'James & Guorong'

import logging
import itertools
import math
from fdr import fdr
import GOLocusParser
#import GOParser
from scipy.stats import hypergeom
from datetime import datetime

def calc_pvalue(gene_list, gene_set, M):
    gene_list = set(gene_list)
    gene_set = set(gene_set)

    N = len(gene_list)
    n = len(gene_set)
    overlap = gene_list & gene_set
    k = len(overlap)

    return hypergeom(M, n, N).sf(k), list(overlap)

def calc_enrichment(gene_list, GO_ID_list, total_unique_gene, GO_Term_list):
    M = total_unique_gene

    enriched_list = []
    for term in GO_ID_list:
        if len(GO_ID_list.get(term)) >= 20 and len(GO_ID_list.get(term)) <= 2000:
            pvalue, overlap = calc_pvalue(gene_list, GO_ID_list.get(term), M)
            if len(overlap) > 1:
                enriched_item = {"go_id": term, "name":GO_Term_list.get(term)[0] ,"description":GO_Term_list.get(term)[1],
                                 "pvalue": pvalue, "overlap": overlap, "genes_from_list": len(gene_list), "genes_from_go": len(GO_ID_list.get(term))}
                enriched_list.append(enriched_item)

    enriched_list.sort(key=lambda it: it['pvalue'])

    for qvalue, it in itertools.izip(fdr([it['pvalue'] for it in enriched_list], presorted=True), enriched_list):
        if math.fabs(qvalue) == 0:
            it['qvalue'] = float("inf")
        else:
            it['qvalue'] = -math.log(qvalue, 10)

    enriched_list.sort(key=lambda it: it['qvalue'], reverse=True)

    return enriched_list

## Main entry
if __name__ == "__main__":
    go_gene_file = "/Users/guorongxu/Desktop/SearchEngine/GO/GO2all_locus.txt"
    gene_info_file = "/Users/guorongxu/Desktop/SearchEngine/GO/Homo_sapiens.gene_info"
    go_term_file = "/Users/guorongxu/Desktop/SearchEngine/GO/go.obo"

    gene_list = ['hsa-mir-3167', 'GUCA2B', 'DONSON', 'CDC37', 'ARPC2', 'HSPH1', 'IL4',
                 'TOR1AIP2', 'SUOX', 'SEZ6', 'RORB', 'RAB33A', 'POPDC3', 'PHF7', 'PARP9',
                 'OR10H4', 'LRTM2', 'ACTL7A', 'ZNF30', 'CDHR3', 'FAM189A1', 'HNRNPA2B1',
                 'CNGA2', 'CHD3', 'CAMKMT', 'ASAH1', 'HEMK1', 'SVIL', 'SCGB2A2', 'OR5I1',
                 'TFE3', 'CTNNBL1', 'PCDHB4', 'RNASE1', 'ZBTB14', 'COPS4', 'NLGN3',
                 'MPHOSPH10', 'NLRP9', 'MMP15', 'NXPE3', 'NOX1', 'TRMT10B', 'KLHL2',
                 'FBXL7', 'KRT38', 'FGFR2', 'GRAMD4', 'TRIO', 'ELAVL3', 'CD6', 'KIAA1586',
                 'CEP104', 'CHST3', 'PCDHA13', 'TNC', 'SCN9A', 'PLEKHA5', 'VIM', 'PDE1C',
                 'SLITRK6', 'CPAMD8', 'PCDHGA10', 'ADCY5', 'TOP2B', 'EPB41L3', 'NOX3',
                 'GPHN', 'KCNT2', 'MME']

    logging.info("parsing GO gene and term files...")

    GO_ID_list, total_unique_gene, GO_Term_list = GOLocusParser.parse(go_gene_file, gene_info_file, go_term_file)
    #GO_ID_list, total_unique_gene, GO_Term_list = GOParser.parse(go_gene_file, go_term_file)
    enriched_list = calc_enrichment(gene_list, GO_ID_list, total_unique_gene, GO_Term_list)

    print enriched_list