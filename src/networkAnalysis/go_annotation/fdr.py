__author__ = 'James'

def fdr(pvalues, presorted=False):
    n = float(len(pvalues))
    if (presorted):
        qvalues = [n / (rank + 1) * pvalue for rank, pvalue in enumerate(pvalues)]
    else:
        pvalues = [(pvalue, k) for k, pvalue in enumerate(pvalues)]
        pvalues = sorted([(pvalue, k) for k, pvalue in enumerate(pvalues)])
        qvalues = sorted([(k, n / (rank + 1) * pvalue) for rank, (pvalue, k) in enumerate(pvalues)])
        qvalues = [qvalue for _, qvalue in qvalues]

    return qvalues