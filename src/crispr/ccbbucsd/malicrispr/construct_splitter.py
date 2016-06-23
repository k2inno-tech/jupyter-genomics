__author__ = 'Amanda Birmingham'
__maintainer__ = "Amanda Birmingham"
__email__ = "abirmingham@ucsd.edu"
__status__ = "prototype"


def expand_gene_info(dataframe, header_mgr, columns_list, neg_control_genes):
    if not (header_mgr.position1_header in columns_list and header_mgr.position2_header in columns_list):
        # TODO: can I blithely "apply" a method or do I have to explicitly account for "self"?
        expanded_tuples_df = dataframe.apply(lambda x: _expand_info_from_construct_id(x,
            header_mgr.construct_id_header, neg_control_genes), axis=1)
        (dataframe[header_mgr.position1_header], dataframe[header_mgr.position2_header],
         dataframe[header_mgr.gene_pair_header]) = zip(*expanded_tuples_df)


def _expand_info_from_construct_id(construct_row, construct_id_header, neg_control_genes):
    num_expected_genes = 2
    pair_separator = "__"
    piece_separator = "_"

    construct_id = construct_row[construct_id_header]
    construct_id_pieces = construct_id.split(pair_separator)

    if len(construct_id_pieces) != num_expected_genes:
        raise ValueError("Construct {0} split into {1} pieces rather than {2}".format(construct_id,
                                                                                      len(construct_id_pieces),
                                                                                      num_expected_genes))

    results = []
    for index in range(0, num_expected_genes):
        gene_pieces = construct_id_pieces[index].split(piece_separator)
        gene = gene_pieces[0]

        # TODO: Double-check with customer if this is correct:
        # if the "gene" name is for one of the non-targeting
        # controls, compress to the control's name
        for control_gene in neg_control_genes:
            if gene.startswith(control_gene):
                gene = control_gene

        results.append(gene)

    # Ensure that gene pairs are treated as combinations not
    # permutations (order not important) by sorting gene list
    # alphabetically
    gene_pair = pair_separator.join(sorted(results))

    results.append(gene_pair)  # UNsorted
    return tuple(results)