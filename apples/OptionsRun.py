from apples.OptionsBasic import OptionsBasic
import logging


def options_config():
    """
    This function configures the options for the given parameters and returns the options and arguments.
    """
    parser = OptionsBasic('jplace')
    parser.add_option('-a', '--database', dest='database_fp', metavar='FILE', help='path to the APPLES database')
    parser.add_option(
        '-d', '--distances', dest='dist_fp', help='path to the table of observed distances', metavar='FILE'
    )
    parser.add_option(
        '-x',
        '--extendedref',
        dest='extended_ref_fp',
        help='path to the extened reference alignment file (FASTA), ' 'containing reference and query sequences',
        metavar='FILE',
    )
    parser.add_option(
        '-q',
        '--query',
        dest='query_fp',
        help='path to the query alignment file (FASTA), containing query sequences',
        metavar='FILE',
    )
    parser.add_option(
        '-m',
        '--method',
        dest='method_name',
        default='FM',
        help='name of the weighted least squares method (OLS, FM, BME, or BE)',
        metavar='METHOD',
    )
    parser.add_option(
        '-c',
        '--criterion',
        dest='criterion_name',
        default='MLSE',
        help='name of the placement selection criterion (MLSE, ME, or HYBRID',
        metavar='CRITERIA',
    )
    parser.add_option(
        '-n',
        '--negative',
        dest='negative_branch',
        action='store_true',
        help='relaxes positivity constraint on new branch lengths, i.e. allows negative branch lengths',
    )
    parser.add_option(
        '-b',
        '--base',
        dest='base_observation_threshold',
        type=int,
        default=25,
        help='minimum number of observations kept for ' 'each query ignoring the filter threshold.',
        metavar='NUMBER',
    )
    parser.add_option(
        '-V',
        '--overlap',
        dest='minimum_alignment_overlap',
        type=float,
        default=0.001,
        help='minimum fraction of nongap sites needed for a valid pairwise distance.',
        metavar='NUMBER',
    )
    parser.add_option(
        '-X',
        '--mask',
        dest='mask_lowconfidence',
        action='store_true',
        default=False,
        help='masks low confidence characters in the alignments indicated by lowercase characters '
        'output by softwares like SEPP.',
    )
    parser.add_option(
        '--exclude',
        dest='exclude_intplace',
        action='store_true',
        default=False,
        help='exclude queries placed on the internal nodes in jplace file.',
    )
    parser.add_option("-S", "--support", dest="find_support", action='store_true', default=False,
                      help="adds support value for the placed queries.")                    
    parser.add_option("-F", "--fast", dest="fast_support", action='store_true', default=False,
                      help="enables fast bootstrapping (True by default if reestimation is turned off).")                    
    parser.add_option("-N", "--sample", dest="sample_count", type=int, default=100,
                      help="number of bootstrapping samples.", metavar="NUMBER")
    parser.add_option("--seed", dest="bootstrapping_seed", type=int, default=56,
                      help="seed for bootstrapping.", metavar="NUMBER") 
    parser.add_option("--keep-factor", dest="keep_factor", type=float, default=0.01,
                      help="throw away anything that has ml_ratio below keep_factor times (best ml_ratio).", metavar="NUMBER")                                                                    
    parser.add_option("--keep-at-most", dest="keep_at_most", type=int, default=5,
                      help="maximum number of placements to be kept.", metavar="NUMBER") 
    parser.add_option("--lse", dest="prioritize_lse", action='store_true', default=False,
                      help="Only keep the placement with the minimum lse.")

    (options, args) = parser.parse()

    if options.dist_fp:
        options.reestimate_backbone = False
        if options.ref_fp:
            raise ValueError('Input should be either an alignment or a distance matrix, but not both!')
        if options.database_fp:
            logging.warning(
                'Input contains both an APPLES database and a distance matrix. Database sequences '
                'will be ignored. Database phylogeny will be used if user did not provide a phylogeny '
                '(using -t option). '
            )

    if options.database_fp:
        if options.ref_fp:
            raise ValueError('Input should be either an alignment or a APPLES database file, but not both!')
        if options.tree_fp:
            logging.warning(
                'Input contains both an APPLES database and a tree file. User provided tree has '
                'higher priority and therefore will be used.'
            )
    if not options.tree_fp and not options.database_fp:
        raise ValueError('No input backbone tree provided by user.')
    if options.query_fp and options.extended_ref_fp:
        raise ValueError('Input should be either an extended alignment or a query alignment, but not both!')

    if options.find_support and options.dist_fp:
        options.find_support = False
        logging.warning("To find support, input should be an alignment.")

    if options.find_support:
        # if not options.fast_support and options.disable_reestimation:
        #     options.fast_support = True

        if options.sample_count <= 0:
            raise ValueError('Sample count has to be a positive integer.')
    else:
        options.fast_support = False

    return options, args
