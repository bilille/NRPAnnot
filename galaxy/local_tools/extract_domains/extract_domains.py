#!/usr/bin/env python3

import os
import sys
import logging
import argparse
import textwrap
import json
import utils
import collections
import re


def geneclusters2fasta(geneclusters_file, a_domains_fasta, ce_domains_fasta):
    """
    Take the geneclusters.json file (from antiSMASH) and return 2 FASTA files:
    A domains fasta
    C/E domains fasta
    """

    logger.info('Parse geneclusters JSON file')
    logger.debug('File location: %s', geneclusters_file)

    if not os.path.isfile(geneclusters_file):
        logger.fatal('The geneclusters JSON file does not exist')
        sys.exit('Missing file')

    json_data_handler = open(geneclusters_file, 'r')
    json_data = json.load(json_data_handler)
    run_id = json_data.get('run_id', 'runid-undef')
    geneclusters_json, details_data = json_data['geneclusters'], json_data['details_data']

    with open(a_domains_fasta, 'w') as a_domains_handler, \
            open(ce_domains_fasta, 'w') as ce_domains_handler:
        for cluster_name, cluster in details_data.items():
            logger.info('Current cluster:%s', cluster_name)
            for orf in cluster['orfs']:
                logger.info('\tCurrent orf:%s', orf['id'])
                domain_func_idx = collections.Counter()
                for domain_idx, domain in enumerate(orf['domains']):
                    logger.info('\t\tCurrent domain:%s [type=%s]', domain_idx, domain['type'])
                    function = utils.domain_function(domain['type'])
                    seq_name = utils.retrieve_sequence_name(cluster_name, geneclusters_json)
                    opt_params = {'seq_name': seq_name, 'run_id': run_id, 'type': domain['type']}
                    fasta_id = utils.build_fasta_header(cluster_name,
                                                        orf['id'],
                                                        domain_idx,
                                                        function, domain_func_idx[function],
                                                        opt_params)
                    domain_func_idx[function] += 1
                    fstring = fasta_string(fasta_id, domain['sequence'])
                    if function == 'A':
                        print(fstring, file=a_domains_handler)
                    elif function == 'C' or function == 'E':
                        print(fstring, file=ce_domains_handler)
    logger.info('Parsing geneclusters JSON file: Done')


def fasta_string(header, sequence, seq_len=70):
    fstring = '>%s\n%s' % (header, "\n".join(textwrap.wrap(sequence, width=seq_len)))
    return fstring


if __name__ == '__main__':
    logging.basicConfig(level=logging.INFO)
    logger = logging.getLogger()

    # parse command line
    parser = argparse.ArgumentParser()

    parser.add_argument(
        '--geneclusters',
        help='Geneclusters JSON file exported by antiSMASH',
        required=True,
        type=argparse.FileType('r')
    )
    parser.add_argument(
        '--a_fasta',
        help='A domains (fasta format)',
        required=True,
        type=argparse.FileType('w')
    )
    parser.add_argument(
        '--ce_fasta',
        help='c/E domains (fasta format)',
        required=True,
        type=argparse.FileType('w')
    )

    args = parser.parse_args()

    geneclusters2fasta(args.geneclusters.name,
                       args.a_fasta.name,
                       args.ce_fasta.name)
