#!/usr/bin/env python3

import sys
import logging
import argparse
import re
import os
import utils
import collections
import json

def nrps2tsv(nrps_dir):
    """
    Format nrpspredictor2 results from antiSMASH as TSV file
    """
    if not os.path.isdir(nrps_dir):
        logger.fatal("The NRPS directory does not exists: %s", nrps_dir)
        sys.exit()
    nrps_files = [os.path.join(nrps_dir,f) for f in os.listdir(nrps_dir) if f.endswith('nrpspredictor2_svm.txt')]

    tsv = []
    header = []
    for nrps_file in nrps_files:
        with open(nrps_file, 'r') as nrps_handler:
            lines = nrps_handler.readlines()
            if not header:
                header = lines[0].replace('<tab>', '\t')
                tsv.append(header)
            tsv.extend(lines[1:])
    return tsv


def replace_antismash_id(tsv, json_data_handler):
    """
    Replace NRPSPredictor2 id from antiSMASH by the id used by us.
    (Enable the capability to match the prediciton to a given domain.)
    """
    json_data = json.load(json_data_handler)
    run_id = json_data.get('run_id', 'runid-undef')
    geneclusters_json, details_data = json_data['geneclusters'], json_data['details_data']

    if not tsv:
        logger.info('No NRPSPredictor2 predicitons found')
        sys.exit()
    updated_tsv = [tsv[0]]
    for cluster_name, cluster in details_data.items():
        for orf in cluster['orfs']:
            domain_func_idx = collections.Counter()
            for domain_idx, domain in enumerate(orf['domains']):
                function = utils.domain_function(domain['type'])
                if function.upper() != 'A':
                    continue
                seq_name = utils.retrieve_sequence_name(cluster_name, geneclusters_json)
                opt_params = {'seq_name': seq_name, 'run_id': run_id, 'type': domain['type']}
                fasta_id = utils.build_fasta_header(cluster_name,
                                                        orf['id'],
                                                        domain_idx,
                                                        function, domain_func_idx[function],
                                                        opt_params)

                domain_func_idx[function] += 1
                found = False
                for l in tsv[1:]:
                    fields = l.split('\t')
                    arr = fields[0].split('_')
                    domain_id = arr[-1]
                    orf_id = utils.sanitaze(utils.SEP.join(arr[:(len(arr)-1)]))

                    if orf_id == utils.sanitaze(orf['id']) and \
                       domain_id == function + str(domain_func_idx[function]):
                        found = True
                        fields[0] = fasta_id
                        updated_tsv.append('\t'.join(fields))
                        break
                if not found:
                    logger.error('Cannot retrieve nrpspredictor2 prediction:')
                    logger.error('%s| |%s|' % (orf['id'], domain_func_idx[function]))
                    sys.exit()
    return updated_tsv


if __name__ == '__main__':
    logging.basicConfig(level=logging.INFO)
    logger = logging.getLogger(__name__)

    # parse command line
    parser = argparse.ArgumentParser()

    parser.add_argument(
        '--nrps_dir',
        help='antiSMASH NRPS directory',
        required=True,
        type=str
    )
    parser.add_argument(
        '--geneclusters',
        help='Geneclusters JSON file exported by antiSMASH',
        required=True,
        type=argparse.FileType('r')
    )

    parser.add_argument(
        '--tsv_out',
        required=True,
        type=argparse.FileType('w')
    )

    args = parser.parse_args()
    tsv = nrps2tsv(args.nrps_dir)
    ids = set()
    for l in tsv:
        id_ = l.split('\t')[0]
        if id_ in ids:
            logger.error('Ambiguous identifiant:%s', id_)
            sys.exit()
        ids.add(id_)

    tsv = replace_antismash_id(tsv, args.geneclusters)
    for l in tsv:
        args.tsv_out.write(l)
