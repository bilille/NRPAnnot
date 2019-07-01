#!/usr/bin/env python3

import os
import sys
import logging
import argparse
import itertools
import json
import utils
import collections


def build_report(report_directory):
    """
    Complete JSON data from antiSMASH with NaPDoS/NRPSPredictor2 results
    Ouptput JSON in the current report_directory: geneclusters_updated.js
    Extract predictions summary from geneclusters_updated.js
    """

    if not os.path.isdir(report_directory):
        logger.fatal('The following directory does not exist:%s', report_directory)
        sys.exit('Directory issue')

    geneclusters_file = os.path.join(report_directory, 'geneclusters.json')
    if not os.path.isfile(geneclusters_file):
        logger.fatal('The following file is required to complete the process:%s',
                     geneclusters_file)
        sys.exit('Missing file')

    # parse json file
    json_data_handler = open(geneclusters_file, 'r')
    json_data = json.load(json_data_handler)

    run_id = json_data.get('run_id', 'runid-undef')
    geneclusters, details_data = json_data['geneclusters'], json_data['details_data']

    # complete json with napdos results
    napdos_table_file = os.path.join(report_directory, 'napdos', 'napdos.tab')
    if not os.path.isfile(napdos_table_file):
        logger.warning('NaPDoS: the following file is missing:%s',
                       napdos_table_file)
    else:
        add_napdos_data(details_data, napdos_table_file)

    # complete json with nrpspredictor results
    nrpspredictor_table_file = os.path.join(report_directory,
                                            'nrpspredictor',
                                            'nrpspredictor.tab')
    if not os.path.isfile(nrpspredictor_table_file):
        logger.warning('NRPSPredictor: the following file is missing:%s',
                       nrpspredictor_table_file)
    else:
        add_nrpspredictor_data(details_data, nrpspredictor_table_file)

    # write updated json file
    updated_geneclusters_file = os.path.join(report_directory, 'geneclusters_updated.js')
    with open(updated_geneclusters_file, 'w') as u_geneclusters_handler:
        print('var run_id = "' + run_id + '";',
              file=u_geneclusters_handler)
        print('var geneclusters = ' + json.dumps(geneclusters, indent=4) + ';',
              file=u_geneclusters_handler)
        print('var details_data = ' + json.dumps(details_data, indent=4) + ';',
              file=u_geneclusters_handler)

    # extract predictions summary
    extract_predictions_summary(run_id, geneclusters, details_data, report_directory)


def add_napdos_data(details_data, napdos_table_file):
    """
    Only the cluster_name, orf_id, domain_id are needed
    to identify corresponding napdos results
    """

    with open(napdos_table_file, 'r') as napdos_table_handler:
        napdos_table = [l.strip().split('\t') for l in napdos_table_handler.readlines() if l.strip()]
        for cluster_id, cluster in _groupby(napdos_table[1:], lambda row: row[0]):
            fasta_data = utils.split_fasta_header(cluster_id)
            cluster_name = fasta_data['cluster_name']
            orf_id = fasta_data['orf_id']
            domain_id = fasta_data['domain_id']

            header = napdos_table[0]
            cluster_table = list(cluster)
            cluster_table.insert(0, header)
            for orf in details_data[cluster_name]['orfs']:
                if utils.sanitaze(orf['id']) == orf_id:
                    orf['domains'][int(domain_id)]['napdos_details'] = cluster_table
                    break


def add_nrpspredictor_data(details_data, nrpspredictor_table_file):
    """
    Only the cluster_name, orf_id, domain_id are needed
    to identify corresponding nrpspredictor results
    """

    with open(nrpspredictor_table_file, 'r') as nrpspredictor_table_filehandler:
        nrpspredictor_table = [l.strip().split('\t') for l in nrpspredictor_table_filehandler.readlines() if l.strip()]
        for cluster_id, cluster in _groupby(nrpspredictor_table[1:], lambda row: row[0]):
            fasta_data = utils.split_fasta_header(cluster_id)
            cluster_name = fasta_data['cluster_name']
            orf_id = fasta_data['orf_id']
            domain_id = fasta_data['domain_id']

            header = nrpspredictor_table[0]
            cluster_table = list(cluster)
            cluster_table.insert(0, header)
            for orf in details_data[cluster_name]['orfs']:
                if utils.sanitaze(orf['id']) == orf_id:
                    orf['domains'][int(domain_id)]['nrpspredictor_details'] = cluster_table
                    break


def extract_predictions_summary(run_id, geneclusters, details_data, report_directory):
    summary_path = os.path.join(report_directory, 'predictions_summary.tsv')
    with open(summary_path, 'w') as summary_handler:
        header = []
        for cluster_name, cluster in details_data.items():
            for orf in cluster['orfs']:
                domain_func_idx = collections.Counter()
                for domain_idx, domain in enumerate(orf['domains']):
                    function = utils.domain_function(domain['type'])
                    if function.upper() != 'A':
                        continue
                    seq_name = utils.retrieve_sequence_name(cluster_name, geneclusters)
                    opt_params = {'seq_name': seq_name, 'run_id': run_id, 'type': domain['type']}
                    fasta_id = utils.build_fasta_header(cluster_name,
                                                        orf['id'],
                                                        domain_idx,
                                                        function, domain_func_idx[function],
                                                        opt_params)

                    domain_func_idx[function] += 1
                    predictions = domain['predictions']
                    if not header:
                        header = [values[0] for values in predictions]
                        header.insert(0, '')
                        print('\t'.join(header), file=summary_handler)
                    row = [values[1] for values in predictions]
                    row.insert(0, fasta_id)
                    print('\t'.join(row), file=summary_handler)


def _groupby(table, operator):
    return itertools.groupby(sorted(table, key=operator), operator)


if __name__ == '__main__':
    logging.basicConfig(level=logging.INFO)
    logger = logging.getLogger(__name__)

    # parse command line
    parser = argparse.ArgumentParser()

    parser.add_argument(
        '--report_dir',
        help='Report directory',
        required=True,
    )

    args = parser.parse_args()
    logger.info('Start the building process')
    build_report(args.report_dir)
    logger.info('Building process completed')
