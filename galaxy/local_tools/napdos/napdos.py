#!/usr/bin/env python3

import os
import logging
import argparse
import requests
import bs4
import itertools
import time
import utils


def _groupby(table, operator):
    return itertools.groupby(sorted(table, key=operator), operator)


def run_napdos(fasta, output_dir,
               domain_type='C',
               query_type='aa',
               ks_hmm_evalue=1e-5, ks_min_matchlength=200,
               c_hmm_evalue=1e-5, c_min_matchlength=200,
               path_blast_evalue=1e-5, max_hits=1):
    """
    As we can't perform a local installation of napdos.
    Send a request to napdos server
    To analyse our sequence with napdos:
    1_ Send a request to process_request.cgi
    2_ Send a request to domain_blast2.cgi (actually run the job)
    """
    if not os.path.isdir(output_dir):
        logger.info("Create output directory")
        try:
            os.makedirs(output_dir)
        except OSError as e:
            logger.fatal("Can't create output directory: %s" % output_dir)
            logger.exception(e)

    process_request_payload = {
        'config_filename': 'pks_03_sdsc.cfg',
        'query_type': query_type,
        'c_hmm_evalue': c_hmm_evalue,
        'c_min_matchlength': c_min_matchlength,
        'ks_hmm_evalue': ks_hmm_evalue,
        'ks_min_matchlength': ks_min_matchlength,
        'path_blast_evalue': path_blast_evalue,
        'max_hits': max_hits,
        'Sequence': ''  # .join(open(fasta, 'r').readlines()),
        # 'seqfile': ''
    }

    if domain_type == 'C':
        process_request_payload['ref_seq_file'] = 'all_C_public_12062011.faa'
    else:
        process_request_payload['ref_seq_file'] = 'all_KS_public_12062011.faa'

    process_request_url = 'http://napdos.ucsd.edu/cgi-bin/process_request.cgi'
    logger.info('Send the request to %s', process_request_url)
    logger.debug('Parameters:%s', process_request_payload)
    process_request_response = requests.post(process_request_url,
                                             data=process_request_payload,
                                             files={'seqfile': open(fasta, 'rb')})
    assert_response(process_request_response)

    # At least one minute between requests is required
    if 'Time since previous user query' in process_request_response.text:
        logger.fatal('Please allow at least one minute between requests')
        exit('Wait a minute')

    # retrieve parameters to submit job
    process_request_soup = bs4.BeautifulSoup(process_request_response.content, 'html.parser')
    submit_job_payload = {input_tag['name']: input_tag['value'] for input_tag in process_request_soup.find_all('input')}
    submit_job_url = 'http://napdos.ucsd.edu/cgi-bin/domain_blast2.cgi'
    logger.info('Send the request to %s', submit_job_url)
    logger.debug('Parameters: %s', submit_job_payload)
    submit_job_response = requests.post(submit_job_url, data=submit_job_payload)
    assert_response(submit_job_response)

    # manage the case where no matches occured or there is an input problem
    submit_job_soup = bs4.BeautifulSoup(submit_job_response.content, 'html.parser')
    if 'Sorry, your request could not be processed' in submit_job_response.text:
        error_message = 'Sorry, your request could not be processed'
        div_soup = submit_job_soup.find(id="maininner")
        if div_soup:
            error_message = '\n'.join(div_soup.stripped_strings)
        logger.fatal(error_message)
        exit("Can't processed the request")
    elif 'No matches found' in submit_job_response.text:
        logger.info('No matches found. You may want to try again using a higher e-value cutoff.')
        return

    logger.info('Download NaPDoS table')
    save_napdos_table(submit_job_soup, output_dir)

    logger.info('Download NaPDoS tree')
    save_napdos_tree(submit_job_soup, output_dir)


def assert_response(response):
    """
    assert that a success code is return
    by the server
    """
    try:
        response.raise_for_status()
    except requests.HTTPError:
        logger.fatal('Received error %s', response.status_code)
        logger.debug('Content:%s', response.text)
        exit('HTTP Error')


def save_napdos_table(results_soup, output_dir):
    def download_link(tag):
        return tag.name == 'a' and \
            tag.has_attr('href') and \
            tag.text == 'DOWNLOAD'

    download_soup = results_soup.find_all(download_link)
    if not download_soup:
        logger.fatal('The link to download the results as tab-delimited file is missing')
        exit('Missing download link')

    download_response = requests.get(download_soup[0]['href'])
    assert_response(download_response)

    table_file = os.path.join(output_dir, 'table.tab')
    logger.debug('Destination:%s', table_file)
    with open(table_file, 'wb') as table_handler:
        table_handler.write(download_response.content)


def save_napdos_tree(results_soup, output_dir):
    """
    Extract needed info from NapDoS response to
    actually download the tree with request_napdos_tree
    """
    def checkbox(tag):
        return tag.name == 'input' and \
            (tag.has_attr('name') and tag['name'] == 'list') and \
            (tag.has_attr('type') and tag['type'] == 'checkbox')

    def align_group(tag):
        return tag.name == 'input' and \
            (tag.has_attr('type') and tag['type'] == 'hidden') and \
            (tag.has_attr('name') and tag['name'] == 'align_group')

    def job_id(tag):
        return tag.name == 'input' and \
            (tag.has_attr('type') and tag['type'] == 'hidden') and \
            (tag.has_attr('name') and tag['name'] == 'job_id')

    # Extract the info needed to request svg tree
    # using the previous filter functions
    checkboxes_soup = results_soup.find_all(checkbox)
    checkboxes_val = [c_soup['value'] for c_soup in checkboxes_soup]

    align_group_soup = results_soup.find_all(align_group)
    align_group_val = align_group_soup[0]['value']

    job_id_soup = results_soup.find_all(job_id)
    job_id_val = job_id_soup[0]['value']

    job_id_soup = results_soup.find_all(job_id)
    job_id_val = job_id_soup[0]['value']

    # Determine if the sequences come from antiSMASH run. If so, treat
    # the results by cluster (one svg file per cluster)
    if 'from-geneclusters' in checkboxes_val[0]:
        # antiSMASH
        for cluster_key, cluster in _groupby(checkboxes_val,
                                             lambda v: utils.split_fasta_header(v)['cluster_name']):
            logger.info('Current cluster:%s', cluster_key)
            response = request_napdos_tree(cluster, align_group_val, job_id_val)
            tree_file = os.path.join(output_dir, '%s.svg' % cluster_key)
            logger.debug('Destination:%s', tree_file)
            with open(tree_file, 'wb') as tree_handler:
                tree_handler.write(response.content)
            time.sleep(1)  # do not flood
    else:
        # otherwise, download all in the same tree
        response = request_napdos_tree(checkboxes_val, align_group_val, job_id_val)
        tree_file = os.path.join(output_dir, 'results.svg')
        logger.debug('Destination:%s', tree_file)
        with open(tree_file, 'wb') as tree_handler:
            tree_handler.write(response.content)


def request_napdos_tree(seq_list, align_group, job_id):
    """
    To export svg tree, we need to select candidates sequences
    of interest and send a POST request.
    Then, in the response, a link will be available to
    actually download the tree
    """
    def download_link(tag):
        return tag.name == 'a' and \
            (tag.has_attr('href') and 'svg' in tag['href']) and \
            tag.text.strip() == 'DOWNLOAD'

    request_payload = {
        'list': seq_list,
        'job_id': job_id,
        'align_group': align_group,
        'result_type': 'newick'
        }

    request_url = 'http://npdomainseeker.sdsc.edu/cgi-bin/align_seqs2.cgi'
    logger.info('Send the request to %s', request_url)
    logger.debug('Parameters:%s', request_payload)
    request_response = requests.post(request_url,
                                     data=request_payload)
    assert_response(request_response)

    results_soup = bs4.BeautifulSoup(request_response.content, 'html.parser')
    download_soup = results_soup.find_all(download_link)

    if not download_soup:
        logger.fatal('The link to download the results as svg file is missing')
        exit('Missing download link')

    download_response = requests.get(download_soup[0]['href'])
    assert_response(download_response)
    return download_response

# def extract_table(table_soup):
#     """
#     Transform html_table to tabular table
#     """

#     table = []
#     header = ('query_id', 'db_math_id', 'percent_identity', 'align_length',
#               'e-value', 'pathway_product', 'domain_class')
#     table.append('\t'.join(header))
#     # retrieve all lines
#     for line_soup in table_soup.find_all('tr'):
#         if 'Query id' in str(line_soup):
#             continue  # ignore header

#         val_list = []
#         for td_soup in line_soup.find_all('td'):
#             # ignore checkbox
#             if not td_soup.contents or td_soup.find('input'):
#                 continue
#             link_soup = td_soup.find('a')
#             if link_soup:
#                 val_list.append(str(link_soup.contents[0]))
#             else:
#                 val_list.append(str(td_soup.contents[0]))
#         table.append('\t'.join(val_list))

#     return '\n'.join(table)


if __name__ == '__main__':
    logging.basicConfig(level=logging.INFO)
    logger = logging.getLogger()

    # parse command line
    parser = argparse.ArgumentParser()

    parser.add_argument(
        '--fasta',
        help='Sequence(s) in FASTA format. Please avoid using sequence'
        ' id numbers containing non-standard characters'
        ' (e.g. commas, slashes, colons, parentheses, ampersands, etc),'
        ' as well as ambiguity codes'
        ' "R" and "Y" for nucleic acid sequence files.',
        required=True,
        type=argparse.FileType('r')
    )

    parser.add_argument(
        '--output_dir',
        default="napdos"
    )

    parser.add_argument(
        '--domain_type',
        help='Domain type. Default is %(default)s',
        choices=['KS', 'C'],
        default='C'
    )

    parser.add_argument(
        '--query_type',
        help="Input type. Default is %(default)s. "
        ' aa - Predicted protein sequences (amino acid).'
        ' cds - Predicted coding sequences or PCR products (DNA).'
        ' genome - Genome or metagenome contigs (DNA).',
        choices=['aa', 'cds', 'genome'],
        default='aa'
    )

    # KS domain detection
    parser.add_argument(
        '--ks_hmm_evalue',
        help='KS domain detection - HMM e-value cutoff',
        choices=['1e-5', '1e-3', '1e-1', '1', '10', '1e-50'],
        default='1e-5',
        )
    parser.add_argument(
        '--ks_min_matchlength',
        help='KS domain detection - minimum match length',
        choices=['200', '150', '100'],
        default='200'
        )

    # C domain detection
    parser.add_argument(
        '--c_hmm_evalue',
        help='C domain detection - HMM e-value cutoff',
        choices=['1e-5', '1e-3', '1e-1', '1', '10', '1e-50'],
        default='1e-5',
        )
    parser.add_argument(
        '--c_min_matchlength',
        help='C domain detection - minimum match length',
        choices=['200', '400', '300', '100'],
        default='200'
        )

    # Pathway assignment
    parser.add_argument(
        '--path_blast_evalue',
        help='Pathway assignment - BLASTP e-value cutoff',
        choices=['1e-5', '1e-20', '1e-10', '1e-1', '1'],
        default='1e-5',
        )
    parser.add_argument(
        '--max_hits',
        help='Pathway assignment - max db hits per domain ',
        choices=['1', '3', '5', '10'],
        default='1'
        )

    args = parser.parse_args()
    fasta_file = args.fasta.name
    output_dir = args.output_dir

    params = vars(args)
    params.pop('fasta', None)
    params.pop('output_dir', None)

    logger.info('Start requests to NaPDoS')
    run_napdos(fasta_file, output_dir, **params)
    logger.info('Requests to NapDoS finised')
