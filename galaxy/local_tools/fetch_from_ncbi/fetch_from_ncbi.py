#!/usr/bin/env python

import os
import requests
from httplib import IncompleteRead
import argparse
import logging


def download_from_ncbi(accession, workdir, sequence_type):

    NCBI_URL = 'https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi'
    error_patterns = (
        'Error reading from remote server',
        'Bad gateway',
        'Cannot process ID list',
        'server is temporarily unable to service your request',
        'Service unavailable',
        'Server Error',
        'ID list is empty',
        'Resource temporarily unavailable',
    )
    params = {}

    # delete / characters and as NCBI ignores IDs after #, do the same.
    params['id'] = accession.replace('/', '').split('#', 1)[0]

    if sequence_type == 'nucleic':
        params['db'] = 'nucleotide'
        params['rettype'] = 'gbwithparts'
        params['retmode'] = 'text'
    else:
        params['db'] = 'protein'
        params['rettype'] = 'fasta'

    logger.info('Downloading the input file from NCBI: %s', params['id'])

    try:
        r = requests.get(NCBI_URL, params=params, stream=True)
        logger.info('Url: %s', r.url)
    except (requests.exceptions.RequestException, IncompleteRead) as e:
        raise Exception((str(e), -1))

    if r.status_code != requests.codes.ok:
        raise Exception(("Failed to download file with id {} from NCBI".format(params['id']),
                         r.status_code))

    outfile_name = os.path.join(workdir, 'sequence')

    logger.info('Write file')
    with open(outfile_name, 'wb') as fh:
        first = True
        # use a chunk size of 4k, as that's what most filesystems use these days
        for chunk in r.iter_content(4096):
            if first:
                first = False
                for pattern in error_patterns:
                    if pattern in chunk:
                        raise Exception(("Failed to download file with id {} from NCBI: {}".format(
                            params['id'], pattern), 1))

            fh.write(chunk)
    logger.info('Downloading the input file from NCBI: Done')


if __name__ == '__main__':
    logging.basicConfig(level=logging.INFO)
    logger = logging.getLogger()

    # parse command line
    parser = argparse.ArgumentParser()

    parser.add_argument(
        '--accession',
        help='Accession to fetch from ncbi',
        required=True,
        type=str
    )
    parser.add_argument(
        '--sequence_type',
        action='store',
        choices=['nucleic', 'protein'],
        required=True,
        type=str
    )
    parser.add_argument(
        '--output',
        help='Output directory (must exists)',
        required=True,
        type=str
    )

    args = parser.parse_args()

    download_from_ncbi(args.accession, args.output, args.sequence_type)
