#!/usr/bin/env python3

import sys
import logging
import argparse
import json
import re


def geneclusters2json(geneclusters_js_file):
    """
    Extract geneclusters & details_data as json data
    Return geneclusters_json & details_data_json
    """
    with open(geneclusters_js_file, 'r') as js_handler:
        lines = [l.strip() for l in js_handler.readlines()]
        details_start_idx = 0
        geneclusters_pattern = r'var\s*geneclusters\s*=\s*'
        details_data_pattern = r'var\s*details_data\s*=\s*'
        if not re.match(geneclusters_pattern, lines[0]):
            logger.fatal('The genecluster.js file need to declare geneclusters variable')
            sys.exit('genecluster variable does not exist')

        for idx, line in enumerate(lines):
            if re.match(details_data_pattern, line):
                details_start_idx = idx
                break
        if not details_start_idx:
            logger.fatal('The genecluster.js file need to declare details_data variable')
            sys.exit('details_data variable does not exist')

        geneclusters_json = lines[:details_start_idx]
        details_data_json = lines[details_start_idx:]

        # remove javascript variable assignment and final ';'
        geneclusters_json[0] = re.sub(geneclusters_pattern, '', geneclusters_json[0])
        details_data_json[0] = re.sub(details_data_pattern, '', details_data_json[0])
        geneclusters_json = ''.join(geneclusters_json)
        details_data_json = ''.join(details_data_json)
        geneclusters_json = geneclusters_json.strip(';')
        details_data_json = details_data_json.strip(';')

        # decode json string
        try:
            geneclusters_json = json.loads(geneclusters_json)
        except json.JSONDecodeError as e:
            logger.fatal("Can't decode geneclusters json")
            logger.exception(e)
            sys.exit('JSON not well formed error')
        try:
            details_data_json = json.loads(details_data_json)
        except json.JSONDecodeError as e:
            logger.fatal("Can't decode details_data json")
            logger.exception(e)
            sys.exit('JSON not well formed error')
        return geneclusters_json, details_data_json


if __name__ == '__main__':
    logging.basicConfig(level=logging.INFO)
    logger = logging.getLogger()

    # parse command line
    parser = argparse.ArgumentParser()

    parser.add_argument(
        '--run_id',
        help='antiSMASH run ID',
        required=True,
        type=str
    )
    parser.add_argument(
        '--geneclusters',
        help='Geneclusters JS file exported by antiSMASH',
        required=True,
        type=argparse.FileType('r')
    )
    parser.add_argument(
        '--json_out',
        help='JSON ouput file',
        required=True,
        type=argparse.FileType('w')
    )

    args = parser.parse_args()

    geneclusters_json, details_data_json = geneclusters2json(args.geneclusters.name)
    json_data = {'run_id': args.run_id,
                 'geneclusters': geneclusters_json,
                 'details_data': details_data_json}
    print(json.dumps(json_data, indent=4), file=args.json_out)
