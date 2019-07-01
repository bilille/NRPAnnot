#!/usr/bin/env python3

import sys
import logging
import json
import re

logger = logging.getLogger()


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


def domain_function(domain_type):
    type2label = {
        "AMP-binding": "A",
        "AOX": "A",

        "PCP": "",
        "ACP": "",
        "NRPS-COM_Nterm": "",
        "NRPS-COM_Cterm":  "",
        "PKS_Docking_Nterm": "",
        "PKS_Docking_Cterm": "",
        "Trans-AT_docking": "",
        "Aminotran_1_2": "",
        "Aminotran_3": "",
        "Aminotran_4": "",
        "Aminotran_5": "",
        "Polyketide_cyc2": "",

        "Cglyc": "C",
        "CXglyc": "C",
        "Condensation_DCL": "C",
        "Condensation_LCL": "C",
        "Condensation_Starter": "C",
        "Condensation_Dual": "C",
        "Heterocyclization": "C",

        "Epimerization": "E",

        "Thioesterase": "TE",

        "PKS_KS": "KS",
        
        "PKS_AT": "AT",

        "PKS_KR": "KR",

        "PKS_DH": "DH",
        "PKS_DH2": "DH",

        "PKS_DHt": "DHt",
        "PKS_ER": "ER"
    }
    return type2label.get(domain_type, domain_type.split('_')[0])


def sanitaze(field):
    """
    Replace all non alphanumeric characters
    """
    return re.sub('[^0-9a-zA-Z]+', '-', str(field))
