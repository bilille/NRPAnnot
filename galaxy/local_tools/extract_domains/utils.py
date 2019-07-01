#!/usr/bin/env python3

import re

SEP = '-'


def retrieve_sequence_name(cluster_name, geneclusters):
    """
    Retrieve the sequence name of the sequence used to predict
    the cluster
    As this field is optional, return undef if we can't retrieve it
    """
    seq_pattern = r"id=(.+?)&"
    cluster = geneclusters.get(cluster_name, {})
    orfs = cluster.get('orfs', [])
    first_orf = {}
    try:
        first_orf = orfs[0]
    except IndexError:
        pass
    description = first_orf.get('description', '')
    m = re.search(seq_pattern, description)
    if m:
        return m.group(1)
    return 'seqname-undef'


def domain_function(domain_type):
    """
    Map domain type to a function
    This function (domain_function) have a JS counterpart in
    build_report/js/draw_report.js -- domain_label_type.
    If you perform any changes, update the counterpart.
    """
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
    This function have a JS counterpart in build_report/js/draw_report.js.
    If you perform any changes, update the counterpart.
    """
    return re.sub('[^0-9a-zA-Z]+', SEP, str(field))


def build_fasta_header(cluster_name, orf_id, domain_idx, func, func_idx, opt_params):
    """
    No underscores are permitted within each field of the header.
    Use sanitaze function to this purpose.
    This function have a JS counterpart in build_report/js/draw_report.js -- complete_with_header.
    If you perform any changes on this function, you have to update accordingly
    the split_fasta_header function and the JS counterpart.

    The from-geneclusters flag is injected for downstream napdos tool:
    The napdos results will be treated by cluster. See the tool for more details.
    """
    opt_str = dict2str(opt_params)
    return '{cluster_name}_{orf_id}_{domain_id}_{func}{func_idx}_from-geneclusters {opt}'.format(
        cluster_name=sanitaze(cluster_name),
        orf_id=sanitaze(orf_id),
        domain_id=sanitaze(domain_idx + 1),  # 1-based index
        func=sanitaze(func),
        func_idx=func_idx + 1,   # 1-based index
        opt=opt_str
    )


def split_fasta_header(header):
    """
    This funciton is intented to retrieve data from headers with the build_fasta_header
    function.
    This function is shared between extract_domains/napdos/build_report modules.
    Pay more attention when updating this function as it can inpact 3 modules;
    """

    fields = {}
    cluster_name, orf_id, domain_id, _ = header.split('_', maxsplit=3)
    fields['cluster_name'] = cluster_name
    fields['orf_id'] = orf_id
    fields['domain_id'] = int(domain_id) - 1  # 0-based index
    return fields


def dict2str(d):
    values = []
    for k in sorted(d.keys()):
        values.append('%s=%s' % (k, sanitaze(d[k])))
    return ';'.join(values)
