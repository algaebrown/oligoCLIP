#!/usr/bin/env python

# transitionning to python2/python3 support
# uncomment from this compatibility import list, as py3/py2 support progresses
from __future__ import print_function
from __future__ import division
# from __future__  import absolute_import
# from __future__  import unicode_literals
# from future import standard_library
# from future.builtins import builtins
# from future.builtins import utils
# from future.utils import raise_with_traceback
# from future.utils import iteritems

import os
import argparse



###############################################################################
def parse_cutadapt_file(report, paired_end):
    #############################
    #print("parse_cutadapt_file report:", report)
    if os.path.getsize(report) == 0:
        return
    old_cutadapt = get_cutadapt_version(report) <= 8
    if paired_end:
        if old_cutadapt:
            return parse_old_cutadapt_file_pe(report)
        else:
            return parse_new_cutadapt_file(report)
    else:
        return parse_new_cutadapt_file(report)

###############################################################################
def get_cutadapt_version(report):
    ##############################
    with open(report) as file_handle:
            version = next(file_handle)
    try:
        version = version.split()[-4]
    except:
        1

    return int(version.split(".")[1])


###############################################################################
def parse_old_cutadapt_file_pe(report):
    ################################
    report_dir = {}
    try:
        with open(report) as report:
            next(report) #header
            next(report) #paramaters
            next(report) #max error rate
            next(report) #adapters (known already)
            processed_reads = [x.strip() for x in next(report).strip().split(":")]
            processed_bases = [x.strip() for x in next(report).strip().split(":")]
            trimmed_reads   = [x.strip() for x in next(report).strip().split(":")]
            quality_trimmed = [x.strip() for x in next(report).strip().split(":")]
            trimmed_bases   = [x.strip() for x in next(report).strip().split(":")]
            too_short_reads = [x.strip() for x in next(report).strip().split(":")]
            too_long_reads  = [x.strip() for x in next(report).strip().split(":")]
            total_time      = [x.strip() for x in next(report).strip().split(":")]
            time_pre_read   = [x.strip() for x in next(report).strip().split(":")]
            report_dir[processed_reads[0]] = int(processed_reads[1])
            report_dir[processed_bases[0]] = int(processed_bases[1].split()[0])
            report_dir[trimmed_reads[0]] = int(trimmed_reads[1].split()[0])
            report_dir[quality_trimmed[0]] = int(quality_trimmed[1].split()[0])
            report_dir[trimmed_bases[0]] = int(trimmed_bases[1].split()[0])
            report_dir[too_short_reads[0]] = int(too_short_reads[1].split()[0])
            report_dir[too_long_reads[0]] = int(too_long_reads[1].split()[0])
            report_dir[trimmed_bases[0]] = int(trimmed_bases[1].split()[0])
    except:
            print(report)
    return report_dir


###############################################################################

def parse_new_cutadapt_file(report):
    ################################
    report_dict = {}
    try:
        with open(report) as file_handle:
            remove_header(file_handle)
            processed_reads = get_number(next(file_handle))
            paired_file = processed_reads[0] == 'Total read pairs processed'
            if paired_file:
                r1_adapter = get_number_and_percent(next(file_handle))
                r2_adapter = get_number_and_percent(next(file_handle))
            else:
                adapter = get_number_and_percent(next(file_handle))

            too_short = get_number_and_percent(next(file_handle))
            written = get_number_and_percent(next(file_handle))
            next(file_handle)
            bp_processed = get_number(strip_bp(next(file_handle)))
            if paired_file:
                r1_bp_processed = get_number(strip_bp(next(file_handle)))
                r2_bp_processed = get_number(strip_bp(next(file_handle)))

            bp_quality_trimmed = get_number_and_percent(strip_bp(next(file_handle)))
            if paired_file:
                r1_bp_trimmed = get_number(strip_bp(next(file_handle)))
                r2_bp_trimmed = get_number(strip_bp(next(file_handle)))

            bp_written = get_number_and_percent(strip_bp(next(file_handle)))
            if paired_file:
                r1_bp_written = get_number(strip_bp(next(file_handle)))
                r2_bp_written = get_number(strip_bp(next(file_handle)))

    except Exception as e:
        print("exception occurred in cutadapt parsing: {}".format(e))
        print(report)
        return report_dict

    report_dict['Processed reads'] = processed_reads[1]
    if paired_file:
        report_dict["Read 1 with adapter"] = r1_adapter[1]
        report_dict["Read 1 with adapter percent"] = r1_adapter[2]
        report_dict["Read 2 with adapter"] = r2_adapter[1]
        report_dict["Read 2 with adapter percent"] = r2_adapter[2]
        report_dict['Read 1 basepairs processed'] = r1_bp_processed[1]
        report_dict['Read 2 basepairs processed'] = r2_bp_processed[1]
        report_dict['Read 1 Trimmed bases'] = r1_bp_trimmed[1]
        report_dict['Read 2 Trimmed bases'] = r2_bp_trimmed[1]
        report_dict['Read 1 {}'.format(bp_written[0])] = r1_bp_written[1]
        report_dict['Read 2 {}'.format(bp_written[0])] = r2_bp_written[1]
    else:
        report_dict['Reads with adapter'] = adapter[1]
        report_dict['Reads with adapter percent'] = adapter[2]


    report_dict['Too short reads'] = too_short[1]
    report_dict['Reads that were too short percent'] = too_short[2]
    report_dict['Reads Written'] = written[1]
    report_dict['Reads Written percent'] = written[2]
    report_dict['Processed bases'] = bp_processed[1]
    report_dict['Trimmed bases'] = bp_quality_trimmed[1]
    report_dict['Trimmed bases percent'] = bp_quality_trimmed[2]
    report_dict[bp_written[0]] = bp_written[1]
    report_dict["{} percent".format(bp_written[0])] = bp_written[2]
    return report_dict


###############################################################################
## parse_new_cutadapt_file utilities
######################################
def get_number_and_percent(line):
    """
    Parses cutadapt line containing a number (string) and returns
    number typecasted to int, as well as a percentage (float), as a list.

    :param line: basestring
    :return line: list
    """
    line = [x.strip() for x in line.strip().split(":")]

    line = [line[0]] + line[1].split()
    line[2] = float(line[2][1:-2])
    line[1] = int(line[1].replace(",", ""))
    return line


def get_number(line):
    """
    Parses cutadapt line containing a number (string) and returns
    number typecasted to int.

    :param line: basestring
    :return number: int
    """
    line = [x.strip() for x in line.strip().split(":")]
    line[1] = int(line[1].replace(",", ""))
    return line


def strip_bp(line):
    return line.replace("bp", "")


def remove_header(file_handle):
    """ for both SE and PE output removes header unifromly from cutadapt metrics file"""
    while not next(file_handle).startswith('=== Summary ==='):  # skip all lines up until it hits this Summary line
        continue
    next(file_handle)  # blank line after === Summary === line

###############################################################################


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description="prints a report_dict for a cutadapt report")
    parser.add_argument("--report", help="report", required=True)
    args = parser.parse_args()
    print("args:", args)
    report_dict = parse_cutadapt_file(args.report)
    print("report_dict:", report_dict)