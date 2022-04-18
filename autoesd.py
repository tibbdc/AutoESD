# -*-coding:utf-8 -*-

import os
import sys

import argparse

def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument(
        '--technology', '-t', help='choose the target technology pipeline', choices=["PSS","PDS","FDS","FDD","ODD"], required=True)
    parser.add_argument('--input', '-i', help='input sequences file', required=True)
    parser.add_argument('--ref', '-r', help='reference genome', required=True)
    parser.add_argument('--dir', '-d', help='outputdir', required=True)
    parser.add_argument('--conf', '-c', help='config', required=True)
    parser.add_argument('--plasmid', '-p', help='plasmid file, must be provided in PSS, PDS.', required=False)
    parser.add_argument('--screeningmarker', '-s', help='screening marker seq file, must be provided in PDS, FDS, FDD, ODD.', required=False)
    parser.add_argument('--ustp', '-u', help='user specific test primer', required=False)
    parser.add_argument('--smr', '-m', help='screening maker removal', choices=["No","Yes"],required=False)

    parser.add_argument(
        '--verb', '-v', help='shows program progress', default=False, required=False)

    return parser.parse_args()

def main():
    # read args
    args = parse_args()
    technology = args.technology
    inputfile = args.input
    ref = args.ref
    outputdir = args.dir
    config = args.conf

    plasmidfile = args.plasmid
    if technology in ["PSS","PDS"] and plasmidfile == None :
        raise ValueError('plasmid file, must be provided in PSS, PDS.')
    
    screeningmarker = args.screeningmarker
    if technology in ["PDS","FDS","FDD","ODD"] and screeningmarker == None:
        raise ValueError('screening marker seq file, must be provided in PDS, FDS, FDD, ODD')

    user_specific_test_primer = args.ustp
    screening_maker_removal = args.smr

    if technology == "PSS":
        cmd = "python ./PSS/plasmid_single_single.py -i %s  -p %s -r %s  -c %s -d %s" % (inputfile,plasmidfile,ref,config,outputdir)
    elif technology == "PDS":
        cmd = "python ./PDS/plasmid_double_single.py -i %s -s %s -r %s -p % -c %s -d %s" % (inputfile,screeningmarker,ref,plasmidfile,config,outputdir)
        if user_specific_test_primer:
            cmd += " -u " + user_specific_test_primer
        if screening_maker_removal:
            cmd += " -m " + screening_maker_removal
    elif technology == "FDS":
        cmd = "python ./FDS/fragment_double_single.py -i %s  -s %s -r %s  -c %s -d %s " % (inputfile,screeningmarker,ref,config,outputdir)
        if screening_maker_removal:
            cmd += " -m " + screening_maker_removal 
    elif technology == "FDD":
        cmd = "python ./FDD/fragment_double_double.py -i %s  -s %s -r %s  -c %s -d %s " % (inputfile,screeningmarker,ref,config,outputdir)
        if user_specific_test_primer:
            cmd += " -u " + user_specific_test_primer
        if screening_maker_removal:
            cmd += " -m " + screening_maker_removal
    elif technology == "ODD":
        cmd = "python ./ODD/oligonucleotide_double_double.py -i %s  -s %s -r %s  -c %s -d %s " % (inputfile,screeningmarker,ref,config,outputdir)
        if user_specific_test_primer:
            cmd += " -u " + user_specific_test_primer
        if screening_maker_removal:
            cmd += " -m " + screening_maker_removal

    os.system(cmd)
if __name__ == "__main__":
    main()