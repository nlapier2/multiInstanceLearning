# Author: Nathan LaPierre
# Date: March 12, 2016

import argparse
import shlex
import subprocess
import sys


def parseargs():    # handle user arguments
    parser = argparse.ArgumentParser(description='Bulk NCBI download script')
    parser.add_argument('--start', default='NONE', help='NCBI SRA run number to start on. Required.')
    parser.add_argument('--end', default='NONE', help='NCBI SRA run number to end on. Default: only download start run')
    parser.add_argument('--dlpath', default='NONE', help='Path to fastq-dump. Default: do not perform download.')
    parser.add_argument('--out', default='./', help='Output directory. Default is current directory.')
    parser.add_argument('--dir', default='./', help='Directory of files for fastq_to_fasta and UCLUST options.')
    parser.add_argument('--fastapath', default='NONE', help='Path to fastq_to_fasta tool. Default: do not perform this')
    parser.add_argument('--uclustpath', default='NONE', help='Path to UCLUST executable. Default: do not perform this')
    parser.add_argument('-t', '--test', default=False, action='store_true', help='Do a test run.')
    parser.add_argument('-v', '--verbose', default=False, action='store_true', help='Verbose output')
    args = parser.parse_args()
    return args


def main():
    args = parseargs()
    if args.test:
        if args.dlpath == 'NONE':
            print "Error: to test download, specify location of fastq-dump with --dlpath"
            sys.exit()
        subprocess.call(shlex.split(args.dlpath + ' -X 3 -Z SRR390728'))
        sys.exit()
    if args.start == 'NONE':
        print "Error: SRA run number to start on must be specified with --start"
        sys.exit()
    if args.verbose:
        print "Verbose output requested."
    start = int(args.start.split('SRR')[1])
    if args.end == 'NONE':
        end = start
    else:
        end = int(args.end.split('SRR')[1])

    if args.dlpath != 'NONE':
        for cur in range(start, end + 1):
            if args.verbose:
                print "Calling: " + args.dlpath + ' -O ' + args.out + ' SRR' + str(cur)
            subprocess.call(shlex.split(args.dlpath + ' -O ' + args.out + ' SRR' + str(cur)))
    if args.fastapath != 'NONE':
        for cur in range(start, end + 1):
            if args.verbose:
                print "Calling: " + args.fastapath + ' -Q33 -i ' + args.dir + 'SRR' + str(cur) \
                      + '.fastq -o ' + args.out + 'SRR' + str(cur) + '.fasta'
            subprocess.call(shlex.split(args.fastapath + ' -Q33 -i ' + args.dir + 'SRR' + str(cur)
                                        + '.fastq -o ' + args.out + 'SRR' + str(cur) + '.fasta'))
    if args.uclustpath != 'NONE':
        for cur in range(start, end + 1):
            if args.verbose:
                print "Calling: " + args.uclustpath + ' --input ' + args.dir + 'SRR' + str(cur) \
                      + '.fasta --uc ' + args.out + 'SRR' + str(cur) + '.uc --usersort --id 0.40'
            subprocess.call(shlex.split(args.uclustpath + ' --input ' + args.dir + 'SRR' + str(cur)
                                        + '.fasta --uc ' + args.out + 'SRR' + str(cur) + '.uc --usersort --id 0.40'))


if __name__ == "__main__":
    main()
