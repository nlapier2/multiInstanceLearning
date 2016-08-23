# Author: Nathan LaPierre
# Date: August 21, 2016

import argparse
import sys
import time


def parseargs():    # handle user arguments
    parser = argparse.ArgumentParser(description='Generates area under curve based on predicted vs actual labels')
    parser.add_argument('--fasta', default='NONE', help='Location of fasta file.')
    parser.add_argument('--output', default='seeds.out', help='Where to write output.')
    parser.add_argument('--uc', default='NONE', help='Location of uclust file.')
    parser.add_argument('-v', '--verbose', default=False, action='store_true', help='Verbose output')
    args = parser.parse_args()
    return args


def main():
    start = time.time()
    args = parseargs()
    if args.fasta == 'NONE' or args.uc == 'NONE':
        print 'UCLUST and FASTA files must be specified with --uc and --fasta, respectively. Aborting...'
        sys.exit()

    uc = open(args.uc, 'r')
    fasta = open(args.fasta, 'r')
    out = open(args.output, 'w')
    for line in uc:
        if not line.startswith('#'):
            if line.startswith('S'):
                out.write(line)
                out.write(fasta.next())
                out.write(fasta.next())
            else:
                fasta.next()
                fasta.next()
    uc.close()
    fasta.close()
    out.close()

    if args.verbose:
        print "\nProgram execution time: " + str(time.time() - start) + " seconds"


if __name__ == "__main__":
    main()