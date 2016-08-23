# Author: Nathan LaPierre
# Date: July 25, 2016

import argparse
import os
import sys
import time


def parseargs():    # handle user arguments
    parser = argparse.ArgumentParser(description='Limit the number of reads per cluster in patient files.')
    parser.add_argument('--fasta', default='NONE', help='Location of input FASTA file. Required.')
    parser.add_argument('--uclust', default='NONE', help='Location of input UCLUST file. Required.')
    parser.add_argument('--limit', default=-1, help='Limit the number of selections per cluster. Default: 1.')
    parser.add_argument('--output', default='limited', help='Files to write lines to. Default: limited{.fasta/.uc}')
    parser.add_argument('-v', '--verbose', default=False, action='store_true', help='Verbose output')
    args = parser.parse_args()
    return args


def main():
    start = time.time()
    args = parseargs()
    if args.fasta == 'NONE' or args.uclust == 'NONE':
        print "Location of input FASTA and UCLUST files with data must be specified with --fasta and --uclust."
        sys.exit()
    try:
        limit = int(args.limit)
    except ValueError:
        print "Error: Limit must be an integer greater than 0."
        sys.exit()
    if limit < 1:
        print "Error: Limit must be an integer greater than 0."
        sys.exit()

    outfasta = open(args.output + '.fasta', 'w')
    outuclust = open(args.output + '.uc', 'w')
    cluster_counts = {}     # keeps track of counts of reads from each cluster for each patient
    try:
        infasta = open(args.fasta, 'r')
        inuclust = open(args.uclust, 'r')
        for line in inuclust:
            if line.startswith('H') or line.startswith('S'):
                fastaline = infasta.next()
                fields = line.split('\t')       # fields in clustering output lines
                tag = fields[8].split(' ')[0]   # gets the tag that matches a patient
                if tag not in cluster_counts:
                    cluster_counts[tag] = [0]
                counts = cluster_counts[tag]
                clustnum = int(fields[1])    # cluster number
                while len(counts) < clustnum + 1:
                    counts.append(0)
                if counts[clustnum] < limit:
                    counts[clustnum] += 1
                    outfasta.write(fastaline)
                    fastaline = infasta.next()
                    outfasta.write(fastaline)
                    outuclust.write(line)
                else:
                    infasta.next()
        infasta.close()
        inuclust.close()
        outfasta.close()
        outuclust.close()
    except:
        outfasta.close()
        os.remove(args.output + '.fasta')   # remove our temporary file
        outuclust.close()
        os.remove(args.output + '.uc')   # remove our temporary file
        print "Unexpected error: ", sys.exc_info()[0]
        raise

    if args.verbose:
        print "Program execution time: " + str(time.time() - start) + " seconds"


if __name__ == "__main__":
    main()