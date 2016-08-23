# Author: Nathan LaPierre
# Date: June 29, 2016

import argparse
import os
import shlex
import subprocess
import sys


def write_conf(config, path):
    conf = open(config, 'w')
    conf.write('max_rd_len=180\n')
    conf.write('[LIB]\n')
    conf.write('avg_ins=700\n')
    conf.write('reverse_seq=0\n')
    conf.write('asm_flags=3\n')
    conf.write('rank=1\n')
    conf.write('rd_len_cutoff=100\n')
    conf.write('\nq=' + path+'\n')
    conf.close()


def assemble(args, kmer):
    start = int(args.start.split('SRR')[1])
    if args.end == 'NONE':
        end = start
    else:
        end = int(args.end.split('SRR')[1])
    config = args.out + 'config-' + args.start + '-' + args.end + '.config'    # location of config file
    contigs = args.out + 'contigs-' + args.start + '-' + args.end + '.fasta'     # the file that will hold all contigs
    if args.verbose:
        print 'Attempting to open: ' + config + ' and ' + contigs
    conf = open(config, 'w')
    conf.close()

    contig = open(contigs, 'w')
    for cur in range(start, end + 1):
        directory = args.out + 'SRR' + str(cur) + '/'     # directory to output assembly of this file to
        path = args.path + 'SRR' + str(cur) + args.extension         # path to current fastx file
        if args.verbose:
            print 'Calling: mkdir ' + directory
        os.mkdir(directory)
        os.chdir(directory)
        if args.verbose:
            print 'Making config file for SRR' + str(cur)
        write_conf(config, path)
        if args.verbose:
            print 'Calling: ' + args.location + ' all -s ' + config + ' -K ' + str(kmer) + \
                  '-R -V -o graph_prefix 1>all.log 2>all.err'
        subprocess.call(shlex.split(args.location + ' all -s ' + config + ' -K ' + str(kmer) +
                                    '-R -V -o graph_prefix 1>all.log 2>all.err'))
        if args.verbose:
            print 'Dumping scaffolds into bulk contig file...'
        scaffolds = open(directory+'graph_prefix.scafSeq', 'r')
        for line in scaffolds:
            if line.startswith('>'):    # defline, not a read
                contig.write('>SRR' + str(cur) + '.0 ' + line.split('>')[1])
            elif len(line) > 3:     # non-empty line (contains read)
                contig.write(''.join(line.split()) + '\n')
    contig.close()
    os.remove(config)


def parseargs():    # handle user arguments
    parser = argparse.ArgumentParser(description='Script that performs assembly on a number of patient files.')
    parser.add_argument('--start', default='NONE', help='File number to start on. Required.')
    parser.add_argument('--end', default='NONE', help='File number to end on. Default: only download start run')
    parser.add_argument('--extension', default='.fasta', help='File extension (ie .fasta, .fastq, ...) Default: .fasta')
    parser.add_argument('--path', default='./', help='Path to files to be assembled. Default: current directory.')
    parser.add_argument('--location', default='NONE', help='Path to assembler. Required.')
    parser.add_argument('--assembler', default='soap', help='Which assembler to use. Default: "soap" (soapDenovo2).')
    parser.add_argument('--kmer', default=63, help='Kmer value to use for the assembler. Default: 63')
    parser.add_argument('--out', default='./', help='Output directory. Default is current directory.')
    parser.add_argument('-v', '--verbose', default=False, action='store_true', help='Verbose output')
    args = parser.parse_args()
    return args


def main():
    args = parseargs()
    if args.start == 'NONE':
        print "Error: File number to start on must be specified with --start"
        sys.exit()
    if args.location == 'NONE':
        print "Error: Location of assembler executable must be specified with --location"
        sys.exit()
    if args.assembler != 'soap':
        print "Error: Only soapDenovo2 is supported at this time (use default for --assembler)."
        sys.exit()
    try:
        kmer = int(args.kmer)
    except ValueError:
        print "Error: --kmer must be an odd integer between 13 and 127."
        sys.exit()
    if kmer < 13 or kmer > 127 or kmer % 2 != 1:
        print "Error: --kmer must be an odd integer between 13 and 127."
        sys.exit()
    if not args.path.endswith('/'):
        args.path += '/'
    if not args.out.endswith('/'):
        args.out += '/'
    if not args.extension.startswith('.'):
        args.extension = '.' + args.extension
    if args.verbose:
        print "Verbose output requested."

    assemble(args, kmer)


if __name__ == "__main__":
    main()
