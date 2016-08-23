import argparse
import time
import sys


def parseargs():    # handle user arguments
    parser = argparse.ArgumentParser(description='Clean unwanted whitespace out of FASTA files.')
    parser.add_argument('-i', '--infile', default='NONE', help='Location of input file')
    parser.add_argument('-o', '--output', default='NONE', help='Where to write output to')
    parser.add_argument('-v', '--verbose', default=False, action='store_true', help='Verbose output')
    args = parser.parse_args()
    return args


def main():
    start = time.time()
    args = parseargs()
    if args.infile == 'NONE' or args.output == 'NONE':
        print "Location of input and output files must be specified."
        sys.exit()
    infile = open(args.infile, 'r')
    output = open(args.output, 'w')
    read = ''
    for line in infile:
	if line.startswith('>'):
	    if read != '':
		output.write(read + '\n')
		read = ''
	    output.write(line)
	else:
	    read += ''.join(line.split())
    if read != '':
	output.write(read + '\n')
    infile.close()
    output.close()
    if args.verbose:
	print "Program execution time: " + str(time.time() - start) + " seconds"


if __name__ == "__main__":
    main()

