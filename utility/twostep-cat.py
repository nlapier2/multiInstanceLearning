# Author: Nathan LaPierre

import argparse
import glob
import os
import time


def parseargs():    # handle user arguments
    parser = argparse.ArgumentParser(description='Script that performs assembly on a number of patient files.')
    parser.add_argument('--out', default='combined.fasta',help='Where to write the output.')
    parser.add_argument('--path', default='', help='Path to files to concatenate. Default: current directory.')
    parser.add_argument('--testset', default='testset.txt', help='File with testset patients.')
    parser.add_argument('-v', '--verbose', default=False, action='store_true', help='Verbose output')
    args = parser.parse_args()
    return args


def main():
    start = time.time()
    args = parseargs()
    testset = []
    tfile = open(args.testset, 'r')
    for line in tfile:
    	testset.append(line.strip())
    tfile.close()
    out = False    # whether to print the line into the outfile
    outfile = open(args.out, 'w')

    # we now write training set lines to the outfile
    for filename in glob.glob(os.path.join(args.path, 'contigs-*.fasta')):
    	prev = 'NONE'  # previous patient, saves computation of checking if patient in testset
    	cfile = open(filename, 'r')  # the contig file
    	for line in cfile:
    		if line.startswith('>'):  # if a defline, check if patient has changed since last defline
    			patient = line.split('>')[1].split('.')[0]
    			if patient != prev:
    				prev = patient  # patient has changed since last line and we write a new prev
    				out = False
    				if patient not in testset:  # write if patient is not in test set
    					out = True 
    		if out:
    			outfile.write(line)
    	if args.verbose:
    		print "Done with: " + filename
    	cfile.close()

    # now we repeat the same thing, but write test set instead of training set
    for filename in glob.glob(os.path.join(args.path, 'contigs-*.fasta')):
    	prev = 'NONE'
    	cfile = open(filename, 'r')
    	for line in cfile:
    		if line.startswith('>'):
    			patient = line.split('>')[1].split('.')[0]
    			if patient != prev:
    				prev = patient
    				out = False
    				if patient in testset:  # the key difference: now write if in testset
    					out = True 
    		if out:
    			outfile.write(line)
    	if args.verbose:
    		print "Done with: " + filename
    	cfile.close()

    outfile.close()
    if args.verbose:
        print "Program execution time: " + str(time.time() - start) + " seconds"


if __name__ == '__main__':
	main()
