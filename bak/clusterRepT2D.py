# Author: Nathan LaPierre
# Date: March 17, 2016

import argparse
import os
import shlex
import subprocess
import sys
import time


def option1(indir, output, alpha, alphaopt, uclustid, uclustpath, verbose):
    print indir + output + alpha + alphaopt + uclustid + uclustpath + verbose
    print "Option not yet implemented."


def option2(indir, output, alpha, alphaopt, uclustid, uclustpath, verbose):
    print indir + output + alpha + alphaopt + uclustid + uclustpath + verbose
    print "Option not yet implemented."


def option3(indir, output, alpha, alphaopt, uclustid, uclustpath, verbose):
    tempfilename = 'ccclusterreptempfile'
    tempfile = open(tempfilename, 'w')      # temporary helper file, deleted at the end
    try:
        clustlist = []      # list of patients and cluster numbers, i.e. patient A cluster 1 becomes "A/1"
        for filename in os.listdir(indir):
            if filename.endswith('.uc'):
                if verbose:
                    print 'Attempting to open ' + filename + ' and ' + filename.split('.')[0] + '.fasta'
                fasta = open(indir + filename.split('.')[0] + '.fasta')  # open corresponding fasta file for uclust file
                f = open(indir + filename)
                for line in f:
                    if line.startswith('H') or line.startswith('S'):
                        fastaline = fasta.next()
                        if line.startswith('S'):        # cluster seed: write to temp file
                            fields = line.split('\t')
                            clustlist.append(fields[8].split('.')[0] + '/' + fields[1])  # for patient A clust 1, "A/1"
                            tempfile.write(fastaline)   # write the defline
                            fastaline = fasta.next()    # go to read line
                            tempfile.write(fastaline)   # write the read
                        else:
                            fasta.next()    # read nothing
                f.close()
                fasta.close()
        tempfile.close()
        if verbose:
            print 'Calling: ' + uclustpath + ' --input ' + os.curdir + '/' + tempfilename \
                  + ' --uc ' + output + ' --usersort --id ' + str(uclustid)
        subprocess.call(shlex.split(uclustpath + ' --input ' + os.curdir + '/' + tempfilename
                                    + ' --uc ' + output + ' --usersort --id ' + str(uclustid)))  # call uclust
        outfile = open(output)
        index = 0   # used to index the cluster list created earlier
        for line in outfile:    # we now go through uclust output to see which original clusters have been grouped
            if line.startswith('H') or line.startswith('S'):
                clustlist[index] += ('->' + line.split('\t')[1])    # if "A/1" is now cluster 3, it becomes "A/1->3"
                index += 1  # move to next spot in list
        outfile.close()     # close file after iterating through it
        outfile = open(output, 'w')  # now reopen to write our final output to it
        clustdict = {}   # dict version of the list we created so we can constant time access
        for clust in clustlist:
            clustdict[clust.split('->')[0]] = clust.split('->')[1]  # for "A/1->3", "A/1" becomes key with "3" value
        for filename in os.listdir(indir):
            if filename.endswith('.uc'):
                if verbose:
                    print 'Attempting to open ' + filename
                f = open(indir + filename)
                if alphaopt == 1:
                    count = 0   # used for selecting 1 % alpha lines
                    for line in f:
                        if line.startswith('H') or line.startswith('S'):
                            count = 1 + (count % alpha)
                            if count == 1:  # this occurs every (1 / alpha) lines, and on the first line of each cluster
                                fields = line.split('\t')
                                fields[1] = clustdict[fields[8].split('.')[0] + '/' + fields[1]]  # reassign clust num
                                outfile.write('\t'.join(fields))  # write line in original form with new cluster number
                elif alphaopt == 2:
                    counts = [0]
                    for line in f:
                        if line.startswith('H') or line.startswith('S'):
                            fields = line.split('\t')
                            clustnum = int(fields[1])
                            countlen = len(counts)
                            if clustnum > countlen:
                                for j in range(countlen, clustnum):
                                    counts.append(0)
                            counts[clustnum-1] = 1 + (counts[clustnum-1] % alpha)
                            if counts[clustnum-1] == 1:
                                fields[1] = clustdict[fields[8].split('.')[0] + '/' + fields[1]]  # reassign clust num
                                outfile.write('\t'.join(fields))  # write line in original form with new cluster number
                f.close()
        os.remove(tempfilename)   # remove our temporary file
    except:
        tempfile.close()
        os.remove(tempfilename)   # remove our temporary file
        print "Unexpected error: ", sys.exc_info()[0]
        raise


def option4(indir, output, alpha, alphaopt, uclustid, uclustpath, verbose):
    tempfilename = 'ccclusterreptempfile'
    tempfile = open(tempfilename, 'w')      # temporary helper file, deleted at the end
    try:
        for filename in os.listdir(indir):
            if filename.endswith('.uc'):
                if verbose:
                    print 'Attempting to open ' + filename + ' and ' + filename.split('.')[0] + '.fasta'
                fasta = open(indir + filename.split('.')[0] + '.fasta')  # open corresponding fasta file for uclust file
                f = open(indir + filename)
                if alphaopt == 1:
                    count = 0   # used for selecting 1 % alpha lines
                    for line in f:
                        if line.startswith('H') or line.startswith('S'):
                            fastaline = fasta.next()
                            count = 1 + (count % alpha)
                            if count == 1:  # this occurs every (1 / alpha) lines
                                tempfile.write(fastaline)   # write the defline
                                fastaline = fasta.next()    # go to read line
                                tempfile.write(fastaline)   # write the read
                            else:
                                fasta.next()    # skip the read line; write nothing
                elif alphaopt == 2:
                    counts = [0]
                    for line in f:
                        if line.startswith('H') or line.startswith('S'):
                            fastaline = fasta.next()
                            countlen = len(counts)
                            clustnum = int(line.split('\t')[1])
                            if clustnum > countlen:
                                for j in range(countlen, clustnum):
                                    counts.append(0)
                            counts[clustnum-1] = 1 + (counts[clustnum-1] % alpha)
                            if counts[clustnum-1] == 1:  # this occurs every (1 / alpha) lines per cluster
                                tempfile.write(fastaline)   # write the defline
                                fastaline = fasta.next()    # go to read line
                                tempfile.write(fastaline)   # write the read
                            else:
                                fasta.next()    # skip the read line; write nothing
                f.close()
                fasta.close()
        if verbose:
            print 'Calling: ' + uclustpath + ' --input ' + os.curdir + '/' + tempfilename \
                  + ' --uc ' + output + ' --usersort --id ' + str(uclustid)
        subprocess.call(shlex.split(uclustpath + ' --input ' + os.curdir + '/' + tempfilename
                                    + ' --uc ' + output + ' --usersort --id ' + str(uclustid)))  # call uclust
        tempfile.close()
        os.remove(tempfilename)   # remove our temporary file
    except:
        tempfile.close()
        os.remove(tempfilename)   # remove our temporary file
        print "Unexpected error: ", sys.exc_info()[0]
        raise


def create_representation(indir, output, alpha, alphaopt, option, uclustid, uclustpath, verbose):
    # Given the input cluster file(s), generate output file given alpha value
    if verbose:
        print 'Opening output file specified: ' + output
        print 'Input directory specified: ' + indir
    if not indir.endswith('/'):
        indir += '/'

    if option == 1:
        option1(indir, output, alpha, alphaopt, uclustid, uclustpath, verbose)
    elif option == 2:
        option2(indir, output, alpha, alphaopt, uclustid, uclustpath, verbose)
    elif option == 3:
        option3(indir, output, alpha, alphaopt, uclustid, uclustpath, verbose)
    elif option == 4:
        option4(indir, output, alpha, alphaopt, uclustid, uclustpath, verbose)


def parseargs():    # handle user arguments
    parser = argparse.ArgumentParser(description='Clustering with small representation for large clusters.')
    parser.add_argument('--indir', default='NONE', help='Input Directory with UCLUST and FASTA files. Required.')
    parser.add_argument('--output', default='out.uclust', help='Output File. Default is "out.uc" in current directory.')
    parser.add_argument('--alpha', default='100', help='(1 / alpha) reads are taken from each cluster. Default: 100.')
    parser.add_argument('--alphaopt', default=1, help='Options for selecting (1 / alpha) reads. 1 (default), 2, or 3.')
    parser.add_argument('--option', default=4, help='Clustering option (1,2,3, or 4); default 4')
    parser.add_argument('--uclustpath', default='NONE', help='path to uclust. required for option 4.')
    parser.add_argument('--uclustid', default=0.9, help='id match to put reads in same cluster. 0.0-1.0, default 0.9')
    parser.add_argument('-v', '--verbose', default=False, action='store_true', help='Verbose output')
    args = parser.parse_args()
    return args


def main():
    start = time.time()
    args = parseargs()
    if args.indir == 'NONE':
        print "Error: Input directory (--indir) must be specified."
        sys.exit()
    if args.option == 4 and args.uclustpath == 'NONE':
        print "Error: For option 4 (default), you must specify the path to the uclust program with --uclustpath."
        sys.exit()
    try:
        option = int(args.option)
    except ValueError:
        print "Error: Option must be 1, 2, 3, or 4."
        sys.exit()
    if option < 1 or option > 4:
        print "Error: Option must be 1, 2, 3, or 4."
        sys.exit()
    try:
        alpha = int(args.alpha)
    except ValueError:
        print "Error: Alpha must be a natural number (positive integer)."
        sys.exit()
    if alpha < 1:
        print "Error: Alpha must be a natural number (positive integer)."
        sys.exit()
    try:
        alphaopt = int(args.alphaopt)
    except ValueError:
        print "Error: Alphaopt must be 1, 2, or 3."
        sys.exit()
    if alphaopt < 1 or alphaopt > 3:
        print "Error: Alphaopt must be 1, 2, or 3."
        sys.exit()
    try:
        uclustid = float(args.uclustid)
    except ValueError:
        print "Error: uclustid must be a float between 0.0 and 1.0, inclusive"
        sys.exit()
    if uclustid < 0.0 or uclustid > 1.0:
        print "Error: uclustid must be a float between 0.0 and 1.0, inclusive"
        sys.exit()
    if args.verbose:
        print "Verbose output requested."
    create_representation(args.indir, args.output, alpha, alphaopt, option,  uclustid, args.uclustpath, args.verbose)
    if args.verbose:
        print "Program execution time: " + str(time.time() - start) + " seconds"


if __name__ == "__main__":
    main()
