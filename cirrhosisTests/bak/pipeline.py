# Author: Nathan LaPierre
# Version: 1.0.0
# Date: 1 May 2015


import argparse
import os
import shlex
import subprocess
import sys
import time


def parseargs():    # handle user arguments
    parser = argparse.ArgumentParser(description='Classify fasta-format sequence reads')
    parser.add_argument('-c', '--cluster', default='NONE', choices=['kraken', 'uclust'],
                        help='Which clustering algorithm was used')
    parser.add_argument('-d', '--dir', default='NONE', help='You can specify directory of all parameters if shared')
    parser.add_argument('-i', '--input', default='NONE', help='Path to clustering output file to be classified')
    parser.add_argument('--map', default='NONE', help='Path to mapping file to be used')
    parser.add_argument('--model', default='model', help='Path to place SVM model file')
    parser.add_argument('--negative', default='NONE', help='negative state to be classified (i.e. Encephalopathy, etc.')
    parser.add_argument('-o', '--output', default='NONE', help='Path to output file (stdout if not used)')
    parser.add_argument('--positive', default='NONE', help='positive state to be classified (i.e. Encephalopathy, etc.')
    parser.add_argument('--predictions', default='predictions', help='Path to place SVM predictions file')
    parser.add_argument('--svm', default='', help='Path to svm_light home folder')
    parser.add_argument('--train', default='NONE', help='Path to svm training data set')
    parser.add_argument('--test', default='NONE', help='Path to svm test data set')
    parser.add_argument('-v', '--verbose', default=False, action='store_true', help='Verbose output')
    args = parser.parse_args()
    return args


def parseCluster(verbose, dir, positive, negative, cluster, input, map, train, test):
    os.chdir(dir)

    if verbose:
        print "Reading map: " + map

    vectors = {}     # dict that holds patient vectors to be used for svm

    m = open(map)
    for line in m:
        fields = line.split('\t')
        if fields[7] == positive:   # if this patient is a positive example of the state
            # poscount += 1
            # if poscount < negcount + 2:
            vectors[fields[3]] = ['1 ']
        elif negative == "NONE" or fields[7] == negative:
            # negcount += 1
            # if negcount < poscount + 2:
            vectors[fields[3]] = ['-1 ']

    if verbose:
        print "Reading input: " + input

    f = open(input)

    if cluster == 'kraken':
        for line in f:
            if line.startswith('C'):
                fields = line.split('\t')       # fields in clustering output lines
                tag = fields[1].split('-')[1]   # gets the tag that matches a patient in the mapping file
                vec = vectors.get(tag)          # the vector associated with this patient tag
                if vec is not None:
                    veclen = len(vec)               # current length of vector's feature list
                    seqlen = int(fields[3])-30      # length of the sequence read represented by this line
                    otu = fields[2]                 # the OTU kraken finds this read to represent
                    otunum = int(otu)               # numerical version of otu variable
                    matchnum = 0                    # number of k-mers that match to OTU
                    kmers = fields[4].split(' ')    # the sequence of k-mers in kraken found for this read
                    for k in kmers:
                        subk = k.split(":")         # k-mer split into OTU and length
                        if subk[0] == otu:
                            matchnum += int(subk[1])
                    if veclen < otunum + 1:         # used to expand feature list if needed
                        for j in range(veclen, otunum + 1):
                            vec.append('')
                    value = float(matchnum) / float(seqlen)  # "value" of feature i.e. percentage of k-mers matching OTU
                    if vec[otunum] == '' or value > float(vec[otunum].split(":")[1]):
                        vec[otunum] = otu + ":" + str(value) + ' '
                    # else:   # add to current value if it already exists
                    #    vec[otunum] = otu + ":" + str(value + float(vec[otunum].split(":")[1])) + ' '

    elif cluster == 'uclust':
        for line in f:
            if line.startswith('H'):            # not cluster seed
                fields = line.split('\t')       # fields in clustering output lines
                tag = fields[8].split('-')[1]   # gets the tag that matches a patient in the mapping file
                vec = vectors.get(tag)          # the vector associated with this patient tag
                if vec is not None:
                    veclen = len(vec)               # current length of vector's feature list
                    featnum = int(fields[1]) + 1    # number of feature is 1 + cluster number
                    value = str(float(fields[3])/100.0)   # "value" of feature, i.e. percentage match to cluster seed
                    if veclen < featnum + 1:        # used to expand feature list if needed
                        for j in range(veclen, featnum + 1):
                            vec.append('')
                    if vec[featnum] == '' or float(value) > float(vec[featnum].split(":")[1]):
                        vec[featnum] = str(featnum) + ":" + value + ' '
                    # else:   # add to current value if it already exists
                    #    vec[featnum] = str(featnum) + ":" + str(float(value) + float(vec[featnum].split(":")[1])) + ' '

            elif line.startswith('S'):          # cluster seed
                fields = line.split('\t')       # fields in clustering output lines
                featnum = int(fields[1]) + 1    # number of feature is 1 + cluster number
                tag = fields[8].split('-')[1]   # gets the tag that matches a patient in the mapping file
                vec = vectors.get(tag)          # the vector associated with this patient tag
                if vec is not None:
                    veclen = len(vec)               # current length of vector's feature list
                    if veclen < featnum + 1:
                        for k in range(veclen, featnum + 1):
                            vec.append('')
                    vec[featnum] = str(featnum) + ":1.0 "

    tr = open(train, 'w')
    te = open(test, 'w')
    poscount = 0            # number of positive examples in training set
    negcount = 0            # number of negative examples in training set
    traincount = 0          # number of training examples
    testcount = 0           # number of test examples
    counter = 0     # every fifth patient will be put into the test set instead of the train set
    for key, value in vectors.iteritems():
        counter = (counter + 1) % 5
        if counter != 0:

            if value[0].startswith('1'):
                for j in range(1, len(value)):      # saturate features of positive examples in training set
                    if value[j] != '':
                        parts = value[j].split(":")
                        value[j] = parts[0] + ":" + str(float(parts[1]) * 0.9) + ' '    # apply multiplier to features

            tr.write((''.join(value)) + '\n')   # write to training set
            '''
            if value[0].startswith('1'):
                poscount += 1
                if poscount < negcount + 2 and negcount < poscount + 2:
                    traincount += 1
                    tr.write((''.join(value)) + '\n')   # write to training set
                else:
                    poscount -= 1
            elif value[0].startswith('-1'):
                negcount += 1
                if poscount < negcount + 2 and negcount < poscount + 2:
                    traincount += 1
                    tr.write((''.join(value)) + '\n')   # write to training set
                else:
                    negcount -= 1
            '''
        else:
            # if testcount <= traincount / 5:
            #    testcount += 1
            te.write((''.join(value)) + '\n')   # write to test set


def executeSVM(verbose, dir, train, test, model, prediction, svm, output):
    # execute call in this format:
    # svm_learn example1/train.dat example1/model
    # svm_classify example1/test.dat example1/model example1/predictions

    if verbose:
        print "Classifying..."

    if svm != '' and not svm.endswith('/'):     # make directory ending proper
        svm += '/'

    os.chdir(dir)
    # print dir
    # print './'+svm+'svm_learn '+train+' '+model+' '
    if output == 'NONE':    # if using stdout
        subprocess.call(shlex.split('./'+svm+'svm_learn -x 1 '+train+' '+model))     # create svm model
        subprocess.call(shlex.split('./'+svm+'svm_classify '+test+' '+model+' '+prediction))    # svm predictions
    else:   # if output file specified
        with open(output, 'w') as outfile:
            subprocess.call(shlex.split('./'+svm+'svm_learn -x 1 '+train+' '+model), stdout=outfile)
            subprocess.call(shlex.split('./'+svm+'svm_classify '+test+' '+model+' '+prediction), stdout=outfile)


def main():
    start = time.time()
    args = parseargs()
    if (args.input == "NONE" or args.map == "NONE") and (args.train == "NONE" or args.test == "NONE"):
            print "You must specify input and map files or training and test files. Use -h for help. Exiting..."
            sys.exit()

    if args.input != "NONE":    # go into input mode, taking clustering output file and breaking into train and test
        if args.cluster == "NONE":
            print "You must specify which clustering algorithm was used. Choices: kraken, uclust. " \
                  "Use -h for help. Exiting..."
            sys.exit()
        if args.positive == "NONE":
            print "You must specify what positive state we are trying to classify. Use -h for help. Exiting..."
            sys.exit()
        if args.train == "NONE":
            args.train = "train.dat"
        if args.test == "NONE":
            args.test = "test.dat"
        if args.verbose:
            print "Entering input mode..."
        parseCluster(args.verbose, args.dir, args.positive, args.negative, args.cluster, args.input, args.map,
                     args.train, args.test)

    # run SVM on given or generated train and test files
    executeSVM(args.verbose, args.dir, args.train, args.test, args.model, args.predictions, args.svm, args.output)
    if args.verbose:
        print "Program execution time: " + str(time.time() - start) + " seconds"


if __name__ == "__main__":
    main()
