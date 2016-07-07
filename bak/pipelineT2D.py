# Author: Nathan LaPierre


import argparse
import os
import shlex
import subprocess
import sys
import time
import numpy
import misvm


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
    parser.add_argument('--result_dir', default='./', help='Path to directory to place results in.')
    parser.add_argument('--split', default='.0', help='Characters to split patient tag on.')
    parser.add_argument('--svm', default='', help='Path to svm_light home folder, OR specify "misvm" ')
    parser.add_argument('--train', default='NONE', help='Path to svm training data set')
    parser.add_argument('--test', default='NONE', help='Path to svm test data set')
    parser.add_argument('-v', '--verbose', default=False, action='store_true', help='Verbose output')
    args = parser.parse_args()
    return args


def parse_kraken(vectors, infile, split):  # generate feature vectors from kraken output
    f = open(infile)
    for line in f:
        if line.startswith('C'):
            fields = line.split('\t')       # fields in clustering output lines
            tag = fields[1].split(split[:1])[int(split[1:])]   # gets the tag that matches a patient in the mapping file
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
    f.close()
    return vectors


def parse_uclust(vectors, infile, split):  # generate feature vectors from uclust output
    f = open(infile)
    for line in f:
        if line.startswith('H') or line.startswith('S'):
            fields = line.split('\t')       # fields in clustering output lines
            tag = fields[8].split(split[:1])[int(split[1:])]   # gets the tag that matches a patient in the mapping file
            vec = vectors.get(tag)          # the vector associated with this patient tag
            if vec is not None:
                veclen = len(vec)               # current length of vector's feature list
                featnum = int(fields[1]) + 1    # number of feature is 1 + cluster number
                if veclen < featnum + 1:        # used to expand feature list if needed
                    for j in range(veclen, featnum + 1):
                        vec.append('')
                if line.startswith('H'):
                    value = str(float(fields[3])/100.0)   # "value" of feature, i.e. percentage match to cluster seed
                    if vec[featnum] == '' or float(value) > float(vec[featnum].split(":")[1]):
                        vec[featnum] = str(featnum) + ":" + value + ' '
                elif line.startswith('S'):
                    vec[featnum] = str(featnum) + ":1.0 "
    f.close()
    return vectors


def parse_cluster(verbose, directory, positive, negative, cluster, infile, mapfile, split):
    # generate feature vectors from cluster output
    os.chdir(directory)
    if verbose:
        print "Reading map: " + mapfile

    vectors = {}     # dict that holds patient vectors to be used for svm
    m = open(mapfile)
    for line in m:
        fields = line.split('\t')
        if fields[7] == positive:   # if this patient is a positive example of the state
            vectors[fields[3]] = ['1 ']
        elif negative == "NONE" or fields[7] == negative:
            vectors[fields[3]] = ['-1 ']

    if verbose:
        print "Reading input: " + infile
    if cluster == 'kraken':
        vectors = parse_kraken(vectors, infile, split)
    elif cluster == 'uclust':
        vectors = parse_uclust(vectors, infile, split)
    return vectors


def write_train_test(vectors, train, test):     # write training and test files for SVM
    tr = open(train, 'w')
    te = open(test, 'w')
    counter = 0     # every Nth patient will be put into the test set instead of the train set
    for key, value in vectors.iteritems():
        counter = (counter + 1) % 2
        if counter != 0:
            if value[0].startswith('1'):
                for j in range(1, len(value)):      # feature engineering loop
                    if value[j] != '':
                        parts = value[j].split(":")
                        value[j] = parts[0] + ":" + str(float(parts[1]) * 0.9) + ' '    # apply multiplier to features

            tr.write((''.join(value)) + '\n')   # write to training set
        else:
            te.write((''.join(value)) + '\n')   # write to test set


def execute_svm(verbose, directory, train, test, model, prediction, svm, output):
    # execute call in this format:
    # svm_learn example1/train.dat example1/model
    # svm_classify example1/test.dat example1/model example1/predictions
    if verbose:
        print "SVM Classifying..."
    if svm != '' and not svm.endswith('/'):     # make directory ending proper
        svm += '/'
    os.chdir(directory)
    if output == 'NONE':    # if using stdout
        subprocess.call(shlex.split('./'+svm+'svm_learn -x 1 '+train+' '+model))     # create svm model
        subprocess.call(shlex.split('./'+svm+'svm_classify '+test+' '+model+' '+prediction))    # svm predictions
    else:   # if output file specified
        with open(output, 'w') as outfile:
            subprocess.call(shlex.split('./'+svm+'svm_learn -x 1 '+train+' '+model), stdout=outfile)
            subprocess.call(shlex.split('./'+svm+'svm_classify '+test+' '+model+' '+prediction), stdout=outfile)


def misvm_classify(verbose, output, vectors, labels):
    # perform the actual misvm classification
    if verbose:
        print "Creating train and test bags and labels..."
    bags = [numpy.array(vectors[v], dtype=float) for v in vectors]      # numpy-format matrix for use in misvm
    labels = numpy.array([labels[l] for l in labels], dtype=float)      # numpy-format labels for use in misvm
    # Spilt dataset into train and test sets
    train_bags = []  # bags[midpoint:]  # []  # bags[10:]
    train_labels = []  # labels[midpoint:]  # []  # labels[10:]
    test_bags = []  # bags[:midpoint]  # []  # bags[:10]
    test_labels = []  # labels[:midpoint]  # []  # labels[:10]
    for i in range(len(labels)):
        if i % 2 == 0:
            train_bags.append(bags[i])
            train_labels.append(labels[i])
        else:
            test_bags.append(bags[i])
            test_labels.append(labels[i])

    if verbose:
        print "MISVM Classifying..."
    if output != 'NONE':
        sys.stdout = open(output, 'w')
    # establish classifiers
    classifiers = {
        'sbMIL': misvm.sbMIL(kernel='linear', eta=0.1, C=1.0),
        'SIL': misvm.SIL(kernel='linear', C=1.0),
        'MISVM': misvm.MISVM(kernel='linear', C=1.0, max_iters=100),
    }
    # Train/Evaluate classifiers
    accuracies = {}
    for algorithm, classifier in classifiers.items():
        classifier.fit(train_bags, train_labels)
        predictions = classifier.predict(test_bags)
        accuracies[algorithm] = numpy.average(test_labels == numpy.sign(predictions))
    for algorithm, accuracy in accuracies.items():
        print '\n%s Accuracy: %.1f%%' % (algorithm, 100 * accuracy)
    if output != 'NONE':
        sys.stdout = sys.__stdout__     # reset stdout to normal


def execute_misvm(verbose, directory, output, positive, negative, cluster, infile, mapfile, split):
    # execute misvm, assuming it has been properly installed
    if cluster == 'kraken':
        print 'Kraken MISVM option not currently supported.'
        sys.exit()
    os.chdir(directory)
    if verbose:
        print "Reading map: " + mapfile
    vectors = {}        # dict that holds patient matrices to be used for misvm
    labels = {}         # dict that holds labels for patients ("bags")
    m = open(mapfile)
    for line in m:
        fields = line.split('\t')
        if fields[7] == positive:       # if this patient is a positive example of the state
            vectors[fields[3]] = []     # initialize matrix for patient
            labels[fields[3]] = 1.0     # assign label for patient
        elif negative == "NONE" or fields[7] == negative:
            vectors[fields[3]] = []     # initialize matrix for patient
            labels[fields[3]] = -1.0    # assign label for patient

    if verbose:
        print "Reading input: " + infile
    f = open(infile)
    for line in f:
        if line.startswith('H') or line.startswith('S'):
            fields = line.split('\t')       # fields in clustering output lines
            tag = fields[8].split(split[:1])[int(split[1:])]   # gets the tag that matches a patient in the mapping file
            vec = vectors.get(tag)          # the vector associated with this patient tag
            if vec is not None:
                vec.append([float(fields[1])])  # for each read, append a 1-element list with cluster number
                if len(vec) < 1:
                    vec.append([float(fields[1])])  # for each read, append a 1-element list with cluster number
    f.close()
    for v in vectors:
        if not vectors[v]:
            vectors[v] = [[0.0]]
    misvm_classify(verbose, output, vectors, labels)    # perform the actual misvm classification


def main():
    start = time.time()
    args = parseargs()
    if (args.input == "NONE" or args.map == "NONE") and (args.train == "NONE" or args.test == "NONE"):
            print "You must specify input and map files or training and test files. Use -h for help. Exiting..."
            sys.exit()
    # make sure all output files go in the specified results directory
    if not args.result_dir.endswith('/'):
        args.result_dir += '/'
    if not args.output.startswith('/'):
        args.output = args.result_dir + args.output
    if not args.model.startswith('/'):
        args.model = args.result_dir + args.model
    if not args.predictions.startswith('/'):
        args.predictions = args.result_dir + args.predictions

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
        if not args.train.startswith('/'):
            args.train = args.result_dir + args.train
        if args.test == "NONE":
            args.test = "test.dat"
        if not args.test.startswith('/'):
            args.test = args.result_dir + args.test
        if args.verbose:
            print "Entering input mode..."
        if args.svm != 'misvm':
            vectors = parse_cluster(args.verbose, args.dir, args.positive, args.negative,
                                    args.cluster, args.input, args.map, args.split)
            write_train_test(vectors, args.train, args. test)

    if args.svm == 'misvm':
        execute_misvm(args.verbose, args.dir, args.output,
                      args.positive, args.negative, args.cluster, args.input, args.map, args.split)
    else:
        # run SVM on given or generated train and test files
        execute_svm(args.verbose, args.dir, args.train, args.test, args.model, args.predictions, args.svm, args.output)
    if args.verbose:
        print "Program execution time: " + str(time.time() - start) + " seconds"


if __name__ == "__main__":
    main()
