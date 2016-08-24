# Author: Nathan LaPierre
# Date: August 23, 2016

import argparse
import sys


def parseargs():    # handle user arguments
    parser = argparse.ArgumentParser(description='Prints ordered predictions for mRMR + SVM for random test sets.')
    parser.add_argument('--labels', default='NONE', help='File with MGWAS labels.')
    parser.add_argument('--patients', default='NONE', help='File with MGWAS patient tags.')
    parser.add_argument('--testset', default='NONE', help='File with patient tags for test set.')
    parser.add_argument('-v', '--verbose', default=False, action='store_true', help='Verbose output')
    args = parser.parse_args()
    return args


def main():
    args = parseargs()
    if args.labels == 'NONE' or args.patients == 'NONE' or args.testset == 'NONE':
        print '--labels, --patients, and --testset must be specified. See -h for help. Aborting...'
        sys.exit()

    labels = []     # holds patient labels
    tags = []       # holds patient tags
    labelfile = open(args.labels, 'r')
    patientfile = open(args.patients, 'r')
    for line in labelfile:
        line = line.strip()
        if not line.startswith('#'):
            patienttag = patientfile.next().strip().split(' = ')[1]
            labels.append(int(line))
            tags.append(patienttag)
    labelfile.close()
    patientfile.close()

    pred_act = []   # holds predicted vs actual labels
    testtags = open(args.testset, 'r')
    for line in testtags:
        line = line.strip()
        if args.verbose:
            print line
        for i in range(len(tags)):
            if line == tags[i]:
                if args.verbose:
                    print "Matched: " + line + ' --> ' + tags[i]
                pred_act.append([len(tags) - i, labels[i]])
    testtags.close()

    pred_act.sort(key=lambda x: x[0])   # sort by predicted label
    for i in range(len(pred_act)):
        print pred_act[len(pred_act)-i-1][1]      # print actual labels in order of predicted labels
    if args.verbose:
        print len(pred_act)

if __name__ == "__main__":
    main()