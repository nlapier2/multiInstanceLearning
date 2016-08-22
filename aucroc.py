# Author: Nathan LaPierre
# Date: August 18, 2016


import argparse
import sys
import time


def accuracy(tp, fp, tn, fn):
    if tp == 0 and fp == 0 and tn == 0 and fn == 0:
        return 0.0
    else:
        return float(tp + tn) / float(tp + fp + tn + fn)


def f1score(precis, recall):
    if precis == 0.0 and recall == 0.0:
        return 0.0
    else:
        return 2.0 * precis * recall / (precis + recall)


def true_positive_rate(tp, fn):
    if tp == 0 and fn == 0:
        return 0.0
    else:
        return float(tp) / float(tp + fn)


def false_positive_rate(fp, tn):
    if fp == 0 and tn == 0:
        return 0.0
    else:
        return float(fp) / float(fp + tn)


def precision(tp, fp):
    if tp == 0 and fp == 0:
        return 0.0
    else:
        return float(tp) / float(tp + fp)


def parseargs():    # handle user arguments
    parser = argparse.ArgumentParser(description='Generates area under curve based on predicted vs actual labels')
    parser.add_argument('--labels', default='NONE', help='Location of file with predicted labels')
    parser.add_argument('-v', '--verbose', default=False, action='store_true', help='Verbose output')
    args = parser.parse_args()
    return args


def main():
    start = time.time()
    args = parseargs()
    if args.labels == 'NONE':
        print 'File with predicted labels must be specified with --labels. Aborting...'
        sys.exit()

    f = open(args.labels, 'r')
    total, positive, negative = 0, 0, 0     # number of total, positive, and negative labels
    labels = []     # will hold the numerical labels
    for line in f:
        line.strip()
        if not line.startswith('#'):
            l = int(line)
            labels.append(l)
            total += 1
            if l == 1:
                positive += 1
            elif l == -1 or l == 0:
                negative += 1
    f.close()
    tp, fp = 0, 0   # number of true and false positives; they start at 0 since we start with highest possible threshold
    tn = negative   # since everything is negatively predicted at first, true negatives = actual negatives
    fn = positive   # since everything is negatively predicted at first, false negatives = actual positives
    tpr = true_positive_rate(tp, fn)    # this will hold the current tpr / recall at the given decision boundary
    fpr = false_positive_rate(fp, tn)   # this will hold the current false positive rate at the given decision boundary
    prec = precision(tp, fp)            # this will hold the current precision and the given decision boundary
    acc = accuracy(tp, fp, tn, fn)      # this will hold the accuracy achieved from optimally-picked decision boundary
    f1 = f1score(prec, tpr)             # this will hold the f1score achieved from an optimally-picked decision boundary
    auc = 0     # area under curve will initially be 0 but we will increment it in the loop below

    for index in range(total):  # we gradually lower the threshold, converting negative predictions to positive ones
        if labels[index] == 1:
            fn -= 1
            tp += 1
        elif labels[index] == 0 or labels[index] == -1:
            tn -= 1
            fp += 1
            auc += (false_positive_rate(fp, tn) - fpr) * tpr    # calculates area under curve for the change in x value
        tpr = true_positive_rate(tp, fn)
        fpr = false_positive_rate(fp, tn)
        prec = precision(tp, fp)
        if accuracy(tp, fp, tn, fn) > acc:
            acc = accuracy(tp, fp, tn, fn)
        if f1score(prec, tpr) > f1:
            f1 = f1score(prec, tpr)

    print "Accuracy: " + str(acc)
    print "F1 Score: " + str(f1)
    print "AUC-ROC:  " + str(auc)

    if args.verbose:
        print "\nProgram execution time: " + str(time.time() - start) + " seconds"

if __name__ == "__main__":
    main()