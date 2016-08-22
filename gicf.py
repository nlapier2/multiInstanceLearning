# Author: Nathan LaPierre
# Date: June 6, 2016


import argparse
import math
import numpy as np
import random
import sys
import time


def similarity(i, j):   # similarity measure (0 to 1) between instances i and j
    norm = 0.0    # the squared 2-norm of the difference between the two instance vectors, computed in for loop below
    for k in range(len(i)):
        norm += (i[k] - j[k]) ** 2
    return math.exp(-norm)  # / len(i) / len(j)      # returns value for kernel described in GICF paper


def instance_penalty(i, j):     # non-negative penalty on the difference between predictions for instances i and j
    return (i - j) ** 2     # using squared loss as penalty, as described in GICF paper


def logistic(w, x):     # implementation of the classifier for the instances; currently logistic regression
    dot = w[0]     # dot product between weights and x, computed in for loop below
    for i in range(len(x)):
        dot += w[i+1]*x[i]
    if dot > 709:    # prevent math overflows
        return 1
    if dot < -709:   # prevent math overflows
        return 0
    return 1 / (1 + math.exp(-dot))


def label_penalty(predicted, actual):   # non-negative penalty on difference between prediction & true label for group
    return (predicted - actual) ** 2     # using squared loss as penalty, as described in GICF paper


def predict_label(group, instance_labels):
    # real scalar representing predicted label for group based on aggregation of instance labels
    total = 0.0     # combined value of all labels in group, computed in for loop below
    for i in range(len(group)):
        total += instance_labels[i]
    return total / len(group)   # return average label value in the group


def d_group_penalty(group, instance_labels, d_instance_labels, group_label):
    # calculate group penalty in derivative of cost function
    total1 = 0.0    # the part of the group penalty using the actual group labels
    total2 = np.array([0.0 for i in range(len(d_instance_labels[0]))])   # part using only predicted instance labels
    for i in range(len(group)):
        total1 += instance_labels[i]
        total2 += d_instance_labels[i]
    total1 = 2.0 * ((total1 / len(group)) - group_label)
    total2 /= len(group)
    return total1 * total2


def cost(num_instances, w, data, l, group_labels):
    # calculate cost function as defined in GICF paper
    instance_cost = 0.0     # first part of cost function, based on instance differences
    group_cost = 0.0    # second part of cost function, based on group labels

    instance_labels = []
    for i in range(len(data)):
        instance_labels.append([])
        for j in range(len(data[i])):
            instance_labels[i].append(logistic(w, data[i][j]))

    for group in range(len(data)):
        group_cost += label_penalty(predict_label(data[group], instance_labels[group]), group_labels[group])
        for i in range(len(data[group])):
            for j in range(len(data[group])):
                instance_cost += similarity(data[group][i], data[group][j]) * \
                    instance_penalty(instance_labels[group][i], instance_labels[group][j])

    instance_cost /= (float(num_instances) ** 2)   # average out instance cost across all instance pairs
    group_cost *= (l / float(len(data)))    # balance group cost using lambda parameter

    return instance_cost + group_cost


def derivative_cost(num_instances, w, data, l, group_labels):      # derivative of cost function above
    instance_labels = []    # labels of instances according to weights and logistic function
    d_instance_labels = []  # labels of instances according to weights and derivative of logistic function
    for i in range(len(data)):
        instance_labels.append([])
        d_instance_labels.append([])
        for j in range(len(data[i])):
            instance_labels[i].append(logistic(w, data[i][j]))
            d_instance_labels[i].append(logistic(w, data[i][j]) * (1 - logistic(w, data[i][j])) * np.array(data[i][j]))

    instance_cost = np.array([0.0 for i in range(1, len(w))])  # 1st part of cost function: instance differences
    group_cost = 0      # second part of cost function, based on group labels
    for group in range(len(data)):
        group_cost += d_group_penalty(data[group], instance_labels[group],
                                      d_instance_labels[group], group_labels[group])
        for i in range(len(data[group])):
            for j in range(len(data[group])):
                instance_cost += similarity(data[group][i], data[group][j]) * \
                    (instance_labels[group][i] - instance_labels[group][j]) * \
                    (d_instance_labels[group][i] - d_instance_labels[group][j])
                '''instance_cost += (instance_labels[group][i] - instance_labels[group][j]) * \
                    (d_instance_labels[group][i] - d_instance_labels[group][j])'''
    instance_cost /= (float(num_instances) ** 2)   # average out instance cost across all instance pairs
    group_cost *= (l / float(len(data)))    # balance group cost using lambda parameter
    return instance_cost + group_cost


def descent(verbose, rate, l, top_k, batch_size, iterations, w, data):
    # mini-batch stochastic gradient descent to learn weights for classifier
    # randomly select a small subset of instances from data to be used as mini batch
    if verbose:
        print "Entering gradient descent..."
    update_vector = []  # this is used to track the previous update vector, used for momentum
    momentum = 0.5      # constance multiplied by the previous update vector for momentum
    for i in range(len(w)):
        update_vector.append(0)  # initialize to same length as weight vector, with all 0s
    for num in range(iterations):
        if verbose:
            print "Gradient descent, instance " + str(num + 1)
        cost_data = []  # the data to be given to the cost function
        labels = []     # group labels
        instances = 0.0       # total number of instances used
        if top_k >= 0.99:
            chance = 0.0    # chance of read being selected for the mini batch
        else:  # randomly select batch_size% of reads for mini batch, from top_k
            chance = 1 - ((1 / (100 * (1 - top_k))) * (batch_size / 0.01))
        d = open(data, 'r')
        for line in d:
            if len(line) < 5:   # deals with blank lines in files
                continue
            cost_data.append([])
            splits = line.split(':')
            temp_label = int(float(splits[0]))
            if temp_label == -1:
                labels.append(0)
            else:
                labels.append(1)
            group = listify(splits[1])
            cutoff = int(top_k * len(group))     # use only top (100 * (1 - top_k))% of reads
            for i in range(cutoff, len(group)):
                if random.random() > chance:
                    cost_data[len(labels)-1].append(group[i])   # append appropriate instance to mini batch dataset
                    instances += 1.0
        d.close()

        # find value of derivative of cost function using weights w, then use learning rate to find adjustment amount
        delta_j = rate * derivative_cost(instances, w, cost_data, l, labels)
        # once slope becomes non-negative (or barely negative), stop; or if we exceed number of iterations
        if verbose:
            print "\nDelta J: " + str(delta_j / rate)
            print "\nRate: " + str(rate)

        # adjust weights according to slope of cost function
        for i in range(1, len(w)):
            adjustment = delta_j[i-1] + (update_vector[i] * momentum)
            w[i] -= adjustment
            update_vector[i] = adjustment
        if verbose:
            print "\nCost: " + str(cost(instances, w, cost_data, l, labels)) + '\n'
    if verbose:
        print "Gradient descent final weights: "
        print w
    return w


def validate(verbose, w, data, k):     # apply the weights determined by gradient descent to a test set
    if verbose:
        print "Validating against test set..."
    predicted_actual = []   # array containing predicted and actual labels for each group
    d = open(data, 'r')
    for line in d:
        if len(line) < 5:   # deals with blank lines in files
            continue
        splits = line.split(':')
        temp_label = int(float(splits[0]))
        if temp_label == -1:
            label = 0
        else:
            label = 1
        group = listify(splits[1])
        total = 0.0     # total of instance label predictions
        # k = 0.99
        top_k = int(k * len(group))  # * 0    # TEMP: * 0
        for i in range(top_k, len(group)):
            total += logistic(w, group[i])
        total /= (len(group) - top_k)
        predicted_actual.append([total, label])     # predicted vs actual label
    d.close()

    # threshold = 0
    average_predicted = 0.0
    for i in predicted_actual:
        average_predicted += i[0]
    average_predicted /= len(predicted_actual)
    threshold = average_predicted
    if verbose:
        print "Average prediction: " + str(average_predicted)

    true_positives = 0
    false_positives = 0
    true_negatives = 0
    false_negatives = 0
    for i in predicted_actual:
        if i[0] >= threshold and i[1] == 1:
            true_positives += 1
        elif i[0] > threshold and i[1] == 0:
            false_positives += 1
        elif i[0] < threshold and i[1] == 0:
            true_negatives += 1
        else:
            false_negatives += 1
    accuracy = float(true_positives + true_negatives) / float(len(predicted_actual))
    if true_positives + false_positives == 0:
        precision = 0
    else:
        precision = float(true_positives) / float(true_positives + false_positives)
    tpr_recall = float(true_positives) / float(true_positives + false_negatives)
    if false_positives == 0 and true_negatives == 0:
        fpr = 0
    else:
        fpr = float(false_positives) / float(false_positives + true_negatives)
    if precision == 0 and tpr_recall == 0:
        f1_score = 0
    else:
        f1_score = 2.0 * precision * tpr_recall / (precision + tpr_recall)
    if verbose:
        print "Accuracy: " + str(accuracy)
        print "Precision and recall: " + str(precision) + " and " + str(tpr_recall)
        print "F1 Score: " + str(f1_score)
        print "True Positive Rate and False Positive Rate: " + str(tpr_recall) + " and " + str(fpr)

    return accuracy, precision, tpr_recall, fpr, f1_score, predicted_actual


def grid_search(verbose, w, data, train, test, res, limit, settings, testset):
    # linear grid search to set learning rate, lambda, and mini batch size parameters
    if verbose:
        print "Entering grid search to select best parameters for gradient descent."
    test_split = 2   # 1 out of every test_split instances are put in test set (i.e. 1/2 if test_split = 2)

    f = open(data, 'r')
    tr = open(train, 'w')
    te = open(test, 'w')
    num = 0     # count of the number of groups
    testitems = []
    if testset != 'NONE':
        setfile = open(testset, 'r')
        for line in setfile:
            testitems.append(int(line.strip()))
        setfile.close()

    for line in f:
        num += 1
        if testset != 'NONE':
            if num in testitems:
                te.write(line + '\n')
            else:
                tr.write(line + '\n')
        else:
            if num % test_split == 0:
                te.write(line + '\n')
            else:
                tr.write(line + '\n')
        if num == limit:
            break
    f.close()
    tr.close()
    te.close()

    if settings == 'NONE':
        rate = [0.0001]  # , 0.001]  # , 0.01, 0.1]
        l = [float(num)]  # [float(num) * 0.1, float(num), float(num) * 10]
        top_k = [0.9]
        batch_size = [0.01]  # [num * 4]  # [num * 4, num * 10, num * 100]
        iterations = [3]  # , 10, 100, 1000]
    else:
        rate = settings[0]
        l = settings[1]
        top_k = settings[2]
        batch_size = settings[3]
        iterations = settings[4]
        for i in range(len(iterations)):
            iterations[i] = int(iterations[i])
    accuracy = -1.0             # best accuracy on test data
    parameters = [0, 0, 0, 0, 0]   # best parameter settings
    results = []                # best results
    new_w = w                   # best weights
    for i in rate:
        for j in l:
            for k in top_k:
                for m in batch_size:
                    for n in iterations:
                        if verbose:
                            print "Linear grid search with rate = " + str(i) + ", lambda = " + str(j) + \
                                  ", top k = " + str(k) + ", mini batch size = " + str(m) + \
                                  ", number of gradient descent iterations= " + str(n)
                        d = descent(verbose, i, j, k, m, n, w, train)
                        v = validate(verbose, d, test, k)
                        if v[0] > accuracy:
                            accuracy = v[0]
                            parameters = [i, j, k, m, n]
                            results = v
                            new_w = d
    if verbose:
        print "\nBest results: "
        print "Accuracy: " + str(results[0])
        print "Precision and recall: " + str(results[1]) + " and " + str(results[2])
        print "F1 Score: " + str(results[4])
        print "True Positive Rate and False Positive Rate: " + str(results[2]) + " and " + str(results[3]) + "\n"
    r = open(res, 'w')
    r.write('[Predicted, Actual]\n')
    for line in results[5]:
        r.write(str(line) + '\n')
    r.write("\nBest weights: " + str(new_w) + "\n")
    r.write("\nBest parameter settings [rate, lambda, top_k, mini batch size, iterations]:\n" + str(parameters))
    r.write("\nBest results: \n")
    r.write("Accuracy: " + str(results[0]) + '\n')
    r.write("Precision and recall: " + str(results[1]) + " and " + str(results[2]) + '\n')
    r.write("F1 Score: " + str(results[4]) + '\n')
    r.write("True Positive Rate and False Positive Rate: " + str(results[2]) + " and " + str(results[3]) + "\n")
    r.close()
    ro = open(res + '_ordered', 'w')
    pred_act = results[5]
    pred_act.sort(key=lambda x: x[0])   # sort by predicted label
    print pred_act
    for line in range(len(pred_act)):  # loop to write actual labels in order of predicted value for that label
        ro.write(str(pred_act[len(pred_act) - line - 1][1]) + '\n')
    ro.close()

    return accuracy, parameters, results, new_w


def listify(line):      # makes a string representing a 2-d list into an actual 2-d list
    line = list(filter(None, line.split('[')))
    for i in range(len(line)):
        line[i] = line[i].strip().strip(',').strip().strip(']').split(', ')
        for j in range(len(line[i])):
            line[i][j] = float(line[i][j])
    return line


def parseargs():    # handle user arguments
    parser = argparse.ArgumentParser(description='Implementation of methods presented in GICF paper.')
    parser.add_argument('-i', '--infile', default='NONE', help='Location of input file with data')
    parser.add_argument('--limit', default=-1, help='Limit the number of bags to use. Default: No limit.')
    parser.add_argument('--out_dir', default='./', help='Directory to write files to.')
    parser.add_argument('--results', default='results', help='Where to write results to.')
    parser.add_argument('--settings', default='NONE', help='Parameter settings can be input.')
    parser.add_argument('--test', default='test.dat', help='Where to write test data to.')
    parser.add_argument('--testset', default='NONE', help='Path to file with list of patients to be put in test set.')
    parser.add_argument('--train', default='train.dat', help='Where to write training data to.')
    parser.add_argument('-v', '--verbose', default=False, action='store_true', help='Verbose output')
    args = parser.parse_args()
    return args


def main():
    start = time.time()
    args = parseargs()
    if args.infile == 'NONE':
        print "Location of input file with data must be specified."
        sys.exit()
    if args.settings != 'NONE':
        try:
            settings = listify(args.settings)
        except ValueError:
            "Error: settings input could not be parsed."
            sys.exit()
    else:
        settings = 'NONE'
    try:
        limit = int(args.limit)
    except ValueError:
        print "Error: Limit must be an integer greater than 1 to ensure non-empty training and test sets."
        sys.exit()
    if limit != -1 and limit < 2:
        print "Error: Limit must be an integer greater than 1 to ensure non-empty training and test sets."
        sys.exit()
    if not(args.out_dir.endswith('/')):
        args.out_dir += '/'
    if not(args.test.startswith('/')):
        args.test = args.out_dir + args.test
    if not(args.train.startswith('/')):
        args.train = args.out_dir + args.train
    if not(args.results.startswith('/')):
        args.results = args.out_dir + args.results
    if args.verbose:
        print "Reading input file: " + args.infile

    w = [0]
    infile = open(args.infile)     # read specified input file
    for line in infile:
        instance = listify(line.split(':')[1])[0]
        randomizer = float(len(instance) * len(instance))
        for i in range(len(instance)):
            w.append((random.random() * 2 - 1) / randomizer)     # initialize weight vector
        break
    infile.close()

    grid_search(args.verbose, w, args.infile, args.train, args.test, args.results, limit, settings, args.testset)
    print "GICF program execution time: " + str(time.time() - start) + " seconds"


if __name__ == "__main__":
    main()