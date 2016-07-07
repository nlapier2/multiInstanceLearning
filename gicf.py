# Author: Nathan LaPierre
# Date: June 6, 2016


import argparse
import math
import random
import sys
import time


def similarity(i, j):   # similarity measure (0 to 1) between instances i and j
    norm = 0.0    # the squared 2-norm of the difference between the two instance vectors, computed in for loop below
    for k in range(len(i)):
        norm += (i[k] - j[k]) ** 2
    return math.exp(-norm)      # returns value for kernel described in GICF paper


def instance_penalty(i, j):     # non-negative penalty on the difference between predictions for instances i and j
    return (i - j) ** 2     # using squared loss as penalty, as described in GICF paper


def logistic(w, x):     # implementation of the classifier for the instances; currently logistic regression
    dot = 0.0     # dot product between weights and x, computed in for loop below
    for i in range(len(x)):
        dot += w[i]*x[i]
    if dot > 20:    # prevent math overflows
        return 1
    if dot < -20:   # prevent math overflows
        return 0
    return ((1 / (1 + math.exp(-dot))) - 0.5) * 2  # return normalized logistic regression formula using dot


def derivative_logistic(w, x):   # derivative of logistic regression function
    dot = 0.0     # dot product between weights and x, computed in for loop below
    for i in range(len(x)):
        dot += w[i]*x[i]
    if dot > 20 or dot < -20:    # prevent math overflows
        return 0
    e = math.exp(dot)   # e ^ (dot product of w and x)
    norm = 0.0        # norm of vector x times scalar e, calculated below
    for k in range(len(x)):
        norm += (x[k] * e) ** 2
    norm = math.sqrt(norm)
    value = math.fabs(norm / ((1 + e) ** 2))    # value for derivative of logistic function
    return ((1 / (value + 1)) - 0.5) * 2    # normalized result


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
    total1 = 0.0    # the part of the group penalty using only predicted instance labels
    total2 = 0.0    # the part of the group penalty using the actual group label
    for i in range(len(group)):
        total1 += (instance_labels[i] * d_instance_labels[i])
        total2 += d_instance_labels[i]
    total1 = total1 * 2.0 / len(group)
    total2 = total2 * 2.0 * group_label / len(group)
    return total1 - total2


def cost(num_instances, w, data, l, group_labels):
    # calculate cost function as defined in GICF paper
    instance_cost = 0.0     # first part of cost function, based on instance differences
    group_cost = 0.0    # second part of cost function, based on group labels
    # l = float(len(data))  # "lambda" in GICF cost function, to be changed later

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
    # print instance_cost
    # print group_cost
    instance_cost /= (num_instances ** 2)   # average out instance cost across all instance pairs
    group_cost *= (l / float(len(data)))    # balance group cost using lambda parameter

    return instance_cost + group_cost


def derivative_cost(num_instances, w, data, l, group_labels):      # derivative of cost function above
    instance_cost = 0.0     # first part of cost function, based on instance differences
    group_cost = 0.0    # second part of cost function, based on group labels

    instance_labels = []    # labels of instances according to weights and logistic function
    d_instance_labels = []  # labels of instances according to weights and derivative of logistic function
    for i in range(len(data)):
        instance_labels.append([])
        d_instance_labels.append([])
        for j in range(len(data[i])):
            instance_labels[i].append(logistic(w, data[i][j]))
            d_instance_labels[i].append(derivative_logistic(w, data[i][j]))

    for group in range(len(data)):
        group_cost += d_group_penalty(data[group], instance_labels[group],
                                      d_instance_labels[group], group_labels[group])
        for i in range(len(data[group])):
            for j in range(len(data[group])):
                instance_cost += similarity(data[group][i], data[group][j]) * \
                    (instance_labels[group][i] - instance_labels[group][j]) * \
                    (d_instance_labels[group][i] - d_instance_labels[group][j])
    instance_cost /= (num_instances ** 2)   # average out instance cost across all instance pairs
    group_cost *= (l / float(len(data)))    # balance group cost using lambda parameter
    return instance_cost + group_cost


def descent(verbose, rate, l, batch_size, iterations, w, count, data):
    # mini-batch stochastic gradient descent to learn weights for classifier
    # randomly select a small subset of instances from data to be used as mini batch
    if verbose:
        print "Entering gradient descent..."
    per_group = batch_size / count    # amount of instances to use per group
    threshold = 1    # when slope crosses this threshold, stop
    final_cost = 0  # the cost function value when gradient descent has finished
    update_vector = []  # this is used to track the previous update vector, used for momentum
    momentum = 0.5      # constance multiplied by the previous update vector for momentum
    for i in range(len(w) - 1):
        update_vector.append(0)  # initialize to same length as instance, with all 0s
    for num in range(iterations):
        if verbose:
            print "Gradient descent, instance " + str(num + 1)
        cost_data = []  # the data to be given to the cost function
        average_instance = []   # the average value of the instances in cost_data, used for weight adjustment
        labels = []     # group labels
        for i in range(len(w) - 1):
            average_instance.append(0)  # initialize to same length as instance, with all 0s
        d = open(data, 'r')
        for line in d:
            if len(line) < 5:   # deals with blank lines in files
                continue
            cost_data.append([])
            instances = []  # instances from this group to use
            splits = line.split(':')
            labels.append(int(float(splits[0])))
            group = listify(splits[1])
            for i in range(int(per_group)):
                num = int(random.random() * (len(group) - 1))     # randomly select index of instance from group
                while num in instances:     # make sure we're picking a unique instance
                    num = int(random.random() * (len(group) - 1))
                instances.append(num)
            for j in range(len(group)):
                if j in instances:
                    cost_data[len(labels)-1].append(group[j])   # append appropriate instance to mini batch dataset
                    for k in range(len(group[j])):
                        average_instance[k] += group[j][k]      # add this instance to our total
        for i in range(len(average_instance)):
            average_instance[i] /= (per_group * count)  # take the average of the instance totals
        d.close()

        # find value of derivative of cost function using weights w, then use learning rate to find adjustment amount
        delta_j = rate * derivative_cost(per_group * count, w, cost_data, l, labels)
        # delta_j = rate * cost(per_group * count, w, cost_data, l, labels)
        # once slope becomes non-negative (or barely negative), stop; or if we exceed number of iterations
        if verbose:
            print "Slope of cost function calculated: " + str(delta_j / rate)
        if math.fabs(delta_j / rate) <= threshold:
            final_cost = cost(per_group * count, w, cost_data, l, labels)
            break
        # adjust weights according to slope of cost function
        w[0] -= delta_j
        for i in range(len(average_instance)):
            # adjust weights using delta_j, instance values, momentum
            adjustment = (delta_j * average_instance[i]) + (update_vector[i] * momentum)
            w[i+1] -= adjustment
            update_vector[i] = adjustment   # update the update vector for future momentum
        # for [w0 w1 w2] & [x1 x2], w0 = w0 - rate*deltaJ; w1 = w1 - rate*x1*deltaJ; w2 = w2 - rate*x2*deltaJ
        if num == iterations - 1:   # we are at the end of the loop
            final_cost = cost(per_group * count, w, cost_data, l, labels)

    if verbose:
        print "Final cost calculated: " + str(final_cost)
        # print "With weights: " + str(w)
    return final_cost, w


def validate(verbose, w, data):     # apply the weights determined by gradient descent to a test set
    if verbose:
        print "Validating against test set..."
    predicted_actual = []   # array containing predicted and actual labels for each group
    d = open(data, 'r')
    for line in d:
        if len(line) < 5:   # deals with blank lines in files
            continue
        splits = line.split(':')
        label = int(float(splits[0]))
        group = listify(splits[1])
        total = 0.0     # total of instance label predictions
        for instance in group:
            total += logistic(w, instance)
        total /= len(group)     # group label prediction is average of instance label predictions
        predicted_actual.append([total, label])     # predicted vs actual label
    d.close()

    # threshold = 0
    average_predicted = 0.0
    for i in predicted_actual:
        average_predicted += i[0]
    average_predicted /= len(predicted_actual)
    threshold = average_predicted
    print "Average prediction: " + str(average_predicted)

    true_positives = 0
    false_positives = 0
    true_negatives = 0
    false_negatives = 0
    for i in predicted_actual:
        if i[0] >= threshold and i[1] == 1:
            true_positives += 1
        elif i[0] > threshold and i[1] == -1:
            false_positives += 1
        elif i[0] < threshold and i[1] == -1:
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
        '''print "[Predicted, Actual]"
        for line in predicted_actual:
            print str(line)'''
        print "Accuracy: " + str(accuracy)
        print "Precision and recall: " + str(precision) + " and " + str(tpr_recall)
        print "F1 Score: " + str(f1_score)
        print "True Positive Rate and False Positive Rate: " + str(tpr_recall) + " and " + str(fpr)

    return accuracy, precision, tpr_recall, fpr, f1_score, predicted_actual


def grid_search(verbose, w, data, train, test, res):
    # linear grid search to set learning rate, lambda, and mini batch size parameters
    if verbose:
        print "Entering grid search to select best parameters for gradient descent."
    test_split = 2   # 1 out of every test_split instances are put in test set (i.e. 1/2 if test_split = 2)

    f = open(data, 'r')
    tr = open(train, 'w')
    te = open(test, 'w')
    num = 0     # count of the number of groups
    for line in f:
        num += 1
        if num % test_split == 0:
            te.write(line + '\n')
        else:
            tr.write(line + '\n')
        if num == 30:
            break
    f.close()
    tr.close()
    te.close()

    rate = [0.0001]  # , 0.01]  # , 0.1]
    l = [float(num)]  # , float(num), float(num)*2]
    batch_size = [num*4]  # , num*4]  # , num_groups*100]
    iterations = [30]  # , 100]  # , 10000]
    accuracy = -1.0             # best accuracy on test data
    parameters = [0, 0, 0, 0]   # best parameter settings
    results = []                # best results
    new_w = w                   # best weights
    for i in rate:
        for j in l:
            for k in batch_size:
                for h in iterations:
                    if verbose:
                        print "Linear grid search with rate = " + str(i) + ", lambda = " + str(j) + \
                              ", mini batch size = " + str(k) + ", number of gradient descent iterations = " + str(h)
                    d = descent(verbose, i, j, k, h, w, num, train)
                    v = validate(verbose, d[1], test)
                    if v[0] > accuracy:
                        accuracy = v[0]
                        parameters = [i, j, k, h]
                        results = v
                        new_w = d[1]
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
    r.write("Best weights: " + str(new_w) + "\n")
    r.write("Best parameter settings [learning rate, lambda, mini batch size, descent iterations]: " + str(parameters))
    r.write("\n\n\n\nBest results: \n")
    r.write("Accuracy: " + str(results[0]) + '\n')
    r.write("Precision and recall: " + str(results[1]) + " and " + str(results[2]) + '\n')
    r.write("F1 Score: " + str(results[4]) + '\n')
    r.write("True Positive Rate and False Positive Rate: " + str(results[2]) + " and " + str(results[3]) + "\n")
    r.close()

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
    parser.add_argument('--out_dir', default='./', help='Directory to write files to.')
    parser.add_argument('--results', default='results', help='Where to write results to.')
    parser.add_argument('--test', default='test.dat', help='Where to write test data to.')
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
        for i in range(len(instance)):
            w.append(0)     # initialize weight vector, with length of instance (+ 1 for bias)
        break
    infile.close()

    grid_search(args.verbose, w, args.infile, args.train, args.test, args.results)
    # grid_search(args.verbose, w, labels, data)

    print "Program execution time: " + str(time.time() - start) + " seconds"


if __name__ == "__main__":
    main()