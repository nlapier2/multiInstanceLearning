# Author: Nathan LaPierre
# Date: July 29, 2016

import argparse
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import os
import shlex
import subprocess
import sys
import time


def graph(title, xaxis, yaxis, names, values, results, outdir, verbose):    # graph the plots
    if verbose:
        print "Plotting graph: " + title
    plt.figure(figsize=(8, 6), dpi=100)
    plt.subplot(111)
    plt.title(title)
    plt.xlabel(xaxis[2])
    plt.xlim(xaxis[0], xaxis[1])
    plt.ylabel(yaxis[2])
    plt.ylim(yaxis[0], yaxis[1])
    for i in range(len(results)):
        plt.plot(values[i], results[i], label=names[i], linewidth=2.5, linestyle='-', marker='o')
    plt.legend(loc='upper left', frameon=False)
    plt.savefig(outdir + title + '.png')


def read_results(result_file, metric, verbose):     # read the results from the program runs
    if verbose:
        print "Reading results file " + result_file + " for metric " + metric
    f = open(result_file, 'r')
    for line in f:
        if line.startswith(metric):
            result = float(line.strip().split(metric + ': ')[1])  # retrieve the value of the desired metric
            f.close()
            return result


def execute(series, outdir, verbose):  # execute necessary program runs
    if verbose:
        print "Executing runs for " + str(series)
    resultdir = outdir + 'results/temp/'
    params = {'rate': 0, 'lambda': 1, 'top-k': 2, 'mini-batch-size': 3, 'iterations': 4}    # gicf parameters
    variable = params[series[1].split(' ')[0]]  # variable in the parameter list to vary
    values = series[1].split(' ')[1].split(',')  # possible values for the variable
    series[0] += ' --out_dir ' + resultdir  # line to run
    line = series[0]
    metrics = series[2].split(',')
    valuelist = []  # contains full list of x values needed for graphs
    results = []    # contains full list of y values needed for graphs
    for m in range(len(metrics)):
        valuelist.append([])
        results.append([])
        for v in range(len(values)):
            valuelist[m].append(float(values[v]))
    for i in range(len(values)):
        # next three lines: isolate settings array, change the variable to desired value, then place back in line
        paramlist = listify(line.split("--settings '")[1].split("'")[0])
        paramlist[variable] = [float(values[i])]
        line = series[0].split("--settings '")[0] + "--settings '" + str(paramlist) + "'" + series[0].split("'")[2]
        run = []
        for m in range(len(metrics)):
            run.append(0.0)
        for j in range(int(series[3])):
            subprocess.call(shlex.split(line))
            for m in range(len(metrics)):
                run[m] += read_results(resultdir + 'results', metrics[m], verbose)     # run result for this iteration
        for m in range(len(run)):
            results[m].append(run[m] / int(series[3]))    # average run result over all iterations for this variable
    return valuelist, results


def listify(line):      # makes a string representing a 2-d list into an actual 2-d list
    line = list(filter(None, line.split('[')))
    for i in range(len(line)):
        line[i] = line[i].strip().strip(',').strip().strip(']').split(', ')
        for j in range(len(line[i])):
            line[i][j] = float(line[i][j])
    return line


def parseargs():    # handle user arguments
    parser = argparse.ArgumentParser(description='File for plotting matplotlib graphs of GICF performance.')
    parser.add_argument('--out_dir', default='./', help='Directory to write files to.')
    parser.add_argument('--settings', default='NONE', help='Location of settings file')
    parser.add_argument('-v', '--verbose', default=False, action='store_true', help='Verbose output')
    args = parser.parse_args()
    return args


def main():
    start = time.time()
    args = parseargs()
    if not args.out_dir.endswith('/'):
        args.out_dir += '/'
    if not os.path.isdir(args.out_dir + 'results/temp'):
        if args.verbose:
            print "Attempting to make results and results/temp directories in " + args.out_dir
        os.makedirs(args.out_dir + 'results/temp')

    if args.verbose:
        print "Reading settings file: " + args.settings
    settings = open(args.settings, 'r')
    linecount = 1   # which line in the settings file we're in
    line = settings.readline()
    while line:
        title = 'notset'
        xaxis = 'notset'
        yaxis = 'notset'
        series = []
        names = []
        if len(line) > 2 and not line.startswith('#'):  # not empty line or comment
            while not line.startswith('End'):
                splits = line.strip().split('#')[0].strip().split(': ')
                try:
                    if splits[0] == 'Title':
                        if len(splits) != 2:
                            print 'Error on line ' + str(linecount) + ': malformed line'
                            sys.exit()
                        title = splits[1]
                    elif splits[0] == 'X-Axis':
                        if len(splits) != 3:
                            print 'Error on line ' + str(linecount) + ': malformed line - incorrect number of arguments'
                            sys.exit()
                        limits = splits[1].split(' to ')
                        xaxis = [float(limits[0]), float(limits[1]), splits[2]]
                    elif splits[0] == 'Y-Axis':
                        if len(splits) != 3:
                            print 'Error on line ' + str(linecount) + ': malformed line - incorrect number of arguments'
                            sys.exit()
                        limits = splits[1].split(' to ')
                        yaxis = [float(limits[0]), float(limits[1]), splits[2]]
                    elif splits[0].startswith('Series'):
                        if len(splits) != 6:
                            print 'Error on line ' + str(linecount) + ': malformed line - incorrect number of arguments'
                            sys.exit()
                        names.extend(splits[1].split(','))
                        series.append([splits[2], splits[3], splits[4], int(splits[5])])
                except ValueError:
                    print 'Error on line ' + str(linecount) + ': malformed line - ' + sys.exc_info()[0]
                    raise
                line = settings.readline()
                linecount += 1
            if title == 'notset' or xaxis == 'notset' or yaxis == 'notset' or len(series) == 0:
                print 'Error on line ' + str(linecount) + ': title, x-axis, y-axis, and series must all be defined'
                sys.exit()
            values = []
            results = []
            for i in range(len(series)):
                value, result = execute(series[i], args.out_dir, args.verbose)
                values.extend(value)
                results.extend(result)
            graph(title, xaxis, yaxis, names, values, results, args.out_dir, args.verbose)
        line = settings.readline()
        linecount += 1
    settings.close()

    if args.verbose:
        print "Program execution time: " + str(time.time() - start) + " seconds"


if __name__ == "__main__":
    main()