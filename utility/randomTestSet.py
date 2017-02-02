# Author: Nathan LaPierre
# Date: August 23, 2016

import argparse
import random
import sys


def parseargs():    # handle user arguments
    parser = argparse.ArgumentParser(description='Script that picks a random test set from a set.')
    parser.add_argument('--map', default='NONE', help='File number to start on. Required.')
    parser.add_argument('--size', default=23, help='Size of test set to pick.')
    args = parser.parse_args()
    return args


def main():
    args = parseargs()
    if args.map == 'NONE':
        print 'Map file must be specified with --map. Aborting...'
        sys.exit()
    try:
        size = int(args.size)
    except ValueError:
        print 'Size must be a positive integer. Aborting...'
        sys.exit()
    if size < 1:
        print 'Size must be a positive integer. Aborting...'
        sys.exit()

    tags = []        # holds patient tags from map file
    testnums = []   # holds the line numbers in the map to choose from
    mapfile = open(args.map, 'r')
    for line in mapfile:
        tags.append(line.split('\t')[3])
    mapfile.close()
    while len(testnums) < size:
        num = random.random() * len(tags)   # pick random number in range of number of patients
        if num not in testnums:
            testnums.append(tags[int(num)])    # keep adding unique numbers until we have desired number of test items
    for t in testnums:
        print t

if __name__ == "__main__":
    main()