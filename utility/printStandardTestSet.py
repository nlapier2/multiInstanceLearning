import argparse


def parseargs():    # handle user arguments
    parser = argparse.ArgumentParser(description='Script that performs assembly on a number of patient files.')
    parser.add_argument('--input', default='fakemap.txt', help='Input mapping file to generate test set from.')
    parser.add_argument('--output', default='testset.txt', help='Where to write testset output to.')
    args = parser.parse_args()
    return args


def main():
	args = parseargs()
	infile = open(args.input, 'r')
	outfile = open(args.output, 'w')
	counter = 0
	for line in infile:
		counter = (counter + 1) % 2
		if counter == 0:
			outfile.write(line.split('\t')[3] + '\n')
	infile.close()
	outfile.close()


if __name__ == '__main__':
	main()