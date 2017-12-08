#!/usr/bin/python3
import argparse

def parseCli():
	'''
	Parses commandline arguments
	'''
	parser = argparse.ArgumentParser()
	parser.add_argument("--orthologs", "-o", help="file containing orthologs to be analyzed")
	# parser.add_argument("--", "-", help="")
	# parser.add_argument("--", "-", help="")

	return parser.parse_args()

def parseOrtho(filename):
	'''
	Parses a .txt containing orthology information
	'''
	orthologRef = dict()
	orthologOrt = dict()

	orthoFile = open(filename, "rU")

	for line in orthoFile:
		line = line.split("\t")
		print(line[1])


def assignOrthologyType(orthologRefKey, orthologOrtKey):
	'''
	assigns type of orthology to orthologs
	ie 1:1 1:many many:1 many:many
	'''

def main():
	'''
	Controller function
	'''
	args = parseCli()
	parseOrtho(args.orthologs)

if __name__ == "__main__":
	main()