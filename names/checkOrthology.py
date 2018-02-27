#!/usr/bin/python3

import argparse
import glob

def parseCli():
	'''
	Parses commandline arguments
	'''
	parser = argparse.ArgumentParser()
	parser.add_argument("--directory", "-d", help="directory")
	# parser.add_argument("--", "-", help="")
	# parser.add_argument("--", "-", help="")

	return parser.parse_args()

def parseOrthoFile(filename, many2one, rest):
	'''
	parse ortholog file
	'''
	file = open(filename, "rU")
	for line in file:
		line = line.strip().split("\t")
		orthology = line[1]
		if orthology == "many2one":
			many2one.append(line[0])
		else:
			rest.append(line[0])

	file.close()
	return many2one, rest

def main():
	'''
	Controller function
	'''
	args = parseCli()
	many2one = list()
	rest = list()
	for filename in glob.glob(args.directory+"*.txt"):
		many2one, rest = parseOrthoFile(filename, many2one, rest)

	print(len(many2one))
	print(len(rest))
	dupes = list()
	for geneid in many2one:
		if geneid in rest:
			dupes.append(geneid)
	print(dupes)
	print(len(dupes))
if __name__ == "__main__":
	main()