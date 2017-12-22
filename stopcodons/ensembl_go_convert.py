#!/usr/bin/python3

import argparse

def parseCli():
	'''
	Parses commandline arguments
	'''
	parser = argparse.ArgumentParser()
	parser.add_argument("--input", "-i", help="File containing geneID > GO mapping")
	parser.add_argument("--out", "-o", help="File containg geneID [GO ids]")
	# parser.add_argument("--", "-", help="")

	return parser.parse_args()

def parseEnsemblGoFile(filename):
	'''
	parses a file containg go enrichment terms
	returns a dict with geneID : GO1, GO2
	'''

	gofile = open(filename, "rU")
	go = dict()

	for line in gofile:
		line = line.strip().split("\t")
		if len(line) > 1:
			geneID, goID = line
			if geneID in go:
				go[geneID].append(goID)
			else:
				go[geneID] = [goID]

	gofile.close()

	return(go)

def writeFile(filename, go):
	'''
	writes output
	'''
	out = open(filename, "w")

	for ensid in go:
		out.write("{}\t{}\n".format(ensid, str(go[ensid])[1:-1].replace("'","")))


def main():
	'''
	Controller function
	'''
	args =  parseCli()
	go = parseEnsemblGoFile(args.input)
	writeFile(args.out, go)


if __name__ == "__main__":
	main()