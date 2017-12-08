#!/usr/bin/python3
import argparse
import re

def parseCli():
	'''
	Parses commandline arguments
	'''
	parser = argparse.ArgumentParser()
	parser.add_argument("--names", "-n", help="File containing id - gene name mapping")
	parser.add_argument("--CDS", "-c", help="CDS file")
	parser.add_argument("--out", "-o", help="name of gene name file with flags")

	return parser.parse_args()

def parseCDS(filename):
	'''
	Parses a CDS file and counts the
	number of nucleotides per transcript
	IN:
	CDS
	OUT:
	dict containing id : length mapping
	'''
	cdsFile = open(filename, "rU")
	cdsLength = dict()
	lineLen = 0

	for line in cdsFile:
		if line.startswith(">"):
			if lineLen != 0:
				cdsLength[ID] = lineLen
			lineLen = 0
			m = re.search("_([a-zA-Z0-9]+)_\d|(Cau_\w+)", line)
			if m:
				if m.group(1) is not None:
					ID = m.group(1)
				else:
					ID = m.group(2)
		else:
			lineLen += len(line)
	cdsLength[ID] = lineLen

	return(cdsLength)

def geneOccurance(filename):
	'''
	reads a file containing id : genename
	mapping to find occurance numbers
	IN:
	id \t name file
	OUT:
	dict with id : occurance
	'''
	nameFile = open(filename, "rU")
	occurance = dict()

	for line in nameFile:
		orgid, gene = line.strip().split("\t")

		if gene in occurance:
			occurance[gene].append(orgid)
		else:
			occurance[gene] = [orgid]

	return occurance


def writeOutFlags(filename, cds, occurance):
	'''
	writes an output file based on gene
	occurance and cds length with removal
	flags set for shorter transcripts
	'''
	out = open(filename, "w")
	out.write("identifier\tgenename\n")
	for gene in occurance:
		if len(occurance[gene]) > 1:
			longest = 0
			for name in occurance[gene]:
				if cds[name] > longest:
					longestTranscript = name
					longest = cds[name]
			for name in occurance[gene]:
				if name == longestTranscript:
					out.write("{}\t{}\n".format(name, gene))
				else:
					out.write("{}\tREMOVE\n".format(name))
		else:
			out.write("{}\t{}\n".format(occurance[gene][0], gene))

	out.close()


def main():
	'''
	Controller function
	'''
	args = parseCli()
	cds = parseCDS(args.CDS)
	occurance = geneOccurance(args.names)
	writeOutFlags(args.out, cds, occurance)

if __name__ == "__main__":
	main()