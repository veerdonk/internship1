#!/usr/bin/python3

import argparse

from Bio import SeqIO

def parseCli():
	'''
	Parses commandline arguments
	'''
	parser = argparse.ArgumentParser()
	parser.add_argument("--cds", "-c", help="CDS file")
	parser.add_argument("--dndsFile", "-d", help="File containing dnds values/gene names")
	parser.add_argument("--out", "-o", help="outfile")

	return parser.parse_args()

def getStopCodons(cdsFilename):
	'''
	Retrieves stopcodons from cds
	IN:
	CDS file
	OUT:
	id : stopcodon map
	'''
	cds = SeqIO.parse(cdsFilename, "fasta")
	stopCodons = dict()

	for rec in cds:
		stopCodons[rec.id] = str(rec.seq[-3:])

	return stopCodons

def parseDnds(dndsFilename):
	'''
	Parses a file comtaininge dndsvalues and gene 
	names, etc
	returns dict with id : data
	'''
	dndsFile = open(dndsFilename, "rU")
	dnds = dict()

	for line in dndsFile:
		line = line.strip().split("\t")
		dnds[line[2]] = line

	return dnds

def writeOut(filename, stops, dnds):
	'''
	writes outputfile
	'''
	out = open(filename, "w")

	for val in dnds:
		out.write("{},{}\n".format(str(dnds[val])[1:-1].replace("'", ""), stops[val]))


def main():
	'''
	Controller function
	'''
	args = parseCli()
	stops = getStopCodons(args.cds)
	dnds = parseDnds(args.dndsFile)
	writeOut(args.out, stops, dnds)

if __name__ == "__main__":
	main()