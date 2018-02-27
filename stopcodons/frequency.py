#!/usr/bin/python3

import argparse

from Bio import SeqIO
from collections import Counter

def parseCli():
	'''
	Parses commandline arguments
	'''

	parser = argparse.ArgumentParser()
	parser.add_argument("--inFile", "-i", help="input file")

	return parser.parse_args()

def parseCds(filename):
	'''
	parsed CDS files
	'''
	legalStops = ["TAA", "TAG", "TGA"]
	cds = SeqIO.parse(filename, "fasta")
	stopCodons = list()

	for rec in cds:
		rec.id = rec.id.split(".")[0]
		stop = str(rec.seq[-3:])
		if stop in legalStops:
			stopCodons.append(stop)

	return stopCodons


def main():
	'''
	Controller function
	'''
	
	args = parseCli()
	stops = parseCds(args.inFile)
	counts = Counter(stops)
	lenAllStops = len(stops)

	print("TGA:\t{}%\nTAA:\t{}%\nTAG:\t{}%".format(round(counts["TGA"]/lenAllStops*100, 2), round(counts["TAA"]/lenAllStops*100, 2), round(counts["TAG"]/lenAllStops*100, 2)))
	
	
	



if __name__ == "__main__":
	main()