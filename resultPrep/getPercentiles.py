#!/usr/bin/python3
'''
Title	:	getPercentiles.py
Author	:	David van de Veerdonk
Purpose	:	Calculate mean and median 
			percentiles for each gene 
			in a group
'''

import argparse
import glob

def parseCli():
	'''
	Parses commandline arguments
	'''
	parser = argparse.ArgumentParser()
	parser.add_argument("--dir", "-d", help="directory with files to analyze")
	# parser.add_argument("--", "-", help="")
	# parser.add_argument("--", "-", help="")

	return parser.parse_args()

def parseDndsFile(filename):
	'''

	'''
	genes = dict()
	file = open(filename, "rU")
	for line in file:
		line = line.split("\t")
		geneName, ref, ort, dnds, dn, ds = line.split("\t")
		entry = dNdSEntry(ref, ort, float(dnds), float(dn), float(ds), geneName)

	return genes

def getPercentile():
	'''

	'''
	pass

def main():
	'''
	Controller function
	'''
	args = parseCli()
	filenames = glob.glob(args.dir+"*.tsv")
	allFiles = list()
	for filename in filenames:
		allFiles.append(parseDndsFile(filename))


if __name__ == "__main__":
	main()