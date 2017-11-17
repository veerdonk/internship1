#!/usr/bin/python3

import argparse
import glob

def parseCli():
	'''
	Parses commandline arguments
	'''
	parser = argparse.ArgumentParser()
	parser.add_argument("--inf", "-i", help="indir")
	parser.add_argument("--ef", "-e", help="ensfile")
	parser.add_argument("--ag", "-g", help="ageing genes")

	return parser.parse_args()

def getKnownGenes(filename):
	'''
	retrieves a set of known genes and their 
	aliases from a file
	IN:
	file containing gene names to intersect
	OUT:
	set containing gene names
	'''
	knownGenes = set()
	geneFile = open(filename, "rU")
	for line in geneFile:
		line = line.split(",")
		knownGenes.add(line[1])
		for alias in line[2].split(" "):
			knownGenes.add(alias)

	return knownGenes

def main():
	'''
	Controller function
	'''

	#TODO clean this up/separate into functions
	args = parseCli()

	files = glob.glob(args.inf+"*.tsv")


	knownGenes = getKnownGenes(args.ag)

	geneNames = set()
	ensFile = open(args.ef, "rU")
	namesForFile = dict()
	for line in ensFile:
		line = line.strip().split("\t")
		if len(line) > 1:
			geneNames.add(line[1])
			namesForFile[line[1]] = [0]*(len(files)+1)
		else:
			geneNames.add(line[0])
			namesForFile[line[0]] = [0]*(len(files)+1)
	ensFile.close()

	for gene in knownGenes:
		if gene in namesForFile:
			namesForFile[gene][0] = 1

	i = 1
	for file in files:
		curSet = set()
		data = open(file, "rU")
		for line in data:
			geneName, ref, ort, dnds, dn, ds = line.split("\t")
			if float(dnds) > 1:
				curSet.add(geneName)
				namesForFile[geneName][i] = 1
		geneNames = curSet & geneNames
		i += 1
		# print(len(curSet & geneNames))
	
	out = open("test.csv", "w")
	out.write("file, knownGenes, {}\n".format(str(files)[1:-1]))
	for name in namesForFile:
		out.write("{},{}\n".format(name, str(namesForFile[name])[1:-1]))

if __name__ == "__main__":
	main()