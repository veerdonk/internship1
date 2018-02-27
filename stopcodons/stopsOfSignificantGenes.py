#!/usr/bin/python3

import argparse
import glob
from Bio import SeqIO

def parseCli():
	'''
	Parses commandline arguments
	'''
	parser = argparse.ArgumentParser()
	parser.add_argument("--cdsDir", "-c", help="directory containing all CDS files")
	parser.add_argument("--namesDir", "-n", help="directory containing all name files")
	parser.add_argument("--sigGenes", "-s", help="significant geneslist")
	parser.add_argument("--outfile", "-o", help="name of output file")

	return parser.parse_args()

def parseCds(filename):
	'''
	parsed CDS files
	'''
	cds = SeqIO.parse(filename, "fasta")
	stopCodons = dict()
	legalStops = {"TGA", "TAA", "TAG"}

	for rec in cds:
		rec.id = rec.id.split(".")[0]
		stop = str(rec.seq[-3:])
		if stop in legalStops:
			stopCodons[rec.id] = stop

	return stopCodons

def parseNamesFile(filename):
	'''
	parse file containing gene names
	'''
	nameMapping = dict()
	namesFile = open(filename, "rU")
	for line in namesFile:
		line = line.strip().split("\t")
		if line[0] != "NA":
			nameMapping[line[0]] = [line[2].split(".")[0], line[5]]

	return nameMapping

def parseSigGenes(filename):
	'''
	parses file with significant genes
	'''
	sigGenes = dict()
	sigFile = open(filename ,"rU")
	next(sigFile)
	for line in sigFile:
		line = line.strip().split(",")
		if len(line) > 1:
			sigGenes[line[0][1:-1]] = line[1]
		else:
			sigGenes[line[0][1:-1]] = 0

	return sigGenes

def main():
	'''
	Controller function
	'''
	args = parseCli()
	outStops = open(args.outfile + "stops.csv", "w")
	outDnds = open(args.outfile + "dnds.csv", "w")

	cdsFiles = glob.glob(args.cdsDir + "/*")
	cdsFilenameDict = dict()
	for cdsFilename in cdsFiles:
		cdsFilenameDict[cdsFilename.split("/")[-1].split(".")[0]] = cdsFilename

	sigGenes = parseSigGenes(args.sigGenes)
	for gene in sigGenes:
		outStops.write(",{}".format(gene))
		outDnds.write("{},".format(gene))

	for filename in glob.glob(args.namesDir + "/*"):
		cdsFile = cdsFilenameDict[filename.split("/")[-1].split("_")[1]]
		outStops.write("\n{},".format(cdsFile.split("/")[-1].split(".")[0]))
		outDnds.write("\n{},".format(cdsFile.split("/")[-1].split(".")[0]))
		print("parsing CDS file: {}".format(cdsFile.split("/")[-1]))
		cds = parseCds(cdsFile)
		names = parseNamesFile(filename)
		for gene in sigGenes:
			try:
				outStops.write("{},".format(cds[names[gene][0]]))
				# print(names[gene][1])
				if float(names[gene][1]) > 0.001 and float(names[gene][1]) < 10:
					outDnds.write("{},".format(names[gene][1]))
				else:
					outDnds.write(",")
			except:
				outStops.write(",")
				outDnds.write(",")

	outStops.close()
	outDnds.close()

if __name__ == "__main__":
	main()