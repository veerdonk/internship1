#!/usr/bin/python3
import argparse

from Bio import SeqIO
from Ortholog import Ortholog

def parseCli():
	'''
	Parses commandline arguments
	'''
	parser = argparse.ArgumentParser()
	parser.add_argument("--cdsOrg1", "-c1", help="CDS file for first organism")
	parser.add_argument("--cdsOrg2", "-c2", help="CDS file for second organism")
	parser.add_argument("--ensemblOrthologs", "-e", help="list of orthologs from ensembl")
	parser.add_argument("--output", "-o", help="name of output file")

	return parser.parse_args()

def parseCds(filename):
	'''
	parsed CDS files
	'''
	cds = SeqIO.parse(filename, "fasta")
	stopCodons = dict()

	for rec in cds:
		stopCodons[rec.id] = str(rec.seq[3:6])

	return stopCodons

def parseOrthologs(filename):
	'''
	parses an orthologs file from ensembl
	expected format:
	ensID \t gene name \t orthology type \t gene name \t ensID
	'''
	orthologFile = open(filename, "rU")

	orthologs = list()
	next(orthologFile)
	for line in orthologFile:
		ort1Ens, ort1Name, orthology, ort2name, ort2Ens = line.strip().split("\t")
		ortholog = Ortholog(ort1Ens, ort1Name, orthology, ort2name, ort2Ens)
		orthologs.append(ortholog)


	orthologFile.close()
	return orthologs

def writeOutput(outfile, orthologs, stopCodonsOrt1, stopCodonsOrt2):
	'''
	writes an outputfile
	'''
	# print(stopCodonsOrt1)
	# print(stopCodonsOrt2)

	out = open(outfile, "w")

	for ortholog in orthologs:
		out.write("{}\t{}\t{}\t{}\n".format(ortholog.id1, ortholog.id2, stopCodonsOrt1[ortholog.id1], stopCodonsOrt2[ortholog.id2]))

	out.close()

def main():
	'''
	Controller function
	'''
	args = parseCli()
	orthologs = parseOrthologs(args.ensemblOrthologs)
	stopCodonsOrt1 = parseCds(args.cdsOrg1)
	stopCodonsOrt2 = parseCds(args.cdsOrg2)
	writeOutput(args.output, orthologs, stopCodonsOrt1, stopCodonsOrt2)

if __name__ == "__main__":
	main()