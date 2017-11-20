#!/usr/bin/python3

'''
Title:	upsetDataConvert.py
Author:	David van de Veerdonk
Pupose:	Convert files containing dnds ratios and gene names
		into a data format compliant with UpSetR
Date:	11/2017
'''

import argparse
import glob

def parseCli():
	'''
	Parses commandline arguments
	'''
	parser = argparse.ArgumentParser()
	parser.add_argument("--dndsDirectory", "-d", help="directory containing files to add to upsetdata")
	parser.add_argument("--knownGenes", "-k", help="csv containing known ageing genes")
	parser.add_argument("--geneNames", "-g", help="file containing all gene names of the reference organism")
	parser.add_argument("--outFile", "-o", help="name of the output file")

	return parser.parse_args()

def parseDndsFile(filename, dndsGeneDict):
	'''
	Parses a .tsv containing dnds ratios and genenames
	IN:
	dnds file with gene names
	OUT:
	filename : names 
	'''

	genes = list()
	file = open(filename, "rU")
	for line in file:
		line = line.split("\t")
		# if float(line[3]) > 1:
		genes.append(line[0])

	dndsGeneDict[filename] = genes
	file.close()

	return dndsGeneDict

def getRefGeneNames(filename):
	'''	
	Collects all genenames in a given file
	IN:
	tsv file with gene names in the 2nd column
	OUT:
	set of gene names
	'''	

	geneNames = set()
	file = open(filename, "rU")
	next(file)
	for line in file:
		line = line.strip().split("\t")
		if len(line) > 1:
			geneNames.add(line[1].split(" ")[0]) 

	file.close()

	return geneNames

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
		knownGenes.add(line[0].lower())
		if len(line) > 1:
			for alias in line[1].split(" "):
				knownGenes.add(alias.lower())

	geneFile.close()
	return knownGenes

def createUpsetData(refGeneNames, knownGeneNames, dndsGeneDict, outfileName):
	'''
	writes an output file compatible with UpSetR
	IN:
	set of reference gene names
	set of known gene names
	dict containing the gene names of psgs tied to filename
	OUT:
	csv containing containing presence of gene data
	eg:
	 		file1	file2	file3
	dpp6a	1		1		0
	zic5	0		1		1
	'''
	
	out = open(outfileName,"w")
	lineSize = len(dndsGeneDict) + 1 #number of columns (+1 for the known genes)
	out.write("gene,ageing genes,{}\n".format(dndsGeneDict.keys()))
	for gene in refGeneNames:
		lineData = [0] * lineSize
		if gene.lower() in knownGeneNames:
			lineData[0] = 1
		for i, file in enumerate(dndsGeneDict):
			if gene in dndsGeneDict[file]:
				lineData[i+1] = 1 #sets index to 1 if gene is present
		
		out.write("{},{}\n".format(gene, str(lineData)[1:-1]))
	out.close()

def main():
	'''
	Controller function
	'''
	args = parseCli()
	dndsGeneDict = dict()
	files = glob.glob(args.dndsDirectory + "*.tsv") #collects all filenames ending in .tsv in given dir
	refGeneNames = getRefGeneNames(args.geneNames)
	knownGeneNames = getKnownGenes(args.knownGenes)
	print(len(knownGeneNames))
	for file in files:
		dndsGeneDict = parseDndsFile(file, dndsGeneDict)

	createUpsetData(refGeneNames, knownGeneNames, dndsGeneDict, args.outFile)


if __name__ == "__main__":
	main()