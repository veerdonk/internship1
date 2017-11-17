#!/usr/bin/python3
'''
Title	:	intersections.py
Author	:	David van de Veerdonk
Purpose	:	Check what genes are present in 
			a file and all organisms orthologs
			present in a given directory
'''

import argparse
import glob

def parseCli():
	'''
	Parses commandline arguments
	'''
	parser = argparse.ArgumentParser()
	parser.add_argument("--dndsWithGeneNameDir", "-d", help="directory containing output of geneNames.py", required = True)
	parser.add_argument("--knownGenes", "-k", help="file containing known ageing genes", required = True)
	parser.add_argument("--genesToCheck", "-g", help="comma separated list of genes to check the presence of eg CEPBa,SKAP2,ANK2")

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
		knownGenes.add(line[0])
		if len(line) > 1:
			for alias in line[1].split(" "):
				knownGenes.add(alias)

	return knownGenes

def findIntersect(filename, interset, genesToCheck):
	'''
	finds the intersection between the genes in a 
	given file and a set of gene names
	IN:
	name of dnds file
	set of gene names to check
	OUT:
	set containing the intersection of known genes
	and the genes from a dnds file
	'''

	psgCount = 0
	totalCount = 0
	intersect = set()
	psg = set()
	dndsGeneFile = open(filename, "rU")
	for line in dndsGeneFile:
		geneName, ref, ort, dnds, dn, ds = line.split("\t")
		# if geneName.lower() in genesToCheck:
			# print("{} contains {}".format(filename, geneName))
		totalCount += 1
		if float(dnds) > 1:
			psgCount += 1

		if geneName in interset:
			intersect.add(geneName)
			if float(dnds) > 1 and float(dnds) < 99:
				psg.add(geneName)
				
	print("{} contains {} dnds ratios >1 out of a total of {} genes ({}%)".format(filename, psgCount, totalCount, round(psgCount/totalCount*100, 2)))	
	
	return psg, intersect

def main():
	'''
	Controller function
	'''
	args = parseCli()
	knownGenes = getKnownGenes(args.knownGenes)
	if args.genesToCheck != None:
		geneQuery = args.genesToCheck.lower().split(",")
		print(geneQuery)
	else:
		geneQuery = list()
	# intersect = knownGenes
	# allIntersect = knownGenes

	for filename in glob.glob(args.dndsWithGeneNameDir+"*.tsv"):
		psg, intersect = findIntersect(filename, knownGenes, geneQuery)
		# psg, allIntersect = findIntersect(filename, allIntersect)
		print("Known ageing genes found in {}: {}\nOf these {} genes {} were positively selected for".format(filename, len(intersect), len(intersect), len(psg)))
	# for gene in allIntersect:
		# print(gene)


if __name__ == "__main__":
	main()