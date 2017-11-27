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
import re

from dNdSEntry import dNdSEntry

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
		line = line.strip().split(",")
		knownGenes.add(line[0].lower())
		if len(line) > 1:
			for alias in line[1].split(" "):
				knownGenes.add(alias.lower())

	return knownGenes

def findIntersect(filename, knownGenes, genesToCheck):
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

	psg = set()
	intersect = set()
	test = dict()
	total = 0
	dndsGeneFile = open(filename, "rU")
	checked = list()

	for line in dndsGeneFile:
		total += 1
		geneName, ref, ort, dnds, dn, ds = line.split("\t")
		entry = dNdSEntry(ref, ort, float(dnds), float(dn), float(ds), geneName)
		if entry.dnds > 1 and entry.dnds < 99:
			psg.add(entry.name.lower())
			if entry.name.lower() in knownGenes:
				intersect.add(entry.name)
			if entry.name.lower() in genesToCheck:
				checked.append(entry)

	print("\n{}:\nout of {} genes {} have a dnds ratio of >1. That is {}% of the total.\nOut of these {} also apear in the known genes set.".format(filename.split("/")[-1], total, len(psg), round(len(psg)/total*100, 3), len(intersect)))

	return psg, intersect, checked


def main():
	'''
	Controller function
	'''
	args = parseCli()
	knownGenes = getKnownGenes(args.knownGenes)
	allPsg = dict()

	geneOccurance = dict()
	
	if args.genesToCheck != None:
		geneQuery = args.genesToCheck.lower().split(",")
		print("looking for: {}".format(str(args.genesToCheck.split(","))[1:-1]))
	else:
		geneQuery = list()
	filenames = glob.glob(args.dndsWithGeneNameDir+"*.tsv")
	print("analyzing {} files..".format(len(filenames)))
	for filename in filenames:
		psg, intersect, checked = findIntersect(filename, knownGenes, geneQuery)

		name = re.search("_([a-z]{3}[a-z]?)_", filename).group(1)
		for gene in intersect:
			if gene not in geneOccurance:
				geneOccurance[gene] = [name]
			else:
				geneOccurance[gene].append(name)

		if len(checked) > 0:
			print("checked genes found in file:")
			for entry in checked:
				print(entry.name)
		allPsg[filename.split("/")[-1]] = psg
		print(sorted(intersect))
		# print("Known ageing genes found in {}: {}\nOf these {} genes {} were positively selected for".format(filename, len(intersect), len(intersect), len(psg)))
	
	old = set(("hgl", "can", "fda", "cla", "ode"))


	for gene in geneOccurance:
		if len(geneOccurance[gene]) > 2:
			print(set(geneOccurance[gene]) & old)
			print("{} : {} -> {}".format(gene, len(geneOccurance[gene]), geneOccurance[gene]))


	# oldIntersect = allPsg["hsa_ggo_dndsGeneNames.tsv"] & allPsg["hsa_pan_dndsGeneNames.tsv"] & allPsg["hsa_mmu_dndsGeneNames.tsv"]
	# youngIntersect = allPsg["hsa_soe_dndsGeneNames.tsv"] & allPsg["hsa_csa_dndsGeneNames.tsv"] & allPsg["hsa_tsy_dndsGeneNames.tsv"]
	# print(len(oldIntersect))
	# print(len(youngIntersect))
	# print(youngIntersect)

if __name__ == "__main__":
	main()