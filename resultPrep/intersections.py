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
	parser.add_argument("--upsetDataFile", "-u", help="if given an UpSetR compatible .csv will be written.")

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

def writeUpset(filename, geneOccurance, organisms):
	'''
	writes an UpSetR compatible file using the known gene occurance rates
	'''
	out = open(filename, "w")
	out.write("gene,{}\n".format(str(organisms)[1:-1].replace("'", "")))

	for gene in geneOccurance:
		line = [0]*len(organisms)
		for org in geneOccurance[gene]:
			line[organisms.index(org)] = 1
		out.write("{},{}\n".format(gene, str(line)[1:-1]))

	out.close()

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
		allPsg[name] = psg
		print(sorted(intersect))
		# print("Known ageing genes found in {}: {}\nOf these {} genes {} were positively selected for".format(filename, len(intersect), len(intersect), len(psg)))
	
	print("\nKnown ageing genes found in more than 2 organisms:")
	organisms = list()
	geneList = open("geneTestFile.csv", "w")
	geneList.write("Ageing genes\ngene,number,org\n")
	for gene in geneOccurance:
		for org in geneOccurance[gene]:
			if org not in organisms:
				organisms.append(org)
		if len(geneOccurance[gene]) > 2:
			print("{}\t\t: {}/{} -> {}".format(gene, len(geneOccurance[gene]), len(filenames), sorted(geneOccurance[gene])))
			geneList.write("{},{}/{},{}\n".format(gene, len(geneOccurance[gene]), len(filenames), str(sorted(geneOccurance[gene])).replace("'", "")[1:-1]))

	if args.upsetDataFile is not None:
		writeUpset(args.upsetDataFile, geneOccurance, organisms)
		psgNames = dict()
		psgOrgs = list()
		for org in allPsg:
			for gene in allPsg[org]:
				if gene not in psgNames:
					psgNames[gene] = [org]
				else:
					psgNames[gene].append(org)
			if org not in psgOrgs:
				psgOrgs.append(org)

	
		writeUpset(args.upsetDataFile[:-4]+"_allPSG.csv", psgNames, psgOrgs)

	print("\nPositively selected genes found in more than 10 organisms:")
	geneList.write("\n\nALL genes\ngene,number,org\n")
	for gene in psgNames:
		if len(psgNames[gene]) > 9:
			print("{}\t: {}/{} -> {}".format(gene, len(psgNames[gene]), len(filenames), sorted(psgNames[gene])))
			geneList.write("{},{}/{},{}\n".format(gene, len(psgNames[gene]), len(filenames), str(sorted(psgNames[gene])).replace("'", "")[1:-1]))
	
	geneList.close()
if __name__ == "__main__":
	main()