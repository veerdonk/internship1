#!/usr/bin/python3
import argparse
import glob
import re

from dNdSEntry import dNdSEntry

def parseCli():
	'''
	Parses commandline arguments
	'''
	parser = argparse.ArgumentParser()
	parser.add_argument("--dndsDir", "-d", help="directory containing output of getDnds.py that you want to assign reference gene names")
	parser.add_argument("--geneNames", "-g", help="file containing matching identifiers with dndsFile and gene names")
	parser.add_argument("--outputDir", "-o", help="directory where you want the output saved")

	return parser.parse_args()

def getNamesFromFile(filename):
	'''
	parses the geneName file and retrieves a list of genes
	IN:
	genename file (either gff or ensembl txt)
	needs the extensions to be correct
	OUT:
	dictionary with ID and genename eg: ENS0003121 : SKAP
	'''

	if filename[-3:] == "gff":
		print("GFF detected, parsing..")
		geneNames = parseGffGenes(filename)

	elif filename[-3:] == "txt":
		print("txt detected, parsing..")
		geneNames = parseEnsemblGenesFile(filename)

	elif filename[-3:] == "cds":
		print("cds detected, parsing..")
		geneNames = parseGeneNamesInName(filename)

	elif filename[-3:] == "tsv":
		print("tsv detected, parsing..")
		geneNames = parseGeneNamesInExcel(filename)

	return geneNames

def parseGeneNamesInExcel(filename):
	'''
	opens and parses an excel file containing identifiers and gene names
	NOTE: in format has to be this identifier \t homolog <whitespace> other stuff
	eg bmy_00001T0 \t PTCHD4 (Bottlenosed dolphin), MGC148542 (Cow)
	'''
	geneNames = dict()
	geneNameFile =  open(filename, "rU")
	next(geneNameFile)
	for line in geneNameFile:
		line = line.strip().split("\t")
		if len(line) > 1:
			geneNames[line[0][:-2]] = line[1].split(" ")[0] #dont use last 2 chars of identifier
		else:												#Change slice as needed
			geneNames[line[0][:-2]] = line[0]
	return geneNames

def getDndsEntries(filename):
	'''
	retrieves dndsEntry objects from a dnds file
	IN:
	file containing dnds entries
	OUT:
	set of dndsEntries
	'''
	entries = set()

	dndsFile = open(filename, "rU")
	next(dndsFile)

	for line in dndsFile:
		line = line.strip().split("\t")
		m = re.search("(\w+)\.\d", line[0])
		if m:
			line[0] = m.group(1)
		entries.add(dNdSEntry(line[0], line[1], line[2], line[3], line[4]))
	dndsFile.close()
	return entries

def parseGeneNamesInName(filename):
	'''
	retrieves genenames when theyre baked into ortholog 
	identifiers (eg Brandts bat)
	works by splitting the first line up into its component
	parts and checking whether it contains a gene name 
	'''
	geneNames = dict()
	geneNameFile = open(filename, "rU")
	for line in geneNameFile:
		if line[0] == ">":
			line = line[1:].strip().split()
			m = re.search("\w+_([A-Z,a-z,0-9]+)_\d+", line[0])
		
			if m:
				geneName = m.group(1)
				geneNames[line[0]] = geneName
			else:
				geneNames[line[0]] = None

	geneNameFile.close()
	return geneNames

def parseEnsemblGenesFile(filename):
	'''
	Parses a file containing genes retrieved from ensembl
	IN:
	File containing genenames from ensembl
	OUT:
	dict with identifier : gene name
	'''
	geneNames = dict()

	ensFile = open(filename, "rU")
	for line in ensFile:
		line = line.strip().split("\t")
		if len(line) > 1:
			geneNames[line[0]] = line[1]
		else:
			geneNames[line[0]] = line[0]
	ensFile.close()
	return geneNames

def parseGffGenes(filename):
	'''
	Parses GFF file to retrieve gene names
	IN:
	GFF file to parse
	OUT:
	dictionary with identifier : gene name
	'''
	geneNames = dict()

	gff = open(filename, "rU")
	for line in gff:
		m = re.search("ID=(\w*);Source=(\w+);Function=\"(\w*.*?)\";", line)
		if m:
			refID = m.group(1)
			ensID = m.group(2)
			if len(m.group(3)) >= 1:
				geneName = m.group(3)
			else:
				geneName = refID

			geneNames[refID] = geneName
	
		else:
			n = re.search("ID=gene:(\w+);(Name=(\w+))?", line)
			if n:
				refID = n.group(1)
				if n.group(3):
					geneName = n.group(3)
				else:
					geneName = refID

				geneNames[refID] = geneName
	
	gff.close()
	return geneNames

def writeOutputFile(outputDir, filename, geneNames, dndsEntries):
	'''
	Function that puts the gene names and the dndsEntries
	together and writes them to an outputfile
	IN:
	name of outputfile
	dictionary of gene names with identifiers
	set of dndsEntries
	OUT:
	file written on disk
	'''
	filename = "{}/{}GeneNames.tsv".format(outputDir.strip("/"), (filename.split("/")[-1][:-4]))
	
	out = open(filename, "w")
	noName = 0

	for entry in dndsEntries:
		try:
			out.write("{}\t{}".format(geneNames[entry.ref], str(entry)))
			if geneNames[entry.ref] == entry.ref:
				noName += 1
		except KeyError:
			print("WARNING: {} does not have a matching gene name".format(entry.ref))
			out.write("None\t{}".format(str(entry)))

	out.close()
	if noName != 0:
		print("WARNING: no gene name found for {} genes, defaulting to reference ensembl ID".format(noName))

def main():
	'''
	Controller function
	'''
	args = parseCli()
	geneNames = getNamesFromFile(args.geneNames)
	for dndsFile in glob.glob(args.dndsDir+"/*.tsv"):
		print("processing: {}...".format(dndsFile.split("/")[-1]))
		dndsEntries = getDndsEntries(dndsFile)
		writeOutputFile(args.outputDir, dndsFile, geneNames, dndsEntries)

if __name__ == "__main__":
	main()