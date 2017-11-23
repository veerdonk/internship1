#!/usr/bin/python3
import argparse
import glob
import sys
import re

from dNdSEntry import dNdSEntry

def parseCli():
	'''
	Parses commandline arguments
	'''
	parser = argparse.ArgumentParser(formatter_class=argparse.RawDescriptionHelpFormatter, description="Has two run modes. basic usage:\n\
		Provide a directory containing your files with dnds values (--dndsDir)\n\
		Provide a file (can be .txt/.gff/.cds/.tsv) containing gene names of your reference organism (--geneNames)\n\
		Provide a directory to write the output to (--outputDir)\n\
		\nAdvanced usage:\n\
		If the normal runmode skips a large number of genes specify --both_names \n\
		use these arguments instead of the normal ones:\n\
		The file you want to get the gene names for (--dndsFile)\n\
		2 files containing gene names, one for reference, one for ortholog (--namesRef, --namesOrt)\n\
		Name of the file to write the output to (--outFile)")
	parser.add_argument("--dndsDir", "-d", help="directory containing output of getDnds.py that you want to assign reference gene names")
	parser.add_argument("--geneNames", "-g", help="file containing matching identifiers with dndsFile and gene names")
	parser.add_argument("--outputDir", "-o", help="directory where you want the output saved")

	parser.add_argument("--both_names", help="use two files containing gene names, one with reference names\
	 the other with gene names of it's pair. NOTE: requires the --dndsFile, --namesRef, --namesOrt and --outFile (-df, -nr, -no, -of) params\
	 and does nothing with --dndsDir, --geneNames, --outputDir when set", action="store_true")
	parser.add_argument("--dndsFile", "-df", help="file containing IDs to assign geneNames to, only use with --both_names")
	parser.add_argument("--namesRef", "-nr", help="file containing gene names matching those of the reference organism")
	parser.add_argument("--namesOrt", "-no", help="file containing gene names matching those of the orholog")
	parser.add_argument("--outFile", "-of", help="the file to write --both-names mode output to")


	return parser.parse_args(), parser

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
			else:
				o = re.search("Name=(\w+\.?\d?);.*gene=([\w\-]*);", line)
				if o:
					refID = o.group(1)
					if o.group(2):
						geneName = o.group(2)
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
		print("no gene name found for {} genes, defaulting to reference ensembl ID".format(noName))

def doubleModeOutput(outFile, namesRef, namesOrt, entries):
	'''
	writes the output for --both_names mode
	'''

	noName = 0
	moreNames = 0
	out = open(outFile, "w")
	for entry in entries:
		if namesRef[entry.ref] == entry.ref:
			if namesOrt[entry.ort] == entry.ort: # second check to see if the name is available for the other organism
				out.write("{}\t{}".format(namesRef[entry.ref], str(entry)))
				noName += 1
			else:
				out.write("{}\t{}".format(namesOrt[entry.ort], str(entry)))
				moreNames += 1
		else:
			out.write("{}\t{}".format(namesRef[entry.ref], str(entry)))
	out.close()
	print("no gene name found for {} genes, defaulting to reference ensembl ID\n\
		using --both_names found {} more gene names".format(noName, moreNames))

def main():
	'''
	Controller function
	Checks whether the right command line arguments were given
	select the right runmode/functions
	calls parser to display help if something is amiss
	'''

	args, parser = parseCli()
	if args.both_names == False:
		geneNames = getNamesFromFile(args.geneNames)
		for dndsFile in glob.glob(args.dndsDir+"/*.tsv"):
			print("processing: {}...".format(dndsFile.split("/")[-1]))
			dndsEntries = getDndsEntries(dndsFile)
			writeOutputFile(args.outputDir, dndsFile, geneNames, dndsEntries)

	elif args.namesRef is not None or args.namesOrt is not None or args.outFile is not None or args.outFile is not None:
		namesRef = getNamesFromFile(args.namesRef)
		namesOrt = getNamesFromFile(args.namesOrt)
		entries = getDndsEntries(args.dndsFile)

		doubleModeOutput(args.outFile, namesRef, namesOrt, entries)

	else:
		parser.print_help()
		sys.exit(1)

if __name__ == "__main__":
	main()