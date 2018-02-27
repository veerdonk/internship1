#!/usr/bin/python3

import argparse
import glob

from retrieveStops import parseCds

def parseCli():
	'''
	Parses commandline arguments
	'''
	parser = argparse.ArgumentParser()
	parser.add_argument("--groupsfile", "-g", help="orthomcl groupsfile")
	parser.add_argument("--organismCodes", "-c", help="threeletter codes of all the organisms in the groupsfile (separated by ',')")
	parser.add_argument("--CDSdir", "-d", help="directory with CDS files")
	parser.add_argument("--outfile", "-o", help="name for output file")

	return parser.parse_args()

def parseGroups(filename, codes):
	'''
	parses groups file
	'''
	groupsfile = open(filename, "rU")
	groupsToPrint = list()

	for line in groupsfile:
		line = line.strip().split("\t")
		codesOnLine = list()
		for gene in line:
			codesOnLine.append(gene.split("|")[0])
		# print(set(codesOnLine))
		if set(codes) == set(codesOnLine):
			groupsToPrint.append(line)

	
	groupsfile.close()

	return groupsToPrint

def main():
	'''
	Controller function
	'''
	args = parseCli()
	orgcodes = args.organismCodes.strip().split(",")
	geneIDs = parseGroups(args.groupsfile, orgcodes)
	cdsFiles = glob.glob(args.CDSdir +"/*.*")
	stops = dict()
	conserved = {"TAA" : 0, "TGA" : 0, "TAG" : 0}
	allStops = dict()
	for file in cdsFiles:
		print(file)
		stops.update(parseCds(file))

	legalStops = ["TAA", "TAG", "TGA"]
	idout = open("/home/dylan/data/orthologs/young/idToConservedStop.csv", "w")
	for group in geneIDs:
		stopsInGroup = set()
		for gene in group:
			if gene.startswith("min") or gene.startswith("tut"):
				gene = gene.split(".")[0]

			gene = gene.split("|")[1]
			
			if gene.startswith("bmy"):
				gene = gene[:-2]

			if stops[gene] in allStops:
				allStops[stops[gene]] += 1
			else:
				allStops[stops[gene]] = 1

			if stops[gene] in legalStops:
				stopsInGroup.add(stops[gene])

		if len(stopsInGroup) == 1:
			stop = stopsInGroup.pop()
			conserved[stop] += 1
			for geneid in group:
				if geneid.startswith("csa"):
					idout.write("{}\t{}\n".format(geneid.split("|")[1], stop))
			
	
	
	sumConserved = (conserved["TAA"] + conserved["TGA"] + conserved["TAG"])

	out = open(args.outfile, "w")
	out.write(",TGA,TAA,TAG\n")
	out.write("number,{},{},{}\n".format(conserved["TGA"], conserved["TAA"], conserved["TAG"]))
	out.write("fraction,{},{},{}\n".format(conserved["TGA"]/sumConserved, conserved["TAA"]/sumConserved, conserved["TAG"]/sumConserved))
	out.write("fractionOfAll,{},{},{}\n".format(conserved["TGA"]/allStops["TGA"], conserved["TAA"]/allStops["TAA"], conserved["TAG"]/allStops["TAG"]))

	print("Number of groups with conserved stopcodon:\nTGA:\t{}\nTAA:\t{}\nTAG:\t{}\n".format(conserved["TGA"], conserved["TAA"], conserved["TAG"]))

	print("Conserved fraction:\nTGA:\t{}\nTAA:\t{}\nTAG:\t{}\n".format(conserved["TGA"]/sumConserved, conserved["TAA"]/sumConserved, conserved["TAG"]/sumConserved))

	print("Fraction of non-conserved + conserved:\nTGA:\t{}\nTAA:\t{}\nTAG:\t{}".format(conserved["TGA"]/allStops["TGA"], conserved["TAA"]/allStops["TAA"], conserved["TAG"]/allStops["TAG"]))

	out.close()

if __name__ == "__main__":
	main()