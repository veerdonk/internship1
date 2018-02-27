#!/usr/bin/python3
import argparse
import glob
import re

def parseCli():
	'''
	Parses commandline arguments
	'''
	parser = argparse.ArgumentParser()
	parser.add_argument("--inDir", "-i", help="file containing orthologs to be analyzed")
	parser.add_argument("--orgCode", "-c", help="3 letter code of reference organism")
	parser.add_argument("--outDir", "-o", help="output directory")

	return parser.parse_args()


def parseOrtho(filename, orgCode):
	'''
	Parses a .txt containing orthology information. seperates the data
	into two dicts. the dicts are mirrors of each other. one containing
	reference ID - ortholog ID mappings and the other orthologID - 
	reference mappings.

	IN: filename, index of reference (1/0)
	OUT: 2 dicts containing ID mappings
	'''

	orthologsRef = dict()
	orthologsOrt = dict()

	orthoFile = open(filename, "rU")
	first = True
	for line in orthoFile:
		if first:
			if line.startswith(orgCode):
				index = 0
			else:
				# print(line.strip().split("\t")[0])
				m = re.search("\w*({})\w*".format(orgCode), line.strip().split("\t")[0], flags=re.IGNORECASE)
				if m:
					index = 0
				else:
					index = 1
			first = False

		line = line.strip().split("\t")
		if index == 0:
			if line[index] in orthologsRef:
				orthologsRef[line[index]].append(line[1])
			else:
				orthologsRef[line[index]] = [line[1]]
			if line[1] in orthologsOrt:
				orthologsOrt[line[1]].append(line[index])
			else:
				orthologsOrt[line[1]] = [line[index]]
		else:
			if line[index] in orthologsRef:
				orthologsRef[line[index]].append(line[0])
			else:
				orthologsRef[line[index]] = [line[0]]

			if line[0] in orthologsOrt:
				orthologsOrt[line[0]].append(line[index])
			else:
				orthologsOrt[line[0]] = [line[index]]

	return orthologsRef, orthologsOrt


def assignOrthologyType(orthologsRef, orthologsOrt):
	'''
	assigns type of orthology to orthologs
	ie 1:1 1:many many:1 many:many
	
	'''
	counter = 0
	one2one = dict()
	one2many = dict()
	many2one = dict()
	many2many = dict()

	for gene in orthologsRef:
		# if both the reference and the ortholog only have 1 id attached its 1:1		
		if len(orthologsRef[gene]) == 1 and len(orthologsOrt[orthologsRef[gene][0]]) == 1:
			one2one[gene] = orthologsRef[gene][0]
			

		# if the reference has only 1 id attached and the ortholog more than one
		# and all ids attached to the ortholog only have 1 id attached many:1
		elif len(orthologsRef[gene]) == 1 and len(orthologsOrt[orthologsRef[gene][0]]) > 1:
			many2oneflag = True
			
			for refGene in orthologsOrt[orthologsRef[gene][0]]:
				if len(orthologsRef[refGene]) != 1 and refGene != orthologsRef[gene][0]:
					many2oneflag = False

			if many2oneflag:
				many2one[gene] = orthologsRef[gene][0]
				
		elif len(orthologsRef[gene]) > 1:
			one2manyflag = True
			for ortGene in orthologsRef[gene]:
				if len(orthologsOrt[ortGene]) == 1:
					
					if orthologsOrt[ortGene][0] != gene:
						one2manyflag = False
				else:
					one2manyflag = False		

			if one2manyflag:
				one2many[gene] = orthologsRef[gene]
				
			else:
				if not gene in many2many:
					genesToCheckRef = [gene]
					genesToCheckOrt = list()
					for gene in genesToCheckRef:
						for i in orthologsRef[gene]:
							if i not in genesToCheckOrt:
								genesToCheckOrt.append(i)
						for geneOrt in genesToCheckOrt:
							for j in orthologsOrt[geneOrt]:
								if j not in genesToCheckRef:
									genesToCheckRef.append(j)

					for refGene in genesToCheckRef:
						many2many[refGene] = orthologsRef[refGene]

	# for gene in one2one:
	# 	if len(orthologsRef[gene]) != 1 or len(orthologsOrt[orthologsRef[gene][0]]) != 1:
	# 		print("{}:{} - is not one2one".format(gene, orthologsRef[gene][0]))

	# for gene in one2many:
	# 	if not len(orthologsRef[gene]) > 1:
	# 		for ort in orthologsRef[gene]:
	# 			if len(orthologsOrt[ort]) != 1 or orthologsOrt[ort][0] != gene:
	# 				print("{} - is not one2many".format(gene))



	return one2one, one2many, many2one, many2many


def writeOutFile(filename, one2one, one2many, many2one, many2many):
	"""
	writes output file with types of orthology
	also keeps counters of the number of genes
	in each category and outputs this to the 
	commandline

	IN: filename, dictionaries containing ID mappings
	for each type of orthology
	OUT: file containing ID \t orthology \t ID
	"""

	out = open(filename, "w")

	one2oneCounter = 0
	one2manyCounter = 0
	many2oneCounter = 0
	many2manyCounter = 0

	for gene in one2one:
		out.write("{}\tone2one\t{}\n".format(gene, one2one[gene]))
		one2oneCounter += 1

	for gene in one2many:
		for ortGene in one2many[gene]:
			out.write("{}\tone2many\t{}\n".format(gene, ortGene))
			one2manyCounter += 1

	for gene in many2one:
		out.write("{}\tmany2one\t{}\n".format(gene, many2one[gene]))
		many2oneCounter += 1

	for gene in many2many:
		for ort in many2many[gene]:
			out.write("{}\tmany2many\t{}\n".format(gene, ort))
			many2manyCounter += 1

	print("one2one:\t{}\nmany2one\t{}\none2many:\t{}\nmany2many\t{}\ntotal:\t\t{}".format(
		one2oneCounter,
		many2oneCounter,
		one2manyCounter,
		many2manyCounter,
		one2oneCounter+one2manyCounter+many2oneCounter+many2manyCounter))

	out.close()


def main():
	'''
	Main controller function
	walks through given directory and
	uses functions to process files
	in directory
	'''

	args = parseCli()
	# orthologsRef, orthologsOrt = parseOrtho(args.orthologs, args.refindex)
	# one2one, one2many, many2one, many2many = assignOrthologyType(orthologsRef, orthologsOrt)
	# writeOutFile(args.filename, one2one, one2many, many2one, many2many)

	for filename in glob.glob(args.inDir+"/*.txt"):
		print(filename)
		outfile = args.outDir + filename.split("/")[-1][:-4] + "Type.txt"
		print("Parsing ortholog file...")
		orthologsRef, orthologsOrt = parseOrtho(filename, args.orgCode)
		print("Done.\nStarting with orthology type assignment...")
		one2one, one2many, many2one, many2many = assignOrthologyType(orthologsRef, orthologsOrt)
		print("Done.\nWriting output. stats:")
		writeOutFile(outfile, one2one, one2many, many2one, many2many)
		


if __name__ == "__main__":
	main()