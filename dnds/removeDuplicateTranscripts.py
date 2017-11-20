#!/usr/bin/python3
'''
Title:		removeDuplicateTranscripts.py
Author: 	David van de Veerdonk
Purpose:	remove shorter transcripts of the same
			gene and writing the longest to a file
Date:		10/2017
'''
import argparse
import re
from Bio import SeqIO

def parseCli():
	'''
	Parses commandline arguments
	'''
	parser = argparse.ArgumentParser()
	parser.add_argument("--infile", "-i", help="Fasta to remove extra transcripts from")
	parser.add_argument("--outfile", "-o", help="name of the output file")		
	return parser.parse_args()

def removeDupes(infile, outfile):
	'''
	Opens and reads a fasta file using a regular expression it finds the location of the gene
	identifier and saves it to a dictionary. Then the longest transcript for every gene is 
	selected and written to a ne file with the gene identifier as record ID.
	IN:
	Fastafile to be filtered
	OUT:
	fasta with only longest transcripts
	and geneID as record id
	'''
	
	dupes = 0
	genes = dict()
	records = list(SeqIO.parse(infile, "fasta"))#parse input fasta

	print("removing shorter transcripts from a pool of {} total transcripts".format(len(records)))

	for rec in records: # |(\w+_\d+\.?\d?) <- append to regex for files in case regex fails to pick up on ids
		m = re.search(".*gene:(\w*\d*)\.?\d?.*|.*\[protein_id=(\w*\.?\d?)\].*|geneid=(\w+)", rec.description)#matches several different gene/protein names in fastas
		if m:
			if m.group(1) != None:#checking what id this file uses
				gene = m.group(1)
			elif m.group(2) != None:
				gene = m.group(2)
			elif m.group(3) != None:
				gene = m.group(3)
			elif m.group(4) != None:
				gene = m.group(4)

		if gene in genes:#adding gene to a dictionary
			dupes += 1
			genes[gene].append(rec)
		else:
			genes[gene] = [rec]



	outfile = open(outfile, "w")
	
	for gene in genes:
		longest = 0 
		for rec in genes[gene]:
			if len(rec.seq) > longest:#selecting the longest transcript
				longest = len(rec.seq)
				longestRec = rec
		longestRec.id = gene #makes the gene the identifier in the fasta (comment out for original id)

		SeqIO.write(longestRec, outfile, "fasta")#write transcripts to file
	outfile.close()		

	print(str(dupes) + " shorter transcripts removed from fasta.")
	print("out of a total of " + str(len(records)) + " records")

def main():
	'''
	Controller function
	'''
	args = parseCli()
	removeDupes(args.infile, args.outfile)
	

if __name__ == "__main__":
	main()