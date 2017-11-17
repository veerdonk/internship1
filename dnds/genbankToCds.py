#! /usr/local/bin/python3
'''
Title:	genbankToCds.py
Author:	David van de Veerdonk
Pupose:	Retrieve coding dna sequences 
		from a genbank file
'''

import argparse
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq

def parseCli():
	'''
	Parses the commandline arguments
	'''
	parser = argparse.ArgumentParser()
	parser.add_argument("--genbank", "-g", help="name of the genbank file to parse")
	parser.add_argument("--outfile", "-o", help="name of the output file.")

	return parser.parse_args()

def parseGenbank(infilename, outfilename):
	'''
	Parses genbank file and retrieves coding seqs
	'''
	print(infilename)
	print(outfilename)
	genbank = SeqIO.parse(infilename, "genbank")
	print("indexed")
	recs = list()
	cds = 0
	failed  = 0
	for i, seqRecord in enumerate(genbank):

		print("records processed: {}{}".format(i, (i%4)*"."), end = "\r")
		for feature in seqRecord.features:
			if feature.type == "CDS":
				try:
					featureName = feature.qualifiers["protein_id"][0]
					featureSeq = feature.extract(seqRecord.seq)
					
					recs.append(SeqRecord(featureSeq, id = featureName, \
						description = feature.qualifiers["gene"][0]))
					#out.write(">{}\n{}\n".format(featureName, featureSeq))
					cds += 1
				except:
					failed += 1

	print("sequences collected: {}\nNo ID for {} sequences".format(cds, failed))			
	SeqIO.write(recs, outfilename, "fasta")
def main():
	'''
	Controller function
	'''
	args = parseCli()
	parseGenbank(args.genbank, args.outfile)

if __name__ == "__main__":
	main()