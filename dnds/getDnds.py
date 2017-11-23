#! /usr/local/bin/python3
'''
Title:	dndsPipeline.py
Author:	Sharlie van der Heide
		(Edited and optimized by David van de Veerdonk)
Pupose:	Calculate dNdS ratio's for ortholog
		genes found by orthoMCL
Date:	10/2017
'''

#TODO turn back no gap no stopcodon? or not?

import subprocess
import argparse
import shutil
import time
import re
import os
from Bio import SeqIO
from Bio.Seq import Seq
from Bio import AlignIO
from Bio.Phylo.PAML import codeml
from Bio.Alphabet import Alphabet
from Bio.SeqRecord import SeqRecord
from Bio.Align import MultipleSeqAlignment
from subprocess import check_call, DEVNULL, STDOUT

def parseCli():
	'''
	parses the commandline arguments needed by this program 
	and provides rudimentary help functions.
	'''

	parser = argparse.ArgumentParser()
	parser.add_argument("--orthologs", "-o", help="ortholgs file", required=True)
	parser.add_argument("--referenceIndex", "-i", help="index (int) of the refence organism in orthologs (0 based)", required=True)
	parser.add_argument("--cds1", "-c1", help="CDS of the reference organism", required=True)
	parser.add_argument("--cds2", "-c2", help="CDS of the second organism", required=True)
	parser.add_argument("--outfile", "-of", help="name to use for final outputfile.", required=True)

	return parser.parse_args()

def parseOrthologs(orthoFile, refIndex):
	'''	
	parses the file containing ortholog gene IDs needs the 
	index of the reference organism
	IN: file containing orthology information
	OUT: dict with org1ID : org2ID
	'''

	orthologs = dict()

	if int(refIndex) == 0:
		orthologIndex = 1
	else:
		orthologIndex = 0
	for line in open(orthoFile):
		if "|" in line:

			line = line.split()
			#if float(line[2]) >= 1:#remove to stop filtering orthomcl
			m = re.search("transcript:(\w+\d+)\.?\d?", line[orthologIndex])
			n = re.search("transcript:(\w+\d+)\.?\d?", line[int(refIndex)])
			if m:
				line[orthologIndex] = "|"+m.group(1)
			if n:
				line[int(refIndex)] = "|"+n.group(1)
			
			orthologs[line[int(refIndex)].split("|")[1]] = line[int(orthologIndex)].split("|")[1]			
		
		else:
			line = line.split()
			if len(line) <= 2:
				orthologs[line[int(refIndex)]] = line[int(orthologIndex)]			

	return orthologs


def getRecords(orthologPairs, cds1, cds2):
	'''
	indexes CDS files and find sequences belonging to each gene
	'''

	refCDS = SeqIO.index(cds1, "fasta")
	ortCDS = SeqIO.index(cds2, "fasta")

	refRec = dict()
	ortRec = dict()
	mitos = []# mitochondrial genes

	notFound = 0

	for orthologPair in orthologPairs:

		MTmatchref = re.search("chromosome:.*?:MT:", refCDS[orthologPair].description)
		MTmatchort = re.search("chromosome:.*?:MT:", ortCDS[orthologPairs[orthologPair]].description)

		if MTmatchref:
			mitos.append(orthologPair)
		elif MTmatchort:
			mitos.append(orthologPairs[orthologPair])

		refRec[orthologPair] = refCDS[orthologPair]
		ortRec[orthologPairs[orthologPair]] = ortCDS[orthologPairs[orthologPair]]

	return refRec, ortRec, mitos

def calculatedNdS(orthologPairs, refRec, ortRec, mitos, outfile, cwd):
	'''
	uses other functions to determine the dnds ratio for 
	every ortholog pair keeps counters to let the user know
	how far along the program is and to keep track of the
	number of orthologs that cant be aligned. Is a 
	controller function.
	'''

	out = open(outfile, "w")
	out.write("ref\tort\tdN/dS\tdN\tdS\n")#header line in outputfile
	fin = time.time() #current time to keep track of speeds
	counter = 0 #counts number of orthologs processed
	noAlignmentAvailable = 0
	ambiguousCodons = 0


	for pair in orthologPairs:#loop through all ortholog pairs
		counter += 1
		SeqIO.write(refRec[pair], cwd+refRec[pair].id+".fa", "fasta")
		SeqIO.write(ortRec[orthologPairs[pair]], cwd+ortRec[orthologPairs[pair]].id+".fa", "fasta")

		alignment = getLastzAlignment(pair+".fa", orthologPairs[pair]+".fa", cwd)#use lastz to align sequences
		maf = mafToPhylip(alignment, refRec, ortRec, cwd, mitos)#convert lastz output (.maf) to phylip file format

		if maf != None:#check whether alignment exists
			phylipFile, mitochondrial = maf
			writeCtl(phylipFile[:-4]+".ctl",phylipFile, phylipFile[:-4], mitochondrial, cwd)#create control file for codeml
			dndsFile = runCodeml(cwd, phylipFile[:-4]+".ctl", phylipFile[:-4], outfile.split("/")[-1][:-4])#run codeml to find dnds

			dnds, dn, ds = parseDnDsFile(dndsFile)#parse codeml output
			if dnds != "-nan" and dn != "-nan" and ds != "-nan":
				out.write("{}\t{}\t{}\t{}\t{}\n".format(pair, orthologPairs[pair], dnds, dn, ds))
			
			else:
				ambiguousCodons += 1

			os.remove(dndsFile)#cleanup
			os.remove(cwd + phylipFile)
			os.remove(cwd + phylipFile[:-4]+".ctl")
		else:
			noAlignmentAvailable += 1 #if alignment fails counter is upped and printed at the end

		if counter % 10 == 0:#update progress bar for every 10 orthologs processed
			done = "#"*round((counter/len(orthologPairs)*100)/2)#progress bar calculation
			notDone = "-"*round(50 - (counter/len(orthologPairs)*100)/2)
			print("\rruntime({}m) progress: [{}{}]"\
				.format(str(round((time.time()-fin)/60, 2)), done, notDone), end="\r")

	print("No alignment available for {} orthologs.\n\tThis can be due to extra stop codons or gaps in alignment".format(str(noAlignmentAvailable)))
	if ambiguousCodons > 0:
		print("No dNdS available for {} orthologs due to 'synonymous' amino acids with different sequences".format(ambiguousCodons))

def getLastzAlignment(refFasta, ortFasta, cwd):
	'''
	Runs lastZ to get an alignment of 2 provided fasta files
	to keep the directory ordered fastas are cleaned up as 
	theyre no longer needed
	'''

	ref = os.path.splitext(os.path.basename(refFasta))[0]#getting filenames for outputname
	ort = os.path.splitext(os.path.basename(ortFasta))[0]
	alignmentFilename = "{}_{}.maf".format(ref, ort) # output filename
	# --nogapped doesnt allow gaps in the aligned sequeces
	check_call("lastz {} {} --output='{}' --format=maf --ambiguous=iupac --nogapped".format(refFasta, ortFasta, alignmentFilename), shell=True, cwd=cwd)#runs lastz with fastas and outfilename
	
	os.remove(cwd + refFasta)#cleanup of fastas
	os.remove(cwd + ortFasta)

	return alignmentFilename

def mafToPhylip(maf, refRec, ortRec, cwd, mitos):
	'''
	Converts a .maf file to a .phy file. phy files are compatible with PAML
	'''

	stopCodons = ["TAG", "TGA", "TAA"]
	MTstopCodons = ["AGA", "AGG", "TAG", "TAA"]
	mitochondrial = 0
	try:
		alignment = next(AlignIO.parse(cwd+maf, "maf"))
		refAlignment = alignment[0]
		ortAlignment = alignment[1]
		gapLocations = set()
		stopCodonLocations = set()
		
		ungappedSeqs = findAndRemoveGaps([str(refAlignment.seq), str(ortAlignment.seq)])
		if refAlignment.id in mitos or ortAlignment.id in mitos:
			noStopSeqs, isWithoutExtraStops = removeStopCodons(ungappedSeqs, MTstopCodons)
			mitochondrial = 1
		else:
			noStopSeqs, isWithoutExtraStops = removeStopCodons(ungappedSeqs, stopCodons)

		if isWithoutExtraStops:
		
			ortRec = SeqRecord(Seq(noStopSeqs[1]), id = ortAlignment.id[:5])#make 2 seqrecords so they can be realigned
			refRec = SeqRecord(Seq(noStopSeqs[0]), id = refAlignment.id[:5])
			
			msalignment = MultipleSeqAlignment([refRec, ortRec])#realign sequences

			os.remove(cwd + maf)
			
			phylipFilename = maf[:-4]+".phy" #construct filename from .maf filename
			output_handle = open(cwd + phylipFilename, "w")
			AlignIO.write(msalignment, output_handle, "phylip-sequential")#Write phylip alignment file
			output_handle.close()
			
			return phylipFilename, mitochondrial
		
		os.remove(cwd + maf) #also remove alignment if extra stops are there

	except StopIteration:
		os.remove(cwd + maf) #ditto
		return None

def findAndRemoveGaps(seqs):
	'''removes gaps from sequence pairs'''

	gapLocations = set()
	ungappedSeqs = list()
	for seq in seqs:
		[gapLocations.add(x.start()) for x in re.finditer("-", seq)]#find all indexes of gaps
	
	for seq in seqs:
		noGapSeq = "".join([char for idx, char in enumerate(seq) if idx not in gapLocations])
		if len(noGapSeq) % 3 == 1:
			noGapSeq = noGapSeq[:-1]
		elif len(noGapSeq) % 3 == 2:
			noGapSeq = noGapSeq[:-2]
		ungappedSeqs.append(noGapSeq)
	
	return ungappedSeqs


def removeStopCodons(ungappedSeqs, stopCodons):
	'''Removes stopcodons from a pair of aligned sequences'''

	stopCodonLocations = set() 
	for seq in ungappedSeqs:
		codons = [str(seq[i:i+3]) for i in range(0, len(str(seq)), 3)]
		[stopCodonLocations.add(i) for i,codon in enumerate(codons) if codon in stopCodons]

	noStopSeq = list()
	for seq in ungappedSeqs:
		if len(stopCodonLocations) > 1:
			isWithoutExtraStops = False
		else:
			isWithoutExtraStops = True
		if len(stopCodonLocations) >= 1:
			codons = [str(seq[i:i+3]) for i in range(0, len(str(seq)), 3)]
			noStopSeq.append(''.join([codon if i not in stopCodonLocations else "" for i, codon in enumerate(codons)]))
		else:
			noStopSeq.append(seq)

	return noStopSeq, isWithoutExtraStops

def writeCtl(filename, phylipFilename, outfilename, mitochondrial, cwd):
	'''
	writes a ctl control file for PAML this has the path
	to the input files, output names and some settings
	'''

	ctlFile = open(cwd + filename, "w")
	ctlFile.write(\
		"seqfile = {}\n"\
		"outfile = {}\n"\
		"noisy = 0\n"\
		"verbose = 0\n"\
		"runmode = -2\n"\
		"seqtype = 1\n"\
		"CodonFreq = 2\n"\
		"clock = 0\n"\
		"aaDist = 0\n"\
		"aaRatefile = 'wag.dat\n"\
		"model = 0\n"\
		"NSsites = 0\n"\
		"icode = {}\n"\
		"Mgene = 0\n"\
		"fix_kappa = 0\n"\
		"kappa = 2\n"\
		"fix_omega = 0\n"\
		"omega = 1\n"\
		"fix_alpha = 1\n"\
		"alpha = 0.\n"\
		"Malpha = 0\n"\
		"ncatG = 8\n"\
		"getSE = 0\n"\
		"RateAncestor = 1\n"\
		"Small_Diff = .5e-6\n"\
		"cleandata = 1\n"\
		"method = 0\n"\
		.format(phylipFilename, outfilename, mitochondrial))#changed cleandata/omega
	ctlFile.close()

def runCodeml(cwd, ctl, filename, wd):
	'''
	Runs codeml to calculate dNdS ratios. uses the ctl file
	generated by writeCtl() to get all neccesary variables.
	'''
	
	check_call("codeml ./{}".format(ctl), shell=True, stdout=DEVNULL, stderr=STDOUT, cwd=cwd)

	return cwd+filename

def parseDnDsFile(dndsFile):
	'''
	opens the codeml output file and retrieves Dn/Ds ratio, Dn and Ds
	from it. these are returned in a list.
	regex matches line in codeml output where dnds is printed
	eg: t= 0.9626  S=     0.0  N=    57.0  dN/dS=  0.4343  dN = 0.3209  dS =   0.001
	dnds, dn, ds are grouped. codeml has variable spacing in its output
	'''
	
	dndsData = open(dndsFile, "rU")

	dnds = re.search("dN/dS\s*?=\s*?(\d+\.\d+)\s*?dN\s*?=\s*?(\d+\.\d+)\s*?dS\s*?=\s*?(\d+\.\d+|-nan)", dndsData.readlines()[-1])
	dndsData.close()

	return dnds.group(1), dnds.group(2), dnds.group(3)

def createTmpDir(path):
	'''
	Creates a temporary directory to store files
	that are currently being processed
	'''

	if not os.path.exists(directory):
		os.makedirs(directory)

def main():
	'''
	Controller function. parses commandline and starts calculate dnds
	'''

	args = parseCli()

	tmpDir = "{}/tmp/{}".format(os.getcwd(), (args.outfile).split("/")[-1][:-4])
	if not os.path.exists(tmpDir):
		os.makedirs(tmpDir)

	start = time.time()
	print("parsing orthologs file...")
	orthologPairs = parseOrthologs(args.orthologs, args.referenceIndex)
	print("OK.. time taken: {}s\nRetrieving ortholog records from CDS...".format(round(time.time() - start, 2)))
	refRec, ortRec, mitos = getRecords(orthologPairs, args.cds1, args.cds2)
	print("OK.. time taken: {}s\nCalculating dNdS for {} ortholog genes...\nThis will take aproximately {} minutes.".format(round(time.time() - start, 2), len(orthologPairs), round(len(orthologPairs)/162.4)))
	calculatedNdS(orthologPairs, refRec, ortRec, mitos, args.outfile, tmpDir+"/")
	print("OK... time taken: {}\nCleaning up directories".format(round(time.time()-start, 2)))
	
	if os.path.exists(tmpDir):
		shutil.rmtree(tmpDir)


if __name__ == "__main__":
	main()