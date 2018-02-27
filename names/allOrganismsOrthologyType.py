#!/usr/bin/python3

import argparse

def parseCli():
	'''
	Parses commandline arguments
	'''
	parser = argparse.ArgumentParser()
	parser.add_argument("--orthology_type_dir", "-o", help="directory woth orthology type files")
	parser.add_argument("--genenames_dir", "-g", help="directory with gene name files")
	parser.add_argument("--outfile", "-of", help="output file")

	return parser.parse_args()

def parseOrthologyFile(filename):
	'''
	parses orthology file
	'''
	pass

def parseGeneNamesFile(filename):
	'''
	parses gene name file
	'''
	pass

def main():
	'''
	Controller function
	'''
	args = parseCli()

if __name__ == "__main__":
	main()