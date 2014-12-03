import sys
import getopt
import string
import os
import re
import glob

###########################################################################
### Credit
###########################################################################
def credit():
	print "pdb2fasta script, Brian Chen, Honig lab, 2009."

###########################################################################
### Usage
###########################################################################
def usage():
	credit();
	print("========================================================="   );
	print(""                                                          );
	print("python pdb2fasta.py"        );
	print("   Running with no input returns this usage info."  );
	print(""                                                          );
	print("python pdb2fasta.py (-switch) [pdbfile(s)] [output.fasta]"       );
	print("   Primary Usage.  Parses the PDB file(s) and outputs their"  );
	print("   amino acids to fasta files."  );
	print("   SWITCHES: -s or -t"  );
	print("   -s: Separate.  PDB files with several chains are entered"  );
	print("      into the FASTA output as seperate entities."  );
	print("   -t: Together.  PDB files with several chains are entered"  );
	print("      into the FASTA output as a single sequence."  );
	print("   [pdbfile(s)]: a path, including standard unix wildcards");
	print("      that lists a (potentially empty) set of PDB files." );
	print("      No tilde expansion, but ../ ./ * and ? are handled");
	print("   [output.fasta]: the FASTA file that contains the sequences");
	print("");
	print("WARNING: CURRENTLY ONLY -t IS IMPLEMENTED RIGHT NOW");
	print("");
	print("IMPORTANT ISSUES AND FEATURES TO NOTE:");
	print("   -alternate coordinates are (correctly) ignored");
	print("   -ONLY canonical amino acid names are recognized");
	print("   -HETATM cards are not processed");
	print("   -It is assumed that your PDB file names end with a .pdb");
	print("   -It is assumed that ATOMS in the same chain are sequential.");
	print("      i.e. you dont have atoms from chain A, then B, then A again.");
	print("========================================================="   );



###########################################################################
###########################################################################
### PDB to FASTA
###########################################################################
###########################################################################
def aminoAcidLookup(line):
	aaCode = line[17:20];
	
	if aaCode == "ALA":	return "A";
	if aaCode == "ARG":	return "R";
	if aaCode == "ASN":	return "N";
	if aaCode == "ASP":	return "D";
	if aaCode == "CYS":	return "C";
	if aaCode == "GLU":	return "E";
	if aaCode == "GLN":	return "Q";
	if aaCode == "GLY":	return "G";
	if aaCode == "HIS":	return "H";
	if aaCode == "ILE":	return "I";
	if aaCode == "LEU":	return "L";
	if aaCode == "LYS":	return "K";
	if aaCode == "MET":	return "M";
	if aaCode == "PHE":	return "F";
	if aaCode == "PRO":	return "P";
	if aaCode == "SER":	return "S";
	if aaCode == "THR":	return "T";
	if aaCode == "TRP":	return "W";
	if aaCode == "TYR":	return "Y";
	if aaCode == "VAL":	return "V";
	if aaCode == "MSE":	return "M";
	print("ERROR: AMINO ACID CODE [" + aaCode + "] not identified.  Skipping");
	


###########################################################################
###########################################################################
### fixPDB - reads a PDB file, eliminates garbage.  Takes NO chain or A chains only.
#NOTE - this does not fix the HETATM methionines.
###########################################################################
###########################################################################
##COLUMNS      DATA TYPE        FIELD      DEFINITION
##------------------------------------------------------
## 1 -  6      Record name      "ATOM    "
## 7 - 11      Integer          serial     Atom serial number.
##13 - 16      Atom             name       Atom name.
##17           Character        altLoc     Alternate location indicator.
##18 - 20      Residue name     resName    Residue name.
##22           Character        chainID    Chain identifier.
##23 - 26      Integer          resSeq     Residue sequence number.
##27           AChar            iCode      Code for insertion of residues.
##31 - 38      Real(8.3)        x          Orthogonal coordinates for X in 
##                                         Angstroms
##39 - 46      Real(8.3)        y          Orthogonal coordinates for Y in 
##                                         Angstroms
##47 - 54      Real(8.3)        z          Orthogonal coordinates for Z in 
##                                         Angstroms
##55 - 60      Real(6.2)        occupancy  Occupancy.
##61 - 66      Real(6.2)        tempFactor Temperature factor.
##77 - 78      LString(2)       element    Element symbol, right-justified.
##79 - 80      LString(2)       charge     Charge on the atom.
##

###########################################################################
###########################################################################
### PDB to FASTA
###########################################################################
###########################################################################
def pdb2fasta(separate, pdbPath, outputFileName):
	
#	credit();
#	print("========================================================="   );

	for pfile in pdbPath:

		#get the code:
		pdbCode = pfile.split(".")[0].split("/")[-1];

		processFile(separate, pfile, pdbCode, outputFileName);


###########################################################################
###########################################################################
### PDB to FASTA
###########################################################################
###########################################################################
def processFile(separate, pdbFile, pdbCode, outputFileName):

	print("PDBFILE: " + pdbFile);
	print("OUTPUTFILE: " + outputFileName );

	outputFile = open(outputFileName, 'a');

	currentChain = "@"
	numProcessed = 0;

	###open the file
	pdbFile = open(pdbFile, 'rt');
	line = pdbFile.readline()
	
	if not separate:
		outputFile.write(">" + pdbCode.upper() +"\n");
		
	while line:
		if not ( line.startswith("ATOM") and line[12:16].split()[0] == "CA" and (line[16] == " " or line[16] == "A" ) ):
			line = pdbFile.readline();
			continue;

		aaChar = aminoAcidLookup(line)
		outputFile.write(aaChar);
		numProcessed = numProcessed + 1;
		if numProcessed == 100:
			outputFile.write("\n");
			numProcessed = 0;

		line = pdbFile.readline();

	if numProcessed != 0:
		outputFile.write("\n");
	
	pdbFile.close();
	outputFile.close();



###########################################################################
###########################################################################
### Main Method
###########################################################################
###########################################################################
def main():
	
#	testDebug();
	
	#print sys.argv[1]
	#print len(sys.argv)
	#print "Hello, World!"

	################################################ Usage Statement
	if len(sys.argv) == 1:
		usage()
		raise SystemExit, 5

	################################################ Help Statement
	else:
		separate = False;

		if sys.argv[1] == "-s":
			separate = True;
		elif sys.argv[1] == "-t":
			separate = False;
		else:
			print("ERROR: incorrect switch provided.  Quitting.\n");
			raise SystemExit, 5

		pdb2fasta(separate, sys.argv[2:-1], sys.argv[-1]);

		raise SystemExit, 5
	################################################ 


if __name__ == "__main__":
	main()
	
