#!/usr/bin/env python3

####### Author
# Pauline Fourgoux p.f.m.fourgoux@qmul.ac.uk ; pauline.fourgoux@gmail.com
####### Packages call
import pandas as pd
import sys
import os
####### Funtions
def usage():
  print("""
        computing_RNAseq_AF.py takes the ID of the sample to open a vcf file and calculates the real AF for RNA-seq mutation calling.
        
        Example of command line:
        =======================
        
        python computing_RNAseq_AF.py --in-vcf myin.vcf --out-vcf myout.vcf
        
        obligatory:
        ==========
        
        --in-vcf  -> VCF
        --out-vcf -> VCF
        """)

####### Main
# Get arguments
infile=""
outfile=""
#---------------------------------------------------
try:
  infile = sys.argv[sys.argv.index("--in-vcf")+1]
print("Infile: ", infile)
except:
  usage()
print("ERROR : Please enter the input vcf", infile)
sys.exit()
try:
  outfile = sys.argv[sys.argv.index("--out-vcf")+1]
print("Outfile: ", outfile)
except:
  usage()
print("ERROR : please enter the input vcf file", outfile)
sys.exit()


vcf = pd.read_table(infile,skiprows=56)
rows = vcf.shape[0]

for r in list(range(rows)):
  #print(vcf.loc[[r]])
  if vcf.at[r,'FORMAT'] == "GT:AD:DP:GQ:PL":
  # Retrieving AD and DP to calculate real AF
  #print(vcf.loc[[r]])
  sample = vcf.at[r,'sample']
sample = sample.split(':')
DP = float(sample[2])
#print(DP)
if DP>0:
  AD = sample[1]
#print(AD)
# The comma ',' separates the AD of the different allele.
# We want to take the last one that is the one of the mutation allele.
if "," in AD:
  AD = AD.split(',')
AD = AD[-1]
AD = float(AD)
AF = round(AD/DP,4)
# Replacing AF in the INFO column
INFO = vcf.at[r,'INFO']
INFO = INFO.split(';')
INFO[1] = "AF=" + str(AF)
#print(AF)
#print(INFO[1])
INFO = ';'.join(INFO)
vcf.at[r,'INFO'] = INFO
#print(vcf.loc[[r]])

print("Writing output file:",outfile)
vcf.to_csv(outfile, '\t', header=True, index=None)
print("Output written")
