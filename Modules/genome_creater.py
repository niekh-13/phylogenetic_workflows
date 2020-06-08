#!/usr/bin/env python3
#!/home/niekh/miniconda3/envs/snakemake/bin/python

import os
import sys
import argparse
from Bio import SeqIO
from Bio.Seq import Seq
import random


#total time


##This python script is for assembling input fasta files in directorye by the ALF file (fasta) into genomes##

__author__= "N.A.H. Huijsmans"
__version__= "genome_creater v1.1.0"
__date__= "16th of March, 2020"

#making arguments
parser = argparse.ArgumentParser(prog="genome_creater", description="",usage="%(prog)s -i <inputdir> -d <inputdir> -l <inputfile> -q <int> [options]",epilog= "Thank you for using genome_creater!")
parser._optionals.title= "Arguments for genome_creater"
parser.add_argument("-v", "--version", help= "prints program version and exits genome_parser",action= "version", version= __version__+" "+__date__+" by "+__author__)
parser.add_argument("-i", metavar= "[input_alf]", help= "input file with .fa file for creating genome", required= False)
parser.add_argument("-d", metavar= "[input_dir_dawg]", help= "input dir with .fas file for creating genome", required= False)
parser.add_argument("-l", metavar= "[input_list_file]", help= "input .tsv file from output of parsing", required= False)
parser.add_argument("-o", metavar= "[output]", help= "output file", required= False, default= "sample.fasta")
parser.add_argument("--keep", metavar= '', help= "keep all files produced by genome_creater", required= False,type=bool, nargs= "?", const= True, default= False)
argument = parser.parse_args()

#when no arguments printing help
if len(sys.argv) == 1:
    parser.print_help(sys.stderr)
    sys.exit(0)

#signing arguments to functions
alf_file = os.path.abspath(argument.i)
dawg_dir = os.path.abspath(argument.d)
list_file = os.path.abspath(argument.l)
output = os.path.abspath(argument.o)
keep = argument.keep

#def for sorting
def sortSecond(val):
    return int(val[0])
#open list file for creating the list where the file names for IGR and RNA genes can be found
the_list = []
with open(list_file, "r") as l_file:
    list = l_file.readline()
    while list:
        the_list.append([list.split('\t')[0], list.split('\t')[1], int(list.split('\t')[2]), int(list.split('\t')[3]), list.split('\t')[4], list.split('\t')[5].strip('\n')])
        list = l_file.readline()

#making lists for random IGR
alf_order = []
dawg_list = []
for line in the_list:
    if line[1] == "IGR":
        sum = (int(line[3]) - int(line[2]))
        if sum < 100:
            dawg_list.append(line[-1])
    if line[1] == "CDS":
        alf_order.append(line[-1].split('_')[-1])



#opening alf_file
alf_list = alf_file.split('/')[-1]
list = []
a = random.randint(0, (len(dawg_list)-1))
correct = alf_list.split('_')[0]
number = 1
not_list = []
cluster_list = []
normal_list = []
#opening the alf file
with open(alf_file, 'r') as file:
    f = file.read()
    lines = f.split('>')
    for line1 in lines:
        if len(line1.split()) > 1:
            x = line1.split()[-2]
            #if statement for loci that are from strange donor.
            if x not in alf_order:
                not_list.append(line1.split('\n'))
                list_dawg = []
                IGR_name = dawg_dir + '/' + dawg_list[a] + '.fas'
                with open(IGR_name, 'r') as dawg1:
                    f2 = dawg1.read()
                    file3 = f2.split('>')
                    for d in file3:
                        if len(d.strip('\n').split()) > 1:
                            seq0 = ''.join(d.strip('\n').split()[1:])
                            list_dawg.append([d.strip('\n').split()[0], seq0])
                seq1 = [line3[1] for line3 in list_dawg if line3[0] == correct]
                # after some decoding and finding the a random IGR append in dict with number of place in the dict
                list.append([line1.split()[0], number, '+', 'IGR', seq1[0].replace('-', '')])
                list.append([line1.split()[0], number,'+', 'CDS',  line1.split()[-1] + 'TAA'])
                number += 1
            #normal loci
            else:
                y = [location[0] for location in the_list if location[-1].split('_')[-1] == x if location[1] == "CDS"]
                pos = (int(y[0]) - 1)
                IGR = [location[-1] for location in the_list if int(location[0]) == pos if location[1] == 'IGR']
                #if statement for CDS with no IGR in front of them
                if len(IGR) == 0:
                    cluster_list.append(line1.split('\n'))
                    complement = [location[-2] for location in the_list if location[-1].split('_')[-1] == x if location[1] == 'CDS']
                    if complement[0] == '-':
                        list.append([line1.split()[0], number, complement[0], 'CDS', str(Seq(line1.split()[-1] + 'TAA').reverse_complement())])
                    else:
                        list.append([line1.split()[0], number, complement[0], 'CDS', line1.split()[-1] + 'TAA'])
                    number += 1
                # CDS with IGR in front of them
                else:
                    normal_list.append(line1.split('\n'))
                    list_dawg = []
                    IGR_name = dawg_dir + '/' + IGR[0] + '.fas'
                    with open(IGR_name, 'r') as dawg2:
                        f4 = dawg2.read()
                        file4 = f4.split('>')
                        for d4 in file4:
                            if len(d4.strip('\n').split()) > 1:
                                seq3 = ''.join(d4.strip('\n').split()[1:])
                                list_dawg.append([d4.strip('\n').split()[0], seq3])
                    seq2 = [line4[1] for line4 in list_dawg if line4[0] == correct]
                    # after some decoding and finding the right IGR append in dict with number of place in the dict
                    list.append([number, '+', 'IGR', seq2[0].replace('-', '')])
                    complement = [location[-2] for location in the_list if location[-1].split('_')[-1] == x if location[1] == 'CDS']
                    if complement[0] == '-':
                        list.append([line1.split()[0], number, complement[0], 'CDS', str(Seq(line1.split()[-1] + 'TAA').reverse_complement())])
                    else:
                        list.append([line1.split()[0], number, complement[0], 'CDS', line1.split()[-1] + 'TAA'])
                    number += 1

    text0 = (str(alf_list) + ' is written in dictionary. \n')
    sys.stdout.write(text0)

text1 = 'Writting fasta file from dictionary starts now.\n'
sys.stdout.write(text1)


#writing fasta files of genomes in dictionary
TotalCoordinate = []
with open(output, 'w+') as txt:
    txt.write('>%s, complete genome \n' % (alf_file.split('.')[0].split("/")[-1]))
    for category in list:
        txt.write(str(category[-1]))
txt.close()
text3 = (str(alf_file) + ' is written to fasta\n')
sys.stdout.write(text3)
sys.exit(0)




