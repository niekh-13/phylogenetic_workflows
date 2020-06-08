#!/usr/bin/env python3
#!/home/niekh/miniconda3/envs/snakemake/bin/python

import os
import sys
import argparse

##This python script is for assembling input fasta files in directorye by the ALF file (fasta) into genomes##

__author__= "N.A.H. Huijsmans"
__version__= "genome_assessor v1.1.0"
__date__= "3th of April, 2020"

#making arguments
parser = argparse.ArgumentParser(prog="genome_assessor", description="",usage="%(prog)s -i <input_dir> -q <int> o- <file_name> [options]",epilog= "Thank you for using genome_assessor!")
parser._optionals.title= "Arguments for genome_assessor"
parser.add_argument("-v", "--version", help= "prints program version and exits genome_parser",action= "version", version= __version__+" "+__date__+" by "+__author__)
parser.add_argument("-g", metavar= "[input_dir]", help= "input dir with .fasta file of genomes", required= False)
parser.add_argument("-q", metavar= "[quantity_files]", help= "how many fasta files contain the dir?", required= False, type=int)
parser.add_argument("-m", metavar= "[input_dir]", help= "input dir with .paf file of minimap2", required= False)
parser.add_argument("-p", metavar= "[input_dir]", help= "input dir with .tsv file of prokka", required= False)
parser.add_argument("-a", metavar= "[input_dir]", help= "input dir with .fa file of alf", required= False)
parser.add_argument("-o", metavar= "[output]", help= "output_file", required= False, default= "genome_assessed.tsv")
parser.add_argument("-t", metavar= "[threads]", help= "amount of threads [max available up to 256 threads]", required= False, type=int, default= 16)
parser.add_argument("--keep", metavar= '', help= "keep all files produced by genome_assessor", required= False,type=bool, nargs= "?", const= True, default= False)
argument = parser.parse_args()

#when no arguments printing help
if len(sys.argv) == 1:
    parser.print_help(sys.stderr)
    sys.exit(0)

#signing arguments to functions
genome_dir = os.path.abspath(argument.g)
alf = os.path.abspath(argument.a)
map = os.path.abspath(argument.m)
prokka = os.path.abspath(argument.p)
output = argument.o
threads = str(argument.t)
max = argument.q

#calculating genome size and CG content#
#first open file#
query_list = []
genome_dict = {}
for number in range(max):
    paf_list = []
    list = []
    list1 = []
    list3 = []
    if number < 9:
        genome_file = genome_dir + '/' + 'SE00' + str(number+1) + '.fasta'
        prokka_tsv = prokka + "/prokka-SE00" + str(number+1) + "/prokka-SE00" + str(number+1) + ".tsv"
        prokka_txt = prokka + "/prokka-SE00" + str(number + 1) + "/prokka-SE00" + str(number + 1) + ".txt"
        paf_file =  map + '/SE00' + str(number + 1) + '.paf'
        alf_file = alf + '/SE00' + str(number + 1) + '_dna.fa'
        name = "SE00" + str(number+1)

    else:
        genome_file = genome_dir + '/' + 'SE0' + str(number+1) + '.fasta'
        prokka_tsv = prokka + "/prokka-SE0" + str(number+1) + "/prokka-SE0" + str(number+1) + ".tsv"
        prokka_txt = prokka + "/prokka-SE00" + str(number + 1) + "/prokka-SE00" + str(number + 1) + ".txt"
        paf_file =  map + '/SE0' + str(number + 1) + '.paf'
        alf_file = alf + '/SE0' + str(number + 1) + '_dna.fa'
        name = "SE0" + str(number+1)
    genome_list = genome_file.split('/')[-1]
    genome_dict[genome_list] = []
    fi2 = open(genome_file, "r")
    f2 = fi2.read()
    file2 = f2.split('\n')
    G = 0
    C = 0
    total = len(file2[1]) - 1
    for nucleotide in file2[1]:
        if nucleotide == 'G':
            G += 1
        if nucleotide == 'C':
            C += 1
    GC_content = ((G+C) / total * 100)
    file_tsv = open(prokka_tsv, "r")
    tsv = file_tsv.read()
    f_tsv = tsv.split('\n')
    CDS_number = 0
    CDS = 0
    for line in f_tsv:
        if line.startswith('locus_tag'):
            continue
        elif len(line.split('\t')) < 3:
            continue
        else:
            if line.split('\t')[1] == "CDS":
                CDS_number += 1
                CDS += int(line.split('\t')[2])
            else:
                continue
    avg_length = (CDS / CDS_number)
    percentage_cds = ((CDS / total) * 100)
    density = (CDS_number / (total / 1000))
    # file_txt = open(prokka_txt, "r")
    # txt = file_txt.read()
    # f_txt = txt.split('\n')
    # for line1 in f_txt:
    #     if line1.startswith("rRNA"):
    #         # print(line1)
    #         rRNA = line1.split(':')[-1]
    #     if line1.startswith("tRNA"):
    #         # print(line1)
    #         tRNA = line1.split(':')[-1]
    #     if line1.startswith("tmRNA"):
    #         # print(line1)
    #         tmRNA = line1.split(':')[-1]
    #     if line1.startswith("repeat_region"):
    #         # print(line1)
    #         repeat_region = line1.split(':')[-1]
    with open(paf_file,"r") as file_paf:
        paf_f = file_paf.read()
        paf_lines = paf_f.split('\n')
        for lin in paf_lines:
            paf_list.append(lin.split('\t'))
        for lines in paf_list:
            if len(lines) > 1:
                list.append(lines[0].split('_')[0])
    with open(alf_file, "r") as file_alf:
        alf_lines = file_alf.read()
        alf_pre_list = alf_lines.split('>')
        for pre_lines in alf_pre_list:
            if len(pre_lines) > 1:
                list1.append(pre_lines.split()[0].replace(',','').split('_')[0])
    list2 = []
    for i in list:
        if i not in list2:
            list2.append(i)
    res = set(list1) - set(list2)
    genome_dict[genome_list].append([total, float(round(GC_content, 2)), float(round(avg_length, 2)),
                                     float(round(percentage_cds, 2)), float(round(density, 3)), len(list1), len(list2)])
#write one file
with open(output, 'w') as txt:
    txt.write('genome\tgenome_size\tgc-content\tmean_cds_size\tcds(%)\t#cds/1kbp\tCDS_alf(#)\tCDS_genome(#)\n')
    for key in genome_dict:
        txt.write(key.split('.')[0] + '\t')
        for coor in genome_dict[key]:
            txt.write(str(coor[0]) + '\t' + str(coor[1]) + '\t' + str(coor[2]) + '\t' + str(coor[3])+ '\t' +
                          str(coor[4])+ '\t' + str(coor[5]) + '\t' + str(coor[6])+ '\n')
txt.close()
sys.exit(0)