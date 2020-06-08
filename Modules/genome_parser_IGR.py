#!/usr/bin/env python3
#!/home/niekh/miniconda3/envs/snakemake/bin/python

import os
import sys
import argparse



##This python script is for parsing the input genome (fasta) by the annotated file (gff) into three groups;
# Protein coding genes, RNA genes and non coding regions (Intergenenic + repeating regions).##

__author__= "N.A.H. Huijsmans"
__version__= "genome_parser v1.1.0"
__date__= "16th of March, 2020"

#making arguments
parser = argparse.ArgumentParser(prog="genome_parser", description="",usage="%(prog)s -i <inputfasta> -a <inputgff> [options]",epilog= "Thank you for using genome_parser!")
parser._optionals.title= "Arguments for genome_parser"
parser.add_argument("-v", "--version", help= "prints program version and exits genome_parser",action= "version", version= __version__+" "+__date__+" by "+__author__)
parser.add_argument("-i", metavar= "[input_fasta]", help= "input .fa file for parsing", required= False)
parser.add_argument("-a", metavar= "[input_GFF]", help= "input .gff file for parsing filename", required= False)

parser.add_argument("--keep", metavar= '', help= "keep all files produced by genome_parser", required= False,type=bool, nargs= "?", const= True, default= False)
argument = parser.parse_args()

#when no arguments printing help
if len(sys.argv) == 1:
    parser.print_help(sys.stderr)
    sys.exit(0)

#signing arguments to functions
fasta_file = os.path.abspath(argument.i)
GFF_file = os.path.abspath(argument.a)

keep = argument.keep





#Parsing the GFF file provided by the user
gene_type_dict = {}
gene_type_dict['CDS'] = []
# gene_type_dict['rRNA'] = []
# gene_type_dict['tRNA'] = []
# gene_type_dict['tmRNA'] = []
# gene_type_dict['repeat_region'] = []
total_regions = {}
locus = 0
with open(GFF_file, 'r') as file:
    line = file.readline()
    while line:
        if len(line.split('\t')) == 1:
            if line[0] == '#':
                if line[2] == 's':
                    chromosome = line.split()[1]
                    total_regions[chromosome] = [line.split()[1], int(line.split()[2]), int(line.split()[3]) , '+']
                    line = file.readline()
                    continue
                else:
                    line = file.readline()
                    continue
            else:
                line = file.readline()
                continue
        elif len(line.split('\t')) == 9:
            if line.split('\t')[2] == 'CDS':
                locus += 1
                gene_type_dict['CDS'].append([line.split('\t')[0], int(line.split('\t')[3]), int(line.split('\t')[4]), line.split('\t')[6], locus, line.split('\t')[-1]])
                if (int(line.split('\t')[4]) - int(line.split('\t')[3]) + 1) % 3 != 0:
                    print('Error in CDS length not a multiple of 3 in %s' % line)
                line = file.readline()
                continue
            # if line.split('\t')[2] == 'rRNA':
            #     locus += 1
            #     gene_type_dict['rRNA'].append([line.split('\t')[0], int(line.split('\t')[3]), int(line.split('\t')[4]), line.split('\t')[6], locus, line.split('\t')[-1]])
            #     line = file.readline()
            #     continue
            # if line.split('\t')[2] == 'tRNA':
            #     locus += 1
            #     gene_type_dict['tRNA'].append([line.split('\t')[0], int(line.split('\t')[3]), int(line.split('\t')[4]), line.split('\t')[6], locus, line.split('\t')[-1]])
            #     line = file.readline()
            #     continue
            # if line.split('\t')[2] == 'tmRNA':
            #     locus += 1
            #     gene_type_dict['tmRNA'].append([line.split('\t')[0], int(line.split('\t')[3]), int(line.split('\t')[4]), line.split('\t')[6], locus, line.split('\t')[-1]])
            #     line = file.readline()
            #     continue
            # if line.split('\t')[2] == 'repeat_region':
            #     locus += 1
            #     gene_type_dict['repeat_region'].append([line.split('\t')[0], int(line.split('\t')[3]), int(line.split('\t')[4]), line.split('\t')[6], locus, line.split('\t')[-1]])
            #     line = file.readline()
            #     continue
            else:
                # print("An unknown gene or region is identified but can't be used")
                line = file.readline()
        else:
            print("file does not seem be standard GFF format")
            line = file.readline()
file.close()
# gene_type_dict['Total_regions'] = total_regions

# The function to subtract ranges
from itertools import chain

def range_diff(r1, r2):
    s1, e1 = r1
    s2, e2 = r2
    endpoints = sorted((s1, s2, e1, e2))
    result = []
    if endpoints[0] == s1:
        result.append((endpoints[0], endpoints[1]))
    if endpoints[3] == e1:
        result.append((endpoints[2], endpoints[3]))
    return (result)

def multirange_diff(r1_list, r2_list):
    for r2 in r2_list:
        r1_list = list(chain(*[range_diff(r1, r2) for r1 in r1_list]))
    return (r1_list)

    # Preparing sets of all the annotated regions

### In RefSeq coordinates the last coordiante is part of the gene and it is 1-based

annotated_sets = {}
for chromo in total_regions:
    annotated_sets[chromo] = []
    for gene_types in gene_type_dict:
        for coordiantes in gene_type_dict[gene_types]:
            if coordiantes[0] == chromo:
                annotated_sets[chromo].append(tuple([coordiantes[1] - 1, coordiantes[2] + 1]))
#print(annotated_sets)
# total length of the bacterial genome
total_length_chromo = {}
for chromo in total_regions:
    total_length_chromo[chromo] = [tuple([total_regions[chromo][1], total_regions[chromo][2]])]
#print(total_length_chromo)
# Subtracting total length of bacterial genome from annotated sets to get intergenic region
IGR = {}
for chromo in total_regions:
    IGR[chromo] = multirange_diff(total_length_chromo[chromo], annotated_sets[chromo])
# Discarding very short IGR with length less than 10 and standardize the format
IGR_uniform = []
iter_igr = 1
for chromo in IGR:
    for coordinates in IGR[chromo]:
        if coordinates[1] - coordinates[0] >= 0:
            locus += 1
            IGR_uniform.append([chromo, coordinates[0], coordinates[1], '+', locus, 'ID=IGR_%s' % iter_igr])
            iter_igr += 1
gene_type_dict['IGR'] = IGR_uniform

#making list of regions and genes
list= []
for key in gene_type_dict:
    if key == 'IGR':
        for coor in gene_type_dict["IGR"]:
            list.append([key, coor[1], coor[2], coor[3], coor[5].split('=')[1]])
    else:
        for coor in gene_type_dict["CDS"]:
            list.append([key, coor[1], coor[2], coor[3], 'locus_%s' % coor[4]])

def sortSecond(val):
    return val[1]

list.sort(key=sortSecond)
number = 1
Pathie = 'list.tsv'
with open(Pathie, 'w+') as txt:
    for coordinate in list:
        txt.write('%s\t%s\t%s\t%s\t%s\t%s\n' % (number, coordinate[0], coordinate[1], coordinate[2], coordinate[3], coordinate[4]))
        number += 1
txt.close()

# Retrieving FASTA sequence of all the elements
from Bio import SeqIO

genome_seq = {}

with open(fasta_file, "r") as handle:
    for records in SeqIO.parse(handle, "fasta"):
        genome_seq[records.id] = records.seq

# We can define three categories with different substitution matrices:
# 1) Protein coding genes (CDS)
# 2) RNA genes (tRNA,rRNA,ncRNA,'tmRNA',transcript,misc_feature')
# 3) noncoding sequences (Pseudogene,IGR and repeat regions) #enough for today, save eyes
total_types = {'Protein_coding_gene': ['CDS'],
               'Noncoding_region': ['IGR']}

# Retrieving gene_ids,gene_types,locus_number and sequence of each gene
from Bio.Seq import Seq

total_types_seq = {}
crispr  = 0
alfsimID = 1
for genes_types in total_types:
    total_types_seq[genes_types] = {}
    for type_categories in total_types[genes_types]:
        total_types_seq[genes_types][type_categories] = []
        for coordinates in gene_type_dict[type_categories]:
            if coordinates[-1].split('=')[0] == 'note':
                crispr += 1
                gene_id_cur = str('repeat_region_%s' % crispr)
            else:
                gene_id_cur = [item.split('=')[1].strip() for item in coordinates[-1].split(';') if item.split('=')[0] == 'ID'][0]
            chromo, st, en, strand = coordinates[0], coordinates[1], coordinates[2], coordinates[3]
            gene_type_cur = genes_types
            gene_locus_cur = coordinates[4]
            gene_seq_cur = str(genome_seq[coordinates[0]][coordinates[1] - 1:coordinates[2]])
            if genes_types == 'Protein_coding_gene':
                if coordinates[3] == '+':
                    gene_prot_cur = str(genome_seq[coordinates[0]][coordinates[1] - 1:coordinates[2]].translate())
                if coordinates[3] == '-':
                    gene_prot_cur = str(genome_seq[coordinates[0]][coordinates[1] - 1:coordinates[2]].reverse_complement().translate())
                    gene_seq_cur = str(Seq(gene_seq_cur).reverse_complement())
                total_types_seq[genes_types][type_categories].append([chromo, st, en, strand, gene_id_cur,gene_type_cur, str(alfsimID), gene_seq_cur,gene_prot_cur, coordinates[-1]])
                alfsimID += 1
            else:
                if coordinates[3] == '-':
                    gene_seq_cur = str(Seq(gene_seq_cur).reverse_complement())
                total_types_seq[genes_types][type_categories].append([chromo, st, en, strand, gene_id_cur,gene_type_cur, gene_locus_cur, gene_seq_cur,coordinates[-1]])

#Writing sequences in the fasta file
desired_types = "Protein_coding_gene,Noncoding_region"
for types in desired_types.split(','):
    if types not in total_types:
        print('The desired type: %s not defined!' % types)
    else:
        if types == 'Protein_coding_gene':
            TotalCoordinate = []
            for category in total_types_seq[types]:
                for coords in total_types_seq[types][category]:
                    TotalCoordinate.append(coords)
                PathNow =  'nucleotide_' + types + '.fa'
                PathNow1 = 'protein_' + types + '.fa'
                PathNow2 = types + '.db'
                with open(PathNow, 'w+') as txt:
                    for coordinate in TotalCoordinate:
                        txt.write('>%s, sequence type: %s, locus: %s\n%s\n' % (coordinate[4], coordinate[5], coordinate[6], coordinate[7]))
                txt.close()
                with open(PathNow1, 'w+') as txt:
                    for coordinate in TotalCoordinate:
                        txt.write('>%s, sequence type: %s, locus: %s\n%s\n' % (coordinate[4], coordinate[5], coordinate[6], coordinate[8]))
                txt.close()
                with open(PathNow2, 'w+') as txt:
                    for coordinate in TotalCoordinate:
                        txt.write("<E><ID>" + coordinate[4] +  "</ID><TP>" + coordinate[5] + "</TP><LC>" + coordinate[6]
                                  + "</LC><SEQ>" + coordinate[8].replace("*", "") + "</SEQ><DNA>" + coordinate[7] + "</DNA></E>\n")
                txt.close()

        # if types == 'RNA_gene':
        #     path = os.path.join(output, types)
        #     os.mkdir(path)
        #     TotalCoordinate = []
        #     for category in total_types_seq[types]:
        #         for coords in total_types_seq[types][category]:
        #             TotalCoordinate.append(coords)
        #         PathNow = output + '/' + types + '/' + types + '.fa'
        #         with open(PathNow, 'w') as txt:
        #             for coordinate in TotalCoordinate:
        #                 txt.write('>%s, sequence type: %s, locus: %s\n%s\n' % (coordinate[4], coordinate[5], coordinate[6], coordinate[7]))
        #         txt.close()

        if types == 'Noncoding_region':
            TotalCoordinate = []
            for category in total_types_seq[types]:
                for coords in total_types_seq[types][category]:
                    TotalCoordinate.append(coords)
                PathNow = types + '.fa'
                with open(PathNow, 'w+') as txt:
                    for coordinate in TotalCoordinate:
                        txt.write('>%s, sequence type: %s, locus: %s\n%s\n' % (coordinate[4], coordinate[5], coordinate[6], coordinate[7]))
                txt.close()

            ##keep for writing info detail in excel, json and tsv file##
if keep is True: #Then Writing sequences info detail into an excel file#
#     seq_info= output + '/' + 'output_genome_parser.xlsx'
#     workbook=xlsxwriter.Workbook(seq_info)
#     worksheet=workbook.add_worksheet('Seq_count')
#     worksheet.write(0,0,'Sequence Category')
#     worksheet.write(0,1,'Sequence Subcategory')
#     worksheet.write(0,2,'Count')
#     iter=1
#     for item in total_types:
#         worksheet.write(iter,0,item)
#         total_cur=0
#         for seq_type in total_types[item]:
#             total_cur+=len(gene_type_dict[seq_type])
#         worksheet.write(iter,2,total_cur)
#         iter+=1
#         for item2 in total_types[item]:
#             worksheet.write(iter,1,item2)
#             worksheet.write(iter,2,len(gene_type_dict[item2]))
#             iter+=1
#
#     worksheet = workbook.add_worksheet('Protein_coding_genes')
#     iter = 1
#     for item in total_types['Protein_coding_gene']:
#         for item2 in gene_type_dict[item]:
#             gene_id = gene_id_cur = [item3.split('=')[1].strip() for item3 in item2[-1].split(';') if item3.split('=')[0] == 'ID'][0]
#             gene_name = gene_id_cur = [item3.split('=')[1].strip() for item3 in item2[-1].split(';') if item3.split('=')[0] == 'Name']
#             gene_seq_cur = str(genome_seq[item2[0]][item2[1] - 1:item2[2]])
#             if item2[3] == '+':
#                 gene_prot_cur = str(genome_seq[item2[0]][item2[1] - 1:item2[2]].translate())
#             if item2[3] == '-':
#                gene_prot_cur = str(genome_seq[item2[0]][item2[1] - 1:item2[2]].reverse_complement().translate())
#             worksheet.write(iter, 0, gene_name)
#             worksheet.write(iter, 1, gene_id)
#             worksheet.write(iter, 2, item2[0])
#             worksheet.write(iter, 3, item2[1])
#             worksheet.write(iter, 4, item2[2])
#             worksheet.write(iter, 5, item2[3])
#             worksheet.write(iter, 6, item2[4])
#             worksheet.write(iter, 7, gene_seq_cur)
#             worksheet.write(iter, 8, gene_prot_cur)
#             iter += 1
#
#     worksheet = workbook.add_worksheet('RNA_gene')
#     iter = 1
#     for item in total_types['RNA_gene']:
#         for item2 in gene_type_dict[item]:
#             gene_id = gene_id_cur = [item3.split('=')[1].strip() for item3 in item2[-1].split(';') if item3.split('=')[0] == 'ID'][0]
#             gene_name = gene_id_cur = [item3.split('=')[1].strip() for item3 in item2[-1].split(';') if item3.split('=')[0] == 'gene']
#             gene_seq_cur = str(genome_seq[item2[0]][item2[1] - 1:item2[2]])
#             if len(gene_name) > 0:
#                 worksheet.write(iter, 0, gene_name[0])
#             worksheet.write(iter, 1, gene_id)
#             worksheet.write(iter, 2, item2[0])
#             worksheet.write(iter, 3, item2[1])
#             worksheet.write(iter, 4, item2[2])
#             worksheet.write(iter, 5, item2[3])
#             worksheet.write(iter, 6, item2[4])
#             worksheet.write(iter, 7, gene_seq_cur)
#             iter += 1
#
#     worksheet=workbook.add_worksheet('Non_coding_region')
#     iter=1
#     for item in total_types['Non_coding_region']:
#         for item2 in gene_type_dict[item]:
#             gene_id=[item3.split('=')[1].strip() for item3 in item2[-1].split(';') if item3.split('=')[0]=='ID'][0]
#             gene_name=[item3.split('=')[1].strip() for item3 in item2[-1].split(';') if item3.split('=')[0]=='gene']
#             gene_seq_cur=str(genome_seq[item2[0]][item2[1]-1:item2[2]])
#             if len(gene_name) >0:
#                 worksheet.write(iter,0,gene_name[0])
#             worksheet.write(iter,1,gene_id)
#             worksheet.write(iter,2,item2[0])
#             worksheet.write(iter,3,item2[1])
#             worksheet.write(iter,4,item2[2])
#             worksheet.write(iter,5,item2[3])
#             worksheet.write(iter,6,item2[4])
#             worksheet.write(iter,7,gene_seq_cur)
#             iter+=1
#     workbook.close()
# #writing the coordiante data into a JSON file for later use
#     JSON = output + '/' + 'output_genome_parser.json'
#     with open(JSON, 'w+') as f:
#         json.dump(total_types_seq, f)
# #Creat a map to transform bacterial chromosome names to regular 1-22
# #human chromosome names to be used by GCTA for phenotype simulation
#     MapNum=1
#     ChromMap = output + '/' + 'ChromMap.tsv'
#     with open(ChromMap, 'w+') as f:
#         for chrom in total_regions:
#             f.write('%s\t%s\n'%(chrom,MapNum))
#             MapNum+=1
#             if MapNum >22:
#                 print('More than 22 chromosomes counted in the dataset,this will cause error in GCTA genotype simulation')

    sys.exit(0)

else:
    sys.exit(0)