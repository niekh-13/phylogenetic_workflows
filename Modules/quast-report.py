#!/usr/bin/env python3
#!/home/niekh/miniconda3/envs/snakemake/bin/python

import sys
import argparse



__author__= "N.A.H. Huijsmans"
__version__= "quast-report v1.1.0"
__date__= "3th of April, 2020"

#making arguments
parser = argparse.ArgumentParser(prog="quast-report", description="",usage="%(prog)s -q <int> o- <file_name> [options]",epilog= "Thank you for using quast_reporter!")
parser._optionals.title= "Arguments for quast-report"
parser.add_argument("-v", "--version", help= "prints program version and exit",action= "version", version= __version__+" "+__date__+" by "+__author__)
parser.add_argument("-q", metavar= "[quantity_genomes]", help= "how many genomes?", required= False, type=int)
parser.add_argument("-o", metavar= "[output]", help= "output_file", required= False, default= "quast-report")
argument = parser.parse_args()

#when no arguments printing help
if len(sys.argv) == 1:
    parser.print_help(sys.stderr)
    sys.exit(0)

#signing arguments to functions
max = argument.q
output = argument.o

#output file
if output == True:
    print("Error(1): Output file already exists, exiting")
    sys.exit(1)

dict = {}

for number in range(max):
    if number < 9:
        spades = "quast/quast_spades-SE00" + str(number + 1) + "/report.tsv"
        skesa = "quast/quast_skesa-SE00" + str(number + 1) + "/report.tsv"
        shovill_v = "quast/quast_shovill_velvet-SE00" + str(number + 1) + "/report.tsv"
        # shovill_s = "quast/quast_shovill_spades-SE00" + str(number + 1) + "/report.tsv"
        prokka = "prokka/prokka-SE00" + str(number + 1) + "/prokka-SE00" + str(number + 1) + ".tsv"
        list = "SE00" + str(number + 1)
    else:
        spades = "quast/quast_spades-SE0" + str(number + 1) + "/report.tsv"
        skesa = "quast/quast_skesa-SE0" + str(number + 1) + "/report.tsv"
        shovill_v = "quast/quast_shovill_velvet-SE0" + str(number + 1) + "/report.tsv"
        # shovill_s = "quast/quast_shovill_spades-SE0" + str(number + 1) + "/report.tsv"
        prokka = "prokka/prokka-SE0" + str(number + 1) + "/prokka-SE0" + str(number + 1) + ".tsv"
        list = "SE0" + str(number + 1)
    dict[list]= []
    CDS_number = 0
    with open(prokka, 'r') as prokka_gff:
        prokka_line = prokka_gff.readline()
        while prokka_line:
            if prokka_line.startswith('locus_tag'):
                prokka_line = prokka_gff.readline()
                continue
            elif len(prokka_line.split('\t')) < 3:
                prokka_line = prokka_gff.readline()
                continue
            else:
                CDS_number +=1
                prokka_line = prokka_gff.readline()


    with open(spades, 'r') as spades_report:
        spades_line = spades_report.readline()
        while spades_line:
            if spades_line.startswith("# genomic features"):
                if len(spades_line.split()) == 7:
                    complete = int(spades_line.split()[3])
                    part = int(spades_line.split()[5])
                    complete_features = round(((complete/CDS_number)*100),2)
                    total_features = round((((complete+part)/CDS_number)*100),2)
                    dict[list].append(str(complete_features))
                    dict[list].append(str(total_features))
                    spades_line = spades_report.readline()
                    continue
                else:
                    spades_line = spades_report.readline()
                    continue
            if spades_line.startswith("# contigs"):
                if len(spades_line.split()) == 3:
                    dict[list].append(spades_line.split('\t')[-1].strip("\n"))
                    spades_line = spades_report.readline()
                    continue
                else:
                    spades_line = spades_report.readline()
                    continue
            elif spades_line.startswith("NGA50"):
                dict[list].append(spades_line.split('\t')[-1].strip("\n"))
                spades_line = spades_report.readline()
                continue
            elif spades_line.startswith("# misassemblies"):
                dict[list].append(spades_line.split('\t')[-1].strip("\n"))
                spades_line = spades_report.readline()
                continue
            elif spades_line.startswith("Genome"):
                dict[list].append(spades_line.split('\t')[-1].strip("\n"))
                spades_line = spades_report.readline()
                continue
            else:
                spades_line = spades_report.readline()
    text0 = (str(list) + '-spades is written in dictionary, moving to next quast report. \n')
    sys.stdout.write(text0)
    with open(skesa, 'r') as skesa_report:
        skesa_line = skesa_report.readline()
        while skesa_line:
            if skesa_line.startswith("# genomic features"):
                if len(skesa_line.split()) == 7:
                    complete = int(skesa_line.split()[3])
                    part = int(skesa_line.split()[5])
                    complete_features = round(((complete/CDS_number)*100),2)
                    total_features = round((((complete+part)/CDS_number)*100),2)
                    dict[list].append(str(complete_features))
                    dict[list].append(str(total_features))
                    skesa_line = skesa_report.readline()
                    continue
                else:
                    skesa_line = skesa_report.readline()
                    continue
            if skesa_line.startswith("# contigs"):
                if len(skesa_line.split()) == 3:
                    dict[list].append(skesa_line.split('\t')[-1].strip("\n"))
                    skesa_line = skesa_report.readline()
                    continue
                else:
                    skesa_line = skesa_report.readline()
                    continue
            elif skesa_line.startswith("NGA50"):
                dict[list].append(skesa_line.split('\t')[-1].strip("\n"))
                skesa_line = skesa_report.readline()
                continue
            elif skesa_line.startswith("# misassemblies"):
                dict[list].append(skesa_line.split('\t')[-1].strip("\n"))
                skesa_line = skesa_report.readline()
                continue
            elif skesa_line.startswith("Genome"):
                dict[list].append(skesa_line.split('\t')[-1].strip("\n"))
                skesa_line = skesa_report.readline()
                continue
            else:
                skesa_line = skesa_report.readline()
    text1 = (str(list) + '-skesa is written in dictionary, moving to next quast report. \n')
    sys.stdout.write(text1)
    with open(shovill_v, 'r') as shovill_v_report:
        shovill_v_line = shovill_v_report.readline()
        while shovill_v_line:
            if shovill_v_line.startswith("# genomic features"):
                if len(shovill_v_line.split()) == 7:
                    complete = int(shovill_v_line.split()[3])
                    part = int(shovill_v_line.split()[5])
                    complete_features = round(((complete/CDS_number)*100),2)
                    total_features = round((((complete+part)/CDS_number)*100),2)
                    dict[list].append(str(complete_features))
                    dict[list].append(str(total_features))
                    shovill_v_line = shovill_v_report.readline()
                    continue
                else:
                    shovill_v_line = shovill_v_report.readline()
                    continue
            if shovill_v_line.startswith("# contigs"):
                if len(shovill_v_line.split()) == 3:
                    dict[list].append(shovill_v_line.split('\t')[-1].strip("\n"))
                    shovill_v_line = shovill_v_report.readline()
                    continue
                else:
                    shovill_v_line = shovill_v_report.readline()
                    continue
            elif shovill_v_line.startswith("NGA50"):
                dict[list].append(shovill_v_line.split('\t')[-1].strip("\n"))
                shovill_v_line = shovill_v_report.readline()
                continue
            elif shovill_v_line.startswith("# misassemblies"):
                dict[list].append(shovill_v_line.split('\t')[-1].strip("\n"))
                shovill_v_line = shovill_v_report.readline()
                continue
            elif shovill_v_line.startswith("Genome"):
                dict[list].append(shovill_v_line.split('\t')[-1].strip("\n"))
                shovill_v_line = shovill_v_report.readline()
                continue
            else:
                shovill_v_line = shovill_v_report.readline()
    # text2 = (str(list) + '-shovill_velvet is written in dictionary, moving to next quast report. \n')
    # sys.stdout.write(text2)
    # with open(shovill_s, 'r') as shovill_s_report:
    #     shovill_s_line = shovill_s_report.readline()
    #     while shovill_s_line:
    #         if shovill_s_line.startswith("# genomic features"):
    #             if len(shovill_s_line.split()) == 7:
    #                 complete = int(shovill_s_line.split()[3])
    #                 part = int(shovill_s_line.split()[5])
    #                 complete_features = round(((complete/CDS_number)*100),2)
    #                 total_features = round((((complete+part)/CDS_number)*100),2)
    #                 dict[list].append(str(complete_features))
    #                 dict[list].append(str(total_features))
    #                 shovill_s_line = shovill_s_report.readline()
    #                 continue
    #             else:
    #                 shovill_s_line = shovill_s_report.readline()
    #                 continue
    #         if shovill_s_line.startswith("# contigs"):
    #             if len(shovill_s_line.split()) == 3:
    #                 dict[list].append(shovill_s_line.split('\t')[-1].strip("\n"))
    #                 shovill_s_line = shovill_s_report.readline()
    #                 continue
    #             else:
    #                 shovill_s_line = shovill_s_report.readline()
    #                 continue
    #         elif shovill_s_line.startswith("NGA50"):
    #             dict[list].append(shovill_s_line.split('\t')[-1].strip("\n"))
    #             shovill_s_line = shovill_s_report.readline()
    #             continue
    #         elif shovill_s_line.startswith("# misassemblies"):
    #             dict[list].append(shovill_s_line.split('\t')[-1].strip("\n"))
    #             shovill_s_line = shovill_s_report.readline()
    #             continue
    #         elif shovill_s_line.startswith("Genome"):
    #             dict[list].append(shovill_s_line.split('\t')[-1].strip("\n"))
    #             shovill_s_line = shovill_s_report.readline()
    #             continue
    #         else:
    #             shovill_s_line = shovill_s_report.readline()
#     text2 = (str(list) + '-shovill_spades is written in dictionary, moving to next quast report. \n')
#     sys.stdout.write(text2)
# text = ("all quast reports are written in dictionary, moving to writing the reports. \n")
# sys.stdout.write(text)

with open(output + "_spades.tsv", 'w') as txt:
    txt.write("genome\t#contigs\tmisassemblies\tgenome_fraction\tcomplete_features\ttotal_features\tNGA50\n")
    for key in dict:
        txt.write(key + '\t' + dict[key][0] + '\t' + dict[key][1] + '\t' + dict[key][2] + '\t' + dict[key][3]+ '\t' + dict[key][4]+ '\t' + dict[key][5] + '\n')
text3 = ("spades quast report is written, moving to writing next report. \n")
sys.stdout.write(text3)

with open(output + "_skesa.tsv", 'w') as txt1:
    txt1.write("genome\t#contigs\tmisassemblies\tgenome_fraction\tcomplete_features\ttotal_features\tNGA50\n")
    for key in dict:
        txt1.write(key + '\t' + dict[key][6] + '\t' + dict[key][7] + '\t' + dict[key][8] + '\t' + dict[key][9] +'\t' + dict[key][10]+ '\t' + dict[key][11] + '\n')
text4 = ("skesa quast report is written, moving to writing next report. \n")
sys.stdout.write(text4)

with open(output + "_shovill_velvet.tsv", 'w') as txt2:
    txt2.write("genome\t#contigs\tmisassemblies\tgenome_fraction\tcomplete_features\ttotal_features\tNGA50\n")
    for key in dict:
        txt2.write(key + '\t' + dict[key][12] + '\t' + dict[key][13] + '\t' + dict[key][14] + '\t' + dict[key][15] +'\t' + dict[key][16]+ '\t' + dict[key][17] + '\n')
text5 = ("shovill_velvet quast report is written, moving to writing next report. \n")
sys.stdout.write(text5)

# with open(output + "_shovill_spades.tsv", 'w') as txt2:
#     txt2.write("genome\t#contigs\tmisassemblies\tgenome_fraction\tcomplete_features\ttotal_features\tNGA50\n")
#     for key in dict:
#         txt2.write(key + '\t' + dict[key][18] + '\t' + dict[key][19] + '\t' + dict[key][20] + '\t' + dict[key][21] + '\t' + dict[key][22]+ '\t' + dict[key][23] + '\n')
# text7 = ("shovill_spades quast report is written, moving to writing next report. \n")
# sys.stdout.write(text7)

with open(output + ".tsv", 'w') as txt3:
    txt3.write("genome\tcontigs\tmisassemblies\tgenome_fraction\tcomplete_features\ttotal_features\tNGA50\tassembler\n")
    for key in dict:
        txt3.write(key + "\t" + dict[key][0] + '\t' + dict[key][1] + '\t' + dict[key][2] + '\t' + dict[key][3] + '\t'  + dict[key][4] + '\t' + dict[key][5] + "\t" +  "spades" + "\n" +
                   key + "\t" + dict[key][6] + '\t' + dict[key][7] + '\t' + dict[key][8] + '\t' + dict[key][9] + '\t' + dict[key][10] + '\t' + dict[key][11] + '\t' + "skesa" + '\n' +
                   key + "\t" + dict[key][12] +'\t' + dict[key][13] + '\t' + dict[key][14] + '\t' + dict[key][15] + "\t" + dict[key][16] + '\t' + dict[key][17] + "\t" "shovill_velvet" + '\n')
                   # key + "\t" + dict[key][18] +'\t' + dict[key][19] + '\t' + dict[key][20] + '\t' + dict[key][21] + "\t" + dict[key][22] + '\t' + dict[key][23] + "\t"  "shovill_spades" + '\n')

text6 = ("all quast reports are written, exiting. \n")
sys.stdout.write(text6)
sys.exit(0)