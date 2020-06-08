#!/usr/bin/env python3
#!/home/niekh/miniconda3/envs/snakemake/bin/python

import os
import sys
import argparse


##This python script is for building individual dawg files for simulating evolution the input contigs (fasta) are write as a dawg file. ##

__author__= "N.A.H. Huijsmans"
__version__= "dawg_builder v1.1.0"
__date__= "16th of March, 2020"

#making arguments
parser = argparse.ArgumentParser(prog="dawg_builder", description="",usage="%(prog)s -i <inputfasta> -a <inputnwk> [options]",epilog= "Thank you for using dawg_builder!")
parser._optionals.title= "Arguments for dawg_builder"
parser.add_argument("-v", "--version", help= "prints program version and exits dawg_builder",action= "version", version= __version__+" "+__date__+" by "+__author__)
parser.add_argument("-i", metavar= "[input_fasta]", help= "input .fa file for building dawg file", required= False)
parser.add_argument("-a", metavar= "[input_tree]", help= "input .nwk or .tre for building dawg file", required= False)
parser.add_argument("-o", metavar= "[output]", help= "output_dir", required= False, default= "dawg_builder")
argument = parser.parse_args()

#when no arguments printing help
if len(sys.argv) == 1:
    parser.print_help(sys.stderr)
    sys.exit(0)

#signing arguments to functions
fasta_file = os.path.abspath(argument.i)
nwk_file = os.path.abspath(argument.a)
output = os.path.abspath(argument.o)




#nwk file is a one line string used for dawg file
nwk = open(nwk_file, 'r')
nwk_f = nwk.read()

#reading fasta file
total=[]
with open(fasta_file, "r") as file:
    line = file.read()
    lines = line.replace("\n", ",")
    l = lines.split(">")
    for lin in l:
        row = lin.replace(':', ',')
        if len(row) < 1:
            continue
        else:
            total.append([row.split(',')[0], row.split(',')[2], row.split(',')[4], row.split(',')[5]])
file.close()

##writing dawg files with fixed params and models
for coord in total:
    if coord[0].split('_')[0] == 'IGR':
        PathNow = output + '/' + coord[0] + '.dawg'
        with open(PathNow, 'w') as txt:
            txt.write('[Root]\nSeq = %s\n\n[Tree]\nTree = %s\n[Output]\nFile =%s.fas\n\n[Indel]\nModel = US\nRate = 0.00175'
                      '\nParams = 0.320238687,0.140726007,0.134261561,0.085529587,0.055693685,0.050223769,0.032819493,0.026852312,0.032819493,0.019890602,0.018896072,0.020885132,0.012928891,0.007956241,0.008950771,0.008950771,0.006464446,0.006961711,0.004972650,0.003978120\n\n'
                      '[Subst]\nModel = GTR\nFreqs = 0.236,0.249,0.279,0.236\nParams = 0.91770,4.47316,1.10375,0.56499,6.01846,1.00000\nRate.Model = Gamma\nRate.Params = 1.0,0.01\n\n[Sim]\nReps = 1\nseed = 1' % (coord[-1], nwk_f, coord[0]))
        txt.close()
    else:
        PathNow = output + '/locus_' + coord[-2].strip() + '.dawg'
        with open(PathNow, 'w') as txt:
            txt.write('[Root]\nSeq = %s\n\n[Tree]\nTree = %s\n[Output]\nFile =locus_%s.fas\n\n[Indel]\nModel = US\nRate = 0.00175'
                      '\nParams = 0.320238687,0.140726007,0.134261561,0.085529587,0.055693685,0.050223769,0.032819493,0.026852312,0.032819493,0.019890602,0.018896072,0.020885132,0.012928891,0.007956241,0.008950771,0.008950771,0.006464446,0.006961711,0.004972650,0.003978120\n\n'
                      '[Subst]\nModel = GTR\nFreqs = 0.236,0.249,0.279,0.236\nParams = 0.91770,4.47316,1.10375,0.56499,6.01846,1.00000\nRate.Model = Gamma\nRate.Params = 1.0,0.01\n\n[Sim]\nReps = 1\nseed = 1' % (coord[-1], nwk_f, coord[-2].strip()))
        txt.close()