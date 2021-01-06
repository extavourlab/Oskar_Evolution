#!/usr/bin/env python

'''
Author: Savandara Besse

Created date: 04-07-2017
Modified date : 04-19-2017

Description: This script allows run hmmsearch for
all protein files obtained for a specific analysis (TSA - GCF - GCA)

'''

import os, sys, re
import progressbar
from subprocess import Popen
import subprocess
from optparse import OptionParser

parser = OptionParser()
parser.add_option("-a", "--analysisType", dest="analysisType", default="None",
                  help="[Required] Type of analysis, TSA - GCA - GCF")
parser.add_option("-p", "--proteinPath", dest="proteinPath", default="None",
                  help="[Required] Location of the protein folder")
parser.add_option("-o", "--oskModel", dest="oskModel", default="None",
                  help="[Required] Location of osk hmm")
parser.add_option("-l", "--lotusModel", dest="lotusModel", default="None",
                  help="[Required] Location of lotus hmm")
parser.add_option("-r", "--outputPath", dest="outputPath", default="None",
                  help="[Required] Location of oskar_tracker folder result")

(options, args) = parser.parse_args()


proteinPath = options.proteinPath
oskModel = options.oskModel
lotusModel = options.lotusModel
analysis_type = options.analysisType
output_path = options.outputPath

protein_folder = os.listdir(proteinPath)

if not os.path.isdir(output_path):
    os.mkdir(output_path)

bar = progressbar.ProgressBar()
hmmsearch = 'hmmsearch'
for proteome in bar(protein_folder) :
    if 'TSA' in analysis_type:
        proteome_name = proteome.split('_')[1][:6]
    else:
        proteome_name = '{}_{}'.format(proteome.split('_')[0], proteome.split('_')[1])

    input_path = os.path.join(proteinPath, proteome)
    if os.path.isfile(input_path) :
        os.system('{} --cpu 8 --tblout {}/{}_osk_search.txt {} {}'.format(hmmsearch, output_path, proteome_name, oskModel, input_path))
        os.system('{} --cpu 8 --tblout {}/{}_lotus_search.txt {} {}'.format(hmmsearch, output_path, proteome_name, lotusModel, input_path))
