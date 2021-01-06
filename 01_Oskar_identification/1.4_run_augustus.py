#!/usr/bin/env python

'''
Author: Savandara Besse

Created date: 04-26-2017
Modified date: 10-15-2019

Description: Generate sbatch scripts able to run Augustus for all genomes with poor annotations 
with automatic Augustus model assignment based on closest insect order ancestor
'''

import re, os, sys
import pandas as pd
from optparse import OptionParser

def set_augustus_model(order, table):
    if order == 'Diplura' :
        hmm_order = 'frankliniella_occidentalis'
    elif order == 'Archaeognatha':
        hmm_order = 'frankliniella_occidentalis'
    elif order == 'Odonata' :
        hmm_order = 'zootermopsis_nevadensis'
    elif order == 'Ephemeroptera' :
        hmm_order = 'frankliniella_occidentalis'
    elif order == 'Plecoptera' :
        hmm_order = 'zootermopsis_nevadensis'
    elif order == 'Orthoptera' :
        hmm_order = 'zootermopsis_nevadensis'
    elif order == 'Phasmatodea' :
        hmm_order = 'zootermopsis_nevadensis'
    elif order == 'Blattodea':
        hmm_order = 'zootermopsis_nevadensis'
    elif order == 'Thysanoptera' :
        hmm_order = 'frankliniella_occidentalis'
    elif order == 'Hemiptera' :
        hmm_order = 'bemisia_tabaci'
    elif order == 'Phthiraptera' :
        hmm_order = 'pediculus_humanus'
    elif order == 'Hymenoptera' :
        hmm_order = 'apis_mellifera'
    elif order == 'Strepsiptera' :
        hmm_order = 'tribolium_castaneum'
    elif order == 'Coleoptera' :
        hmm_order = 'tribolium_castaneum'
    elif order == 'Trichoptera' :
        hmm_order = 'papilio_xuthus'
    elif order == 'Lepidoptera' :
        hmm_order = 'papilio_xuthus'
    elif order == 'Siphonaptera' :
        hmm_order = 'ctenophalides_felis'
    elif order == 'Diptera' :
        family = table[table['order_name'] == order ]['family_name']
        if 'Culicidae' in family :
            hmm_order = 'aedes'
        elif 'Pteromalidae' in family :
            hmm_order = 'nasonia'
        else :
            hmm_order = 'fly'
    return hmm_order

def set_metadata(csv_file):
    table = pd.read_csv(csv_file)
    GCA = table[table['genome_id'].str.contains('GCA')]
    GCA['augustus_model'] = GCA['order_name'].apply(set_augustus_model, args=(GCA,))
    return GCA

def return_args(table, genome_id):
    augustus_model = table[table['genome_id'] == genome_id]['augustus_model'].values[0]
    return augustus_model

if __name__ == "__main__":
    parser = OptionParser()
    parser.add_option("-t", "--table", dest='table', default="None",
                      help="[Required] Genome table information")
    parser.add_option("-c", "--config_path", dest='config_path', default="None",
                      help="[Required] Location of config folder for Augustus")
    parser.add_option("-g", "--genome_path", dest='genome_path', default="None",
                      help="[Required] Location of the folder which contain all genomes for augustus analysis")
    parser.add_option("-r", "--result_path", dest='result_path', default="None",
                      help="[Required] Location to store sbatch scripts")

    (options, args) = parser.parse_args()
    table = options.table
    config_path = options.config_path
    genome_path = options.genome_path
    result_path = options.result_path

    if not os.path.isdir(result_path):
        os.mkdir(result_path)

    GCA = set_metadata(table)

    for genome in os.listdir(genome_path):
        genome_id = re.findall(r'(^[A-Z]{3}_[0-9]{9}.[0-9]{1})', genome)[0]
        augustus_model = return_args(GCA, genome_id)

        path = '{}/{}_augustus.sh'.format(result_path, genome_id)

        if not os.path.isfile(path) :
            print('Making sbatch file for {}'.format(genome))
            f = open(path,'w')
            f.write('#!/bin/bash\n#\n')
            f.write('''#SBATCH -J {}_training
#SBATCH -o ./{}_augustus.gff
#SBATCH -e ./{}_augustus.err\n'''.format(genome_id, genome_id, genome_id))
            f.write('''#SBATCH -n 1 # Number of cores
#SBATCH -p serial_requeue  # Partition
#SBATCH --mem 3000 # Memory request (Mo)
#SBATCH -t 7-00:00 # Maximum execution time (D-HH:MM)
#
#
#
module load augustus/3.3-fasrc02
module load boost/1.54.0-fasrc03
module load zlib/1.2.8-fasrc09 
module load bamtools/2.3.0-fasrc01
module load samtools/0.1.19-fasrc02
module load htslib/1.5-fasrc02
module load bcftools/1.5-fasrc02
module load tabix/0.2.6-fasrc01
module load lp_solve/5.5.2.5-fasrc01
module load SuiteSparse/4.2.1-fasrc02
module load sqlite/3081101-fasrc01\n\n''')
            f.write('augustus --strand=both --codingseq=on --protein=on --gff3=on --AUGUSTUS_CONFIG_PATH={}  --species={} {}\n'.format(config_path, augustus_model, os.path.join(genome_path, genome)))
            f.close()

        else :
            ('sbatch file is already created for this genome ({})'.format(genome_id))
