#!/usr/bin/env python

'''
Author: Savandara Besse

Created date: 04-13-2017
Modified date: 11-20-2017

Description: - Generate sbatch scripts for running Augustus training and create the input folders needed for 
Augustus models with genomic fasta and gff files of only well-annotated sequences

'''

import pandas as pd
import os, sys, re, subprocess, threading
from subprocess import Popen
from optparse import OptionParser


class modelFolder:
    def __init__(self, genomeFolder, augustusFolderPath):
        self.genomeFolder = genomePath
        self.augustusFolderPath = augustusFolderPath
        self.refseqGenomeList = []


    def getRefseqGenomeList(self):
        for genome in os.listdir(self.genomeFolder):
            if 'GCF' in genome:
                self.refseqGenomeList.append(genome)
        return self.refseqGenomeList


    def getFile(self, path, genome):
        if os.path.isfile(path):
            if not os.path.isdir(os.path.join(self.augustusFolderPath, genome)):
                os.mkdir(os.path.join(self.augustusFolderPath, genome))
            os.system('cp {} {}'.format(path, os.path.join(self.augustusFolderPath, genome)))

    def extractFiles(self,folder):
        contents = os.listdir(os.path.join(self.genomeFolder, folder))
        for element in contents:
            if '_genomic.fna.gz' in element:
                if '_rna_from_genomic.fna.gz' in element:
                    pass
                elif '_cds_from_genomic.fna.gz' in element:
                    pass
                else:
                    path = os.path.join(genomePath, folder, element)
                    self.getFile(path, folder)
            elif '_genomic.gff.gz' in element:
                path = os.path.join(genomePath, folder, element)
                self.getFile(path, folder)


def unzipFile(path):
    if os.path.isfile(path) :
        P = Popen(['gunzip',path])
        ret = P.wait()
        if ret != 0:
            print("Error Gunzipping !")
    path = path.replace('.gz','')
    return path


def writeSbatchFile(augustusFolderPath, name, folder, model):
    augustus_path = '/the/path/to/augustus/install/folder'
    path = '/the/path/to/working/folder'
    config_path = '/the/path/to/augustus/config/files'

    f = open('{}/SLURM/{}_training.sh'.format(augustusFolderPath,name),'w')
    f.write('#!/bin/bash\n#\n')
    f.write('#SBATCH -J {}_training \n#SBATCH -o {}_training.out\n#SBATCH -e {}_training.err\n'.format(name,name,name))
    f.write('''#SBATCH -n 1 # Number of cores
            #SBATCH -p serial_requeue  # Partition
            #SBATCH --mem 2000 # Memory request
            #SBATCH -t 7-00:00 # Maximum execution time (D-HH:MM)
            #
            #
            #
            module load augustus/3.0.3-fasrc02
            module load boost/1.55.0-fasrc01
            module load bamtools/2.3.0-fasrc01
            module load samtools/0.1.19-fasrc01
            module load bcftools/1.0-fasrc01
            module load htslib/1.1-fasrc01
            module load zlib/1.2.8-fasrc02
            module load tabix/0.2.6-fasrc01
            ''')
    f.write('{}/scripts/gff2gbSmallDNA.pl {}/{}{} {}/{}{} 5000 {}/{}/{}.genes.gb\n'.format(augustus_path, path,folder,model['gff'],path,folder,model['fasta'],path,folder,folder))
    f.write('{}/scripts/randomSplit.pl {}/{}/{}.genes.gb 100\n'.format(augustus_path, path,folder,folder))
    f.write('{}/scripts/new_species.pl --species={} --AUGUSTUS_CONFIG_PATH={}\n'.format(augustus_path, name,config_path))
    f.write('{}/bin/etraining --species={} --AUGUSTUS_CONFIG_PATH={} {}/{}/{}.genes.gb.train\n'.format(augustus_path, name,config_path,path,folder,folder))
    f.write('{}/bin/augustus --species={} --AUGUSTUS_CONFIG_PATH={} {}/{}/{}.genes.gb.test\n'.format(augustus_path, name,config_path,path,folder,folder))
    f.close()

if __name__ == '__main__':
    parser = OptionParser()
    parser.add_option("-g", "--genomePath", dest="genomePath", default="None",
                        help="[Required] Location of the directory where are genome folders")
    parser.add_option("-t", "--genomeTable", dest="genomeTable", default="None",
                        help="[Required] Location of the genome database table")
    (options, args) = parser.parse_args()
    genomePath = options.genomePath
    genomeTable = options.genomeTable
    if genomePath == "None":
        print("Genomes folder path must be provided.\n -h for more information")
        sys.exit(1)
    genomeTable = pd.read_csv(genomeTable)

    if not os.path.isdir('./augustus_training/'):
        os.mkdir('./augustus_training/')
    augustusFolderPath = os.path.join('./augustus_training')

    models = modelFolder(genomePath, augustusFolderPath)
    models.getRefseqGenomeList()
    for genome in models.refseqGenomeList:
        print('Building {} folder...'.format(genome))
        models.extractFiles(genome)

    for folder in os.listdir(augustusFolderPath):
        model = {}
        for gzFile in os.listdir(os.path.join(augustusFolderPath, folder)):
            path = os.path.join(augustusFolderPath, folder, gzFile)
            path = unzipFile(path)
            if '_genomic.gff' in path:
                gff = re.findall('/[A-Z]{3}_[0-9]{9}.[0-9]{1}_[A-Za-z0-9._]+_genomic.gff$',path)[0]
                model['gff'] = gff
            else:
                fasta = re.findall('/[A-Z]{3}_[0-9]{9}.[0-9]{1}_[A-Za-z0-9._]+_genomic.fna$',path)[0]
                model['fasta'] = fasta
        for index in genomeTable.index :
            if genomeTable['genome'][index] == folder :
                name = genomeTable['species'][index]
                name = name.lower().replace(' ','_')
        if not os.path.isdir(os.path.join(augustusFolderPath,'SLURM')):
            os.mkdir(os.path.join(augustusFolderPath,'SLURM'))
        writeSbatchFile(augustusFolderPath,name, folder, model)
    print('SBATCH files written')
