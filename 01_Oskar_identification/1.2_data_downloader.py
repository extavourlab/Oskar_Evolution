#!/usr/bin/env python

'''
Author: Savandara Besse

Created date: 03-14-2017
Modified date: 10-06-2019

Description: This script help to download insect genomes and trannscriptome folders by using
accession numbers registered in the insect genome and transcriptome tables built with
the script 1.1_insect_database_builder.py
'''

import pandas as pd
import os, multiprocessing, progressbar, re, time
from ftplib import FTP
from optparse import OptionParser

class Download:
    def __init__(self, table):
        self.table = pd.read_csv(table)
        self.isGenome = False
        self.accessionNumberList = []
        self.urlList = []


    def setTypeTable(self):
        if 'tsa_id' in self.table.columns:
            self.accessionNumberList = list(self.table.tsa_id)
        elif 'genome_id' in self.table.columns:
            self.isGenome = True
            self.accessionNumberList = list(self.table.genome_id)


    def setUrlList(self,List):
        print('Getting data URLs...')
        bar = progressbar.ProgressBar()
        ftp = FTP('ftp.ncbi.nih.gov')
        ftp.login()
        if self.isGenome:
            for accessionNumber in bar(self.accessionNumberList):
                arg = re.findall('^([A-Z]{3})_([0-9]{3})([0-9]{3})([0-9]{3}).[0-9]{1}',accessionNumber)[0]
                url = '/genomes/all/{}/{}/{}/{}/'.format(arg[0],arg[1],arg[2],arg[3])
                try:
                    ftp.cwd(url)
                    for assembly in ftp.nlst():
                        if accessionNumber in assembly:
                            url +='{}/'.format(assembly)
                            self.urlList.append(url)
                except:
                    List.append(accessionNumber)
        else:
            for accessionNumber in bar(self.accessionNumberList):
                arg = re.findall('^([A-Z]{2})([A-Z]{2})[0-9]{8}.([0-9]{1}$)',accessionNumber)[0]
                tsaId = '{}{}0{}'.format(arg[0],arg[1],arg[2])
                try:
                    tmp = '/sra/wgs_aux/{}/{}/{}/'.format(arg[0], arg[1], tsaId)
                    ftp.cwd(tmp)
                    url = '/sra/wgs_aux/{}/{}/{}/{}.1.fsa_nt.gz'.format(arg[0],arg[1],tsaId,tsaId)
                    self.urlList.append(url)
                except:
                    List.append(accessionNumber)
        return List


    def getData(self, folder):
        bar = progressbar.ProgressBar()
        for url in bar(self.urlList):
            if self.isGenome:
                folderName = re.findall('([A-Z]{3}_[0-9]{9}.[0-9]{1})',url)[0]
                path = os.path.join(folder,folderName)
                if not os.path.isdir(path):
                    os.mkdir(path)
            else:
                path = os.path.join(folder)
            print('rsync --copy-links --recursive --times --verbose rsync://ftp.ncbi.nlm.nih.gov{} {}'.format(url,path))
            os.system('rsync --copy-links --recursive --times --verbose rsync://ftp.ncbi.nlm.nih.gov{} {}'.format(url,path))

    def add_download_status(self, table, notFound):
        def isDownloaded(x, notFound):
            if x not in notFound :
                return True
            return False
        DF = pd.read_csv(table)
        if self.isGenome :
            DF['download_status'] = [ isDownloaded(DF['genome_id'][i], notFound) for i in range(len(DF)) ]
        else :
            DF['download_status'] = [ isDownloaded(DF['tsa_id'][i], notFound) for i in range(len(DF)) ]
        DF.to_csv(table, sep=',', header=True, index=False)

if __name__ == '__main__':
    parser = OptionParser()
    parser.add_option("-d", "--download_status", dest="download_status", default=False,
                      help="[Required] Create updated table with download status")
    parser.add_option("-t", "--table", dest="table", default="None",
                      help="[Required] Location of the genome_X_database.csv or transcriptome_X_database.csv")
    parser.add_option("-f", "--folder", dest="path", default="None",
                      help="[Required] Location of the directory where the data will be downloaded")
    (options, args) = parser.parse_args()
    
    download_status = options.download_status
    table = options.table
    folder = options.path

    path = os.path.join(folder)
    if not os.path.isdir(path):
        os.mkdir(path)
    downloader = Download(table)
    downloader.setTypeTable()
    notFound = downloader.setUrlList([])
    if download_status:
        downloader.add_download_status(table, notFound)
    downloader.getData(path)
    print('Not downloaded:', notFound, len(notFound))
