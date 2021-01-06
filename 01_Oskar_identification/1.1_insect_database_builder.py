#!/usr/bin/env python

'''
Author: Savandara Besse
Created: 03-08-2017
Modified: 09-18-2019

Description: Build the insect genome and transcriptome tables containing taxonomic information, 
bioproject and biosample IDs for the available species in NCBI. 
Also working for crustacean genomes and transcriptomes.
'''

import pandas as pd
from Bio import Entrez
import os, progressbar, sys, threading, time
Entrez.email = 'your@email.here' #Use your own email


def collectGenomeInfo(assemblyId):
    '''
    Collecting global information for genome data such as genome accession number,
    bioproject identifier, biosample identifiers and taxonomy identifier.
    '''
    taxId = ''
    genomeId = ''
    bioproject = ''
    biosample = ''
    species = ''
    time.sleep(0.5)
    handle = Entrez.esummary(db='assembly',id=assemblyId)
    record = Entrez.read(handle,validate=False)['DocumentSummarySet']['DocumentSummary'][0]
    taxId = record['Taxid']
    genomeId = record['AssemblyAccession']
    species = record['SpeciesName']
    bioproject = record['GB_BioProjects'][0]['BioprojectAccn']
    biosample = record['BioSampleAccn']
    return taxId, genomeId, species, bioproject, biosample


def collectTranscriptomeInfo(nuccoreId):
    '''
    Collecting global information for transcriptome data such as transcriptome shotgun assembly identifier,
    bioproject identifier, biosample identifiers, single read archive identifiers, and taxonomy identifier.
    '''
    bioproject = ''
    biosample = []
    sra = []
    species = ''
    taxId = ''
    time.sleep(0.5)
    handle = Entrez.efetch(db='nuccore', id=nuccoreId, retmode='xml')
    record = Entrez.read(handle)[0]
    tsaId = record['GBSeq_accession-version']
    if 'GBSeq_xrefs' in record :
        for Dict in record['GBSeq_xrefs']:
            if 'BioProject' in Dict['GBXref_dbname']:
                bioproject = Dict['GBXref_id']
            if 'BioSample' in Dict['GBXref_dbname']:
                biosample.append(Dict['GBXref_id'])
            if 'Sequence Read Archive' in Dict['GBXref_dbname']:
                sra.append(Dict['GBXref_id'])
    species = record['GBSeq_organism']
    taxHandle = Entrez.esearch(db='taxonomy', term=species)
    taxId = Entrez.read(taxHandle)['IdList'][0]
    if len(biosample) == 0:
        biosample = None
    elif len(biosample) == 1 :
        biosample = biosample[0]
    else:
        biosample = ','.join(biosample)
    if len(sra) == 0:
        sra = None
    elif len(sra) == 1 :
        sra = sra[0]
    else:
        sra = ','.join(sra)
    return tsaId, taxId, species, bioproject, biosample, sra


def collectSpeciesInfo(taxId):
    '''
    Collecting global information such as order and family identifiers
    and names based on the taxonomy identifier of a specific organism
    '''
    familyId = ''
    familyName = ''
    orderId = ''
    orderName = ''
    time.sleep(0.5)
    taxHandle = Entrez.efetch(db="taxonomy", id=taxId)
    record = Entrez.read(taxHandle)
    for index in range(len(record)) :
        lineageEx = record[index]["LineageEx"]
        for Dict in lineageEx :
            if Dict["Rank"] == "family" :
                familyId = Dict["TaxId"]
                familyName = Dict["ScientificName"]
            if Dict["Rank"] == "order" :
                orderId = Dict["TaxId"]
                orderName = Dict["ScientificName"]
    if 'familyId' not in locals():
        familyId = None
    if 'familyName' not in locals():
        familyName = None
    if 'orderId' not in locals():
        orderId = None
    if 'orderName' not in locals():
        orderName = None
    return familyId, familyName, orderId, orderName


def getRequestIds(request):
    if 'Properties' in request:
        time.sleep(0.5)
        getCount = Entrez.esearch(db="nuccore", term=request)
        count = int(Entrez.read(getCount)['Count'])
        getList = Entrez.esearch(db="nuccore", term=request, retmax=count)
        record = Entrez.read(getList)['IdList']
    else:
        time.sleep(0.5)
        getCount = Entrez.esearch(db="assembly", term=request)
        count = int(Entrez.read(getCount)['Count'])
        getList = Entrez.esearch(db="assembly", term=request, retmax=count)
        record = Entrez.read(getList)['IdList']
    return record


def getIds(refseqRequest,genbankRequest,transcriptomeRequest):
    print('NCBI requests in progress...')
    refseqIds = getRequestIds(refseqRequest)
    genbankIds = getRequestIds(genbankRequest)
    genomeIds = refseqIds+genbankIds
    transcriptomeIds = getRequestIds(transcriptomeRequest)
    return genomeIds, transcriptomeIds

def createTable(idList, dataType, toEscape, outFile):
    insectGenomeList = []
    insectTranscriptomeList =[]
    bar = progressbar.ProgressBar(widgets=[' [', progressbar.Timer(), '] ', progressbar.Bar(), ' (', progressbar.ETA(), ')', ])
    if 'genome' in dataType :
        print('Collecting information for genomes...')
        for Id in bar(genomeIds):
            taxId, genomeId, species, bioproject, biosample = collectGenomeInfo(Id)
            if taxId not in toEscape:
                toEscape.append(taxId)
                familyId, familyName, orderId, orderName = collectSpeciesInfo(taxId)
                insectGenomeList.append([genomeId, taxId, species, familyId, familyName, orderId, orderName, bioproject, biosample])
        print('Creating genome database...')
        genomeDataframe = pd.DataFrame(data=insectGenomeList,columns=['genome_id', 'tax_id','species','family_id', 'family_name', 'order_id', 'order_name', 'bioproject', 'biosample'])
        genomeDataframe.to_csv(outFile, sep=',', header=True, index=False)
    elif 'transcriptome' in dataType:
        print('Collecting information for transcriptomes...')
        for Id in bar(transcriptomeIds):
            tsaId, taxId, species, bioproject, biosample, sra = collectTranscriptomeInfo(Id)
            familyId, familyName, orderId, orderName = collectSpeciesInfo(taxId)
            insectTranscriptomeList.append([tsaId, taxId, species, familyId, familyName, orderId, orderName, bioproject, biosample, sra])
        print('Creating transcriptome database...')
        transcriptomeDataframe = pd.DataFrame(data=insectTranscriptomeList,columns=['tsa_id', 'tax_id','species','family_id', 'family_name', 'order_id', 'order_name', 'bioproject', 'biosample', 'sequence_read_archive'])
        transcriptomeDataframe.to_csv(outFile, sep=',', header=True, index=False)

if __name__ == '__main__':
    ### Insects
    refseqRequest="latest+refseq[filter] AND txid50557[Orgn]"
    genbankRequest="latest+genbank[filter] AND txid50557[Orgn]"
    transcriptomeRequest="tsa+master[Properties] AND txid50557[Orgn]"
    ## Collecting Ids from NCBI requests
    genomeIds,transcriptomeIds = getIds(refseqRequest, genbankRequest, transcriptomeRequest)
    ## Insect genome and trnascriptome databases
    createTable(genomeIds, 'genome', [], '../Data/01_Oskar_identification/genome_insect_database.csv')
    createTable(transcriptomeIds, 'transcriptome', [], '../Data/01_Oskar_identification/transcriptome_insect_database.csv')

    ### Crustaceans
    refseqRequest="latest+refseq[filter] AND txid6657[Orgn]"
    genbankRequest="latest+genbank[filter] AND txid6657[Orgn]"
    transcriptomeRequest="tsa+master[Properties] AND txid6657[Orgn]"
    # Collecting Ids from NCBI requests
    genomeIds,transcriptomeIds = getIds(refseqRequest, genbankRequest, transcriptomeRequest)
    ## Insect genome and trnascriptome databases
    createTable(genomeIds, 'genome', [], '../Data/01_Oskar_identification/genome_crustacean_database.csv')
    createTable(transcriptomeIds, 'transcriptome', [], '../Data/01_Oskar_identification/transcriptome_crustacean_database.csv')
