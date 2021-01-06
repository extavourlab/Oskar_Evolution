#!/usr/bin/env python
import re, os, sys, json
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Alphabet import generic_protein

class Table:
    def __init__(self, msa_path, json_path):
        self.alignment = msa_path
        self.json = json_path

    def upload_data(self):
        with open(self.json, encoding='utf-8') as data_file:
            data = json.loads(data_file.read())
        validated = []
        for seq_record in SeqIO.parse(self.alignment, 'fasta'):
            validated.append(seq_record.id)
        return data,validated

    def build_table(self,data,validated):
        table = {}
        for Id in data :
            if 'Transcriptome_analysis' in data[Id] :
                for analysis in data[Id]['Transcriptome_analysis']:
                    if len(analysis['result']['oskar']['validated']) != 0 :
                        for seq in analysis['result']['oskar']['validated'] :
                            if seq in validated :
                                table[seq] = [data[Id]['metadata']['order']['name'],data[Id]['metadata']['family']['name'],data[Id]['metadata']['species']['name']]

            if 'Genome_analysis' in data[Id] :
                for analysis in data[Id]['Genome_analysis']:
                    if len(analysis['result']['oskar']['validated']) != 0 :
                        for seq in analysis['result']['oskar']['validated'] :
                            if seq in validated :
                                table[seq] = [data[Id]['metadata']['order']['name'],data[Id]['metadata']['family']['name'],data[Id]['metadata']['species']['name']]
        return table
