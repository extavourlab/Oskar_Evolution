#!/usr/bin/env python

'''
Author: Savandara Besse and Leo Blondel

Created date: 06-28-2017
Modified date: 07-05-2017

Description: Calculate conservation score for an alignment with Valdar algorithms

'''
from optparse import OptionParser
from Bio import AlignIO
import pandas as pd
import sys
import numpy as np
import progressbar


class Scores:
    def __init__(self, path):
        self.alignment = AlignIO.read(path,'fasta')
        self.N = len(self.alignment)
        self.valdar_table = self.get_valdar_table()
        self.hydrophobicity_scale = self.get_hydrophobicity_scale()
        self.W = None

    def distance_seq(self,si,sj):
        Pi = 0.0
        if len(si) != len(sj):
            return False
        for letter in range(len(si)):
            if si[letter] == sj[letter]:
                Pi += 1
        Pi /= len(si)
        return Pi

    def weight(self,i):
        w = 0.0
        for j in range(self.N):
            if j != i :
                w += self.distance_seq(list(self.alignment[i].seq),list(self.alignment[j].seq))
        w /= (self.N-1)
        return w

    def get_weight_sequences(self):
        print("Calculating Weights...")
        W = []
        for i in range(self.N):
            W.append(self.weight(i))
        self.W = W
        return W

    def get_scale_valdar(self,position):
        l = 0.0
        for i in range(self.N):
            for j in range(i+1, self.N):
                Six = self.alignment[i].seq[position]
                Sjx = self.alignment[j].seq[position]
                if (Six != ' ') and (Sjx != ' '):
                    l += (self.W[i] * self.W[j])
        l = 1.0/l
        return l

    def Score_Valdar(self,position):
        M = 0.0
        SP = 0.0
        valdar = 0.0
        if not self.W:
            self.get_weight_sequences()
        l = self.get_scale_valdar(position)
        for i in range(self.N):
            for j in range(i+1, self.N):
                Six = self.alignment[i].seq[position]
                Sjx = self.alignment[j].seq[position]
                if (Six == ' ') or (Sjx == ' '):
                    SP += 0
                elif (Six == '-') or (Sjx == '-'):
                    SP += 0
                elif (Six == 'X') or (Sjx == 'X'):
                    SP += 0
                elif (Six == '*') or (Sjx == '*'):
                    SP += 0
                else:
                    M = self.valdar_table[Six][Sjx]
                    SP += self.W[i] * self.W[j] * M
        valdar = SP * l
        return valdar

    def Score_Elec(self, position, weighted=False):
        "This Calculate the electrostatic conservation score."
        S = 0 #Define the score variable
        #Define the amino acids of interest and their scores
        scoring_table = {
            'D': -1,
            'E': -1,
            'K': 1,
            'R': 1,
            'H': 0.5
        }
        # If we want to weight the alignment Calculate the scores if they don't exists
        if weighted and not self.W:
            self.get_weight_sequences()
        # Calculate the scaling factor for this position
        if weighted:
            l = self.get_scale_valdar(position)
        #Iterate over all the sequences
        for i in range(self.N):
            #Get the letter at position position for the ith sequence
            letter = self.alignment[i].seq[position]
            #Check if letter is an amino acid of interest and assign the score
            if letter in scoring_table:
                # If we want to weight the alignment, normalize by the weight factor of this sequence
                if weighted:
                    S += scoring_table[letter] * self.W[i]
                else:
                    S += scoring_table[letter]
        if weighted:
            S *= l
        else:
            S /= self.N
        return S

    def Score_Hydro(self, position, weighted=False):
        "This Calculate the hydrophobicity conservation at a given position of the alignment"
        # Define the score at this position
        S = 0
        # If we want to weight the alignment Calculate the scores if they don't exists
        if weighted and not self.W:
            self.get_weight_sequences()
        # Calculate the scaling factor for this position
        if weighted:
            l = self.get_scale_valdar(position)
        # Iterate of all sequences
        for i in range(self.N):
            #Get the letter of sequence i at position
            letter = self.alignment[i].seq[position]
            #Add the hydrophobicity score from scoring table
            if weighted:
                S += self.hydrophobicity_scale[letter] * self.W[i]
            else:
                S +=  self.hydrophobicity_scale[letter]
        if weighted:
            S *= l
        else:
            S /= self.N
        return S

    def get_Score(self, score_type, weighted=False):
        # Score_type is a list with all score types ex: ['valdar','elec']
        scores = {}
        nbCol = self.alignment.get_alignment_length()
        bar = progressbar.ProgressBar(widgets=[' [', progressbar.Timer(), '] ', progressbar.Bar(), ' (', progressbar.ETA(), ') ', ])
        print("Computing Scores")
        for position in bar(range(nbCol)):
            bar.update(position)
            if 'valdar' in score_type or score_type == 'all':
                if 'valdar' not in scores:
                    scores['valdar'] = []
                scores['valdar'].append(self.Score_Valdar(position))
            if 'elec' in score_type or score_type == 'all':
                if 'elec' not in scores:
                    scores['elec'] = []
                scores['elec'].append(self.Score_Elec(position, weighted))
            if 'hydro' in score_type or score_type == 'all':
                if 'hydro' not in scores:
                    scores['hydro'] = []
                scores['hydro'].append(self.Score_Hydro(position, weighted))
        return scores

    def save_Score(self, score_type, outpath, weighted=False):
        results = self.get_Score(score_type, weighted)
        res = []
        colums = []
        for key in results:
            res.append(results[key])
            colums.append(key)
        res = np.array(res).T
        scores = pd.DataFrame(data=res,columns=colums)
        scores.to_csv(outpath)

    def get_valdar_table(self):
        "Get the Valdar Table as described in Valdar 2002"
        values = [
           [ 1.   ,  0.267,  0.333,  0.267,  0.267,  0.267,  0.267,  0.4  , 0.2  ,  0.333,  0.267,  0.267,  0.267,  0.133,  0.4  ,  0.4  , 0.467,  0.067,  0.133,  0.4  ],
           [ 0.267,  1.   ,  0.333,  0.267,  0.267,  0.467,  0.333,  0.333, 0.467,  0.133,  0.133,  0.6  ,  0.2  ,  0.067,  0.267,  0.267, 0.267,  0.333,  0.2  ,  0.133],
           [ 0.333,  0.333,  1.   ,  0.467,  0.267,  0.333,  0.4  ,  0.333, 0.4  ,  0.2  ,  0.133,  0.4  ,  0.2  ,  0.133,  0.267,  0.4  , 0.4  ,  0.067,  0.267,  0.2  ],
           [ 0.267,  0.267,  0.467,  1.   ,  0.133,  0.333,  0.6  ,  0.4  , 0.333,  0.133,  0.067,  0.333,  0.133,  0.   ,  0.2  ,  0.333, 0.267,  0.   ,  0.2  ,  0.133],
           [ 0.267,  0.267,  0.267,  0.133,  1.   ,  0.133,  0.067,  0.267, 0.333,  0.2  ,  0.133,  0.133,  0.2  ,  0.333,  0.2  ,  0.4  , 0.267,  0.4  ,  0.467,  0.2  ],
           [ 0.267,  0.467,  0.333,  0.333,  0.133,  1.   ,  0.467,  0.267, 0.533,  0.133,  0.2  ,  0.467,  0.2  ,  0.067,  0.333,  0.267, 0.267,  0.133,  0.267,  0.133],
           [ 0.267,  0.333,  0.4  ,  0.6  ,  0.067,  0.467,  1.   ,  0.4  , 0.333,  0.133,  0.067,  0.4  ,  0.133,  0.   ,  0.2  ,  0.267, 0.267,  0.   ,  0.067,  0.2  ],
           [ 0.4  ,  0.333,  0.333,  0.4  ,  0.267,  0.267,  0.4  ,  1.   , 0.2  ,  0.133,  0.067,  0.267,  0.133,  0.   ,  0.267,  0.4  , 0.333,  0.2  ,  0.067,  0.2  ],
           [ 0.2  ,  0.467,  0.4  ,  0.333,  0.333,  0.533,  0.333,  0.2  , 1.   ,  0.133,  0.2  ,  0.4  ,  0.2  ,  0.333,  0.333,  0.267, 0.267,  0.133,  0.6  ,  0.133],
           [ 0.333,  0.133,  0.2  ,  0.133,  0.2  ,  0.133,  0.133,  0.133, 0.133,  1.   ,  0.467,  0.133,  0.533,  0.333,  0.2  ,  0.267, 0.4  ,  0.067,  0.2  ,  0.6  ],
           [ 0.267,  0.133,  0.133,  0.067,  0.133,  0.2  ,  0.067,  0.067, 0.2  ,  0.467,  1.   ,  0.133,  0.533,  0.467,  0.333,  0.2  , 0.267,  0.2  ,  0.267,  0.467],
           [ 0.267,  0.6  ,  0.4  ,  0.333,  0.133,  0.467,  0.4  ,  0.267, 0.4  ,  0.133,  0.133,  1.   ,  0.2  ,  0.   ,  0.2  ,  0.267, 0.267,  0.133,  0.133,  0.133],
           [ 0.267,  0.2  ,  0.2  ,  0.133,  0.2  ,  0.2  ,  0.133,  0.133, 0.2  ,  0.533,  0.533,  0.2  ,  1.   ,  0.333,  0.2  ,  0.267, 0.333,  0.133,  0.133,  0.467],
           [ 0.133,  0.067,  0.133,  0.   ,  0.333,  0.067,  0.   ,  0.   , 0.333,  0.333,  0.467,  0.   ,  0.333,  1.   ,  0.2  ,  0.2  , 0.2  ,  0.267,  0.667,  0.333],
           [ 0.4  ,  0.267,  0.267,  0.2  ,  0.2  ,  0.333,  0.2  ,  0.267, 0.333,  0.2  ,  0.333,  0.2  ,  0.2  ,  0.2  ,  1.   ,  0.4  , 0.4  ,  0.   ,  0.133,  0.267],
           [ 0.4  ,  0.267,  0.4  ,  0.333,  0.4  ,  0.267,  0.267,  0.4  , 0.267,  0.267,  0.2  ,  0.267,  0.267,  0.2  ,  0.4  ,  1.   , 0.4  ,  0.133,  0.267,  0.267],
           [ 0.467,  0.267,  0.4  ,  0.267,  0.267,  0.267,  0.267,  0.333, 0.267,  0.4  ,  0.267,  0.267,  0.333,  0.2  ,  0.4  ,  0.4  , 1.   ,  0.067,  0.133,  0.333],
           [ 0.067,  0.333,  0.067,  0.   ,  0.4  ,  0.133,  0.   ,  0.2  , 0.133,  0.067,  0.2  ,  0.133,  0.133,  0.267,  0.   ,  0.133, 0.067,  1.   ,  0.333,  0.067],
           [ 0.133,  0.2  ,  0.267,  0.2  ,  0.467,  0.267,  0.067,  0.067, 0.6  ,  0.2  ,  0.267,  0.133,  0.133,  0.667,  0.133,  0.267, 0.133,  0.333,  1.   ,  0.133],
           [ 0.4  ,  0.133,  0.2  ,  0.133,  0.2  ,  0.133,  0.2  ,  0.2  , 0.133,  0.6  ,  0.467,  0.133,  0.467,  0.333,  0.267,  0.267, 0.333,  0.067,  0.133,  1.   ]]
        columns = ['A', 'R', 'N', 'D', 'C', 'Q', 'E', 'G', 'H', 'I', 'L', 'K', 'M', 'F', 'P', 'S', 'T', 'W', 'Y', 'V']
        M = pd.DataFrame(values, columns=columns, index=columns)
        return M

    def get_hydrophobicity_scale(self):
        hydro = {
            'I':-1.56,
            'V':-0.78,
            'L':-1.81,
            'F':-2.20,
            'C':0.49,
            'M':-0.76,
            'A':0,
            'G':1.72,
            'T':1.78,
            'S':1.83,
            'W':-0.38,
            'Y':-1.09,
            'P':-1.52,
            'H':4.76,
            'E':1.64,
            'Q':3.01,
            'D':2.95,
            'N':3.47,
            'K':5.39,
            'R':3.71,
            'X':0,
            '-':0,
            '*':0
        }

        return hydro



if __name__ == '__main__':
    parser = OptionParser()
    parser.add_option("-a", "--alignment_path", dest='alignment_path', default="None",
                        help="[Required] Location of protein alignment in fasta format")
    parser.add_option("-v", "--valdar", dest='valdar', action="store_true", default=False,
                        help="[Optional] Compute valdar score")
    parser.add_option("-e", "--electro", dest='electro', action="store_true", default=False,
                        help="[Optional] Compute electrostatic score")
    parser.add_option("-y", "--hydro", dest='hydro', action="store_true", default=False,
                        help="[Optional] Compute hydrophobicity score")
    parser.add_option("-w", "--weighted", dest='weighted', action="store_true", default=False,
                        help="[Optional] Compute weighted scores for electrostatic and hydrophobicity scores")
    parser.add_option("-o", "--output_name", dest='output', default="None",
                        help="[Required] Name of your output file")

    (options, args) = parser.parse_args()

    path = options.alignment_path
    valdar = options.valdar
    electro = options.electro
    hydro = options.hydro
    weighted = options.weighted
    if not(valdar) and not(electro) and not(hydro):
        score_type = ['valdar','elec','hydro']
    else:
        score_type = []
        if valdar:
            score_type.append('valdar')
        if electro:
            score_type.append('elec')
        if hydro :
            score_type.append('hydro')

    output = options.output

    scores = Scores(path)
    scores.save_Score(score_type, output, weighted)
