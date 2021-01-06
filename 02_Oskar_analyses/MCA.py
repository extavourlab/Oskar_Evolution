import numpy as np
from Bio import AlignIO
import mca
import pandas as pd
import matplotlib.pyplot as plt


# MCA method for protein MSA decomposition.
# Taken from the paper: Protein interactions and ligand binding: From protein subfamilies to functional specificity
# availabe here https://www.ncbi.nlm.nih.gov/pmc/articles/PMC2808218/bin/0908044107_pnas.0908044107_SI.pdf


# Define the matrix
# Given a MSA of N sequences and L positions, a data matrix
# WN×Q of dimensions N × Q (where Q ¼ 21L) is built, representing
# each position “l” in the alignment as a complete disjunctive
# category with 21 different modalities (representing the 20 amino
# acid types plus the gap), encoding the presence of a modality with
# “1” and its absence with “0.”

# Given the alignement
# ATTCG
# -CTCG

# The matrix will be
#     1A 1C 1G 1T 1- 2A 2C 2G 2T 2- ... 5A 5C 5G 5T 5-
# s1  1  0  0  0  0  0  0  0  1  0      0  0  1  0  0
# s2  0  0  0  0  1  0  1  0  0  0      0  0  1  0  0
# etc ...

class MCA:
    def __init__(self, alignpath, exclude=[]):
        self.alphabet = list("ACDEFGHIKLMNPQRSTVWYX*-")
        align = AlignIO.read(alignpath, 'fasta')
        self.alignment = [i for i in align if i.name not in exclude]
        self.legend = [i.name for i in self.alignment if i.name not in exclude]
        self.N = len(self.alignment)
        self.L = len(self.alignment[0])
        self.Q = self.L * len(self.alphabet) # == 21

    def build_W(self):
        """Given a MSA of N sequences and L positions, a data matrix
        WN×Q of dimensions N × Q (where Q ¼ 21L) is built, representing
        each position “l” in the alignment as a complete disjunctive
        category with 21 different modalities (representing the 20 amino
        acid types plus the gap), encoding the presence of a modality with
        “1” and its absence with “0”."""
        W = np.zeros((self.N, self.Q))
        n = 0
        for seq in self.alignment:
            l = 0
            for pos in seq:
                index = self.alphabet.index(pos)
                q = l * 21 + index
                W[n][q] = 1
                l += 1
            n += 1
        return W

    def build_X(self, W):
        """Columns in W without a '1' are removed for subsequent
        consistency without a loss of generality, resulting in a matrix X"""
        X = []
        self.mapback = {}
        j = 0
        for i in range(W.shape[1]):
            if W[:,i].sum() != 0:
                X.append(W[:,i])
                j += 1
                self.mapback[i] = j
        X = np.array(X).T
        return X

    def MCA(self, X, n=2):
        X = pd.DataFrame(X, index=self.legend)
        mca_ben = mca.MCA(X, ncols=self.L)
        print(mca_ben.inertia, mca_ben.L.sum())
        data = np.array([mca_ben.L[:n],
                         mca_ben.expl_var(greenacre=True, N=n) * 100]).T
        df = pd.DataFrame(data=data, columns=['cλ','%c'], index=range(1,n+1))
        # print(df)
        fs, cos, cont = 'Factor score','Squared cosines', 'Contributions x 1000'
        results = pd.DataFrame(columns=X.index, index=pd.MultiIndex
                              .from_product([[fs, cos, cont], range(1, n+1)]))

        results.loc[fs,    :] = mca_ben.fs_r(N=n).T
        results.loc[cos,   :] = mca_ben.cos_r(N=n).T
        results.loc[cont,  :] = mca_ben.cont_r(N=n).T * 1000

        # print(np.round(results.astype(float), 2))
        self.mca_ben = mca_ben
        self.metadata = df
        return results

    def Plot(self, results):
        fs, cos, cont = 'Factor score','Squared cosines', 'Contributions x 1000'
        points = results.loc[fs].values
        labels = results.columns.values

        plt.figure()
        plt.margins(0.1)
        plt.axhline(0, color='gray')
        plt.axvline(0, color='gray')
        plt.xlabel('Factor 1')
        plt.ylabel('Factor 2')
        plt.scatter(*points, s=120, marker='o', c='r', alpha=.5, linewidths=0)
        for label, x, y in zip(labels, *points):
            plt.annotate(label, xy=(x, y), xytext=(x + .03, y + .03))
        plt.show()

    def run(self):
        W = self.build_W()
        X = self.build_X(W)
        res = self.MCA(X)
        self.Plot(res)


    def run2(self):
        W = self.build_W()
        # X = self.build_X(W)
        res = self.MCA(W)
        self.Plot(res)

mapping_colors = {
'Zygentoma':"#fffbc8",
'Ephemeroptera':"#ffedb4",
'Plecoptera':"#f9cdac",
'Orthoptera':"#f3aca2",
'Phasmatodea':"#ee8b97",
'Blattodea':"#e96a8d",
'Thysanoptera':"#db5087",
'Psocoptera':"#b8428c",
'Hymenoptera':"#973490",
'Coleoptera':"#742796",
'Trichoptera':"#5e1f88",
'Lepidoptera':"#4d1a70",
'Mecoptera':"#3d1459",
'Diptera':"#2d0f41",
'uknw':"#b3b3b3"}


mapping_colors_r = {
'Diptera':"#fffbc8",
'Mecoptera':"#ffedb4",
'Lepidoptera':"#f9cdac",
'Trichoptera':"#f3aca2",
'Coleoptera':"#ee8b97",
'Hymenoptera':"#e96a8d",
'Psocoptera':"#db5087",
'Thysanoptera':"#b8428c",
'Blattodea':"#973490",
'Phasmatodea':"#742796",
'Orthoptera':"#5e1f88",
'Plecoptera':"#4d1a70",
'Ephemeroptera':"#3d1459",
'Zygentoma':"#2d0f41",
'uknw':"#b3b3b3"}


#
#
# import pandas as pd
# import numpy as np
# import mca
#
# alignment = [list("ATTCG"),
#              list("ACTCG"),
#              list("ATTC-"),
#              list("ATT-G"),
#              list("CTTCG"),
#              list("CTT-G"),
#              list("CTCCG"),
#              list("GT-CG"),
#              list("GTTCG"),
#              list("GTTAG")]
# alphabet = list("ATCG-")
# N = len(alignment)
# L = len(alignment[0])
# Q = L * len(alphabet)
# W = np.zeros((N, Q))
# n = 0
# for seq in alignment:
#     l = 0
#     for pos in seq:
#         index = alphabet.index(pos)
#         q = l * len(alphabet) + index
#         W[n][q] = 1
#         l += 1
#     n += 1
# X = []
# for i in range(W.shape[1]):
#     if W[:,i].sum() != 0:
#         X.append(W[:,i])
# X = np.array(X).T
# X = pd.DataFrame(X)
# mca_ben = mca.MCA(X, ncols=L)
# print(mca_ben.inertia, mca_ben.L.sum())
# data = np.array([mca_ben.L[:2],
#                  mca_ben.expl_var(greenacre=True, N=2) * 100]).T
# df = pd.DataFrame(data=data, columns=['cλ','%c'], index=range(1,3))
# print(df)
# fs, cos, cont = 'Factor score','Squared cosines', 'Contributions x 1000'
# table3 = pd.DataFrame(columns=X.index, index=pd.MultiIndex
#                       .from_product([[fs, cos, cont], range(1, 3)]))
#
# table3.loc[fs,    :] = mca_ben.fs_r(N=2).T
# table3.loc[cos,   :] = mca_ben.cos_r(N=2).T
# table3.loc[cont,  :] = mca_ben.cont_r(N=2).T * 1000
#
# print(np.round(table3.astype(float), 2))
#
# points = table3.loc[fs].values
# labels = table3.columns.values
#
# plt.figure()
# plt.margins(0.1)
# plt.axhline(0, color='gray')
# plt.axvline(0, color='gray')
# plt.xlabel('Factor 1')
# plt.ylabel('Factor 2')
# plt.scatter(*points, s=120, marker='o', c='r', alpha=.5, linewidths=0)
# for label, x, y in zip(labels, *points):
#     plt.annotate(label, xy=(x, y), xytext=(x + .03, y + .03))
# plt.show()
#
#
# P = X.shape[1]
# xnS = X.sum(axis=1)
# xSp = X.sum(axis=0)
# xSS = X.sum()
# fnS = xnS / xSS
# fSp = xSp / xSS
# fnp = X / xSS
# Y = np.zeros((N,P))
# Z = np.zeros((N,P))
# for n in range(N):
#     for p in range(P):
#         Y[n][p] = fnp[n][p] / (fSp[p] * np.sqrt(fnS[n]))
#         Z[n][p] = fnp[n][p] / np.sqrt(fSp[p] * fnS[n])
