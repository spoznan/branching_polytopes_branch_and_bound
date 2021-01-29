#Behind the scene functions and algorithms for bb_main.sage

from rna_poly import *
from polytope_funs import *
from SuboptFile import *
from NopctFile import *
from accuracy_analysis import *
from FilenameToAccessionNumberTools import *

import subprocess
import csv
import os

## Reads in max accuracy data for sequences 
def FindBounds(AccuracyValuesFileName):
    ## AccuracyValuesFileName.txt format:
        ## Aquifex.aeolicus.VF5_AE000657 = 1
        ## arabidopsis.thaliana_AC006340 = 0.930232558
        ## ...

    AccuracyFile = open(AccuracyValuesFileName)
    RawData = AccuracyFile.readlines()
    AccuracyFile.close()    
    u = []
    for i in RawData:
        u.append(eval(i.split(" = ")[1]))
    
    return u

def InitializePolytopes(AccuracyValuesFileName):
    ## Reads in accuracy data
    AccuracyFile = open(AccuracyValuesFileName)
    RawData = AccuracyFile.readlines()
    AccuracyFile.close()

    PolytopeList = []
    for i in RawData:
        if i.startswith("d.5"):
            print("/projects/rna/rnatope/Data/5S_50/rnapoly/" + i.split(" = ")[0] + ".rnapoly")
            x = RNAPolytope.construct_from_file("./Data/5S_50/rnapoly/" + i.split(" = ")[0] + ".rnapoly")
            PolytopeList.append(x)
        else:
            print("/projects/rna/rnatope/Data/tRNA_50/rnapoly/" + i.split(" = ")[0] + ".rnapoly")
            x = RNAPolytope.construct_from_file("./Data/tRNA_50/rnapoly/" + i.split(" = ")[0] + ".rnapoly")
            PolytopeList.append(x)

    return PolytopeList

def OrderSlices(PolytopeList):
    Slices = []
    for i in range(len(PolytopeList)):
        s = list(PolytopeList[i].d1_slices())

        for i in range(len(s)):
            s[i] = [s[i].representative_point(), s[i]]
        s.sort()
        for i in range(len(s)):
            s[i] = s[i][1]
        
        Slices.append(s)
    return Slices

def compute_slice_accuracy(subopt_file_name, PolytopeName):
        with open(subopt_file_name, 'r') as fileHandle:
                subopt = SuboptFile(fileHandle.read())
                sequence = subopt.sequence
                if PolytopeName.startswith("d.5"):
                    nopct = NopctFile(open("./Data/5S_50/nopct/" + FullSequenceToNopct[PolytopeName + ".sobj"] + ".nopct","r").read())
                else:
                    nopct = NopctFile(open("./Data/tRNA_50/nopct/" + FullSequenceToNopct[PolytopeName + ".sobj"] + ".nopct","r").read())
                F_value_sum = 0
                for struct in subopt.structures:
                        predicted_pairs = dot_bracket_to_base_pairs(struct)
                        true_pairs = nopct.canonicalPairs

                        true_positives = true_pairs.intersection(predicted_pairs)
                        false_negatives = true_pairs.difference(predicted_pairs)
                        false_positives = predicted_pairs.difference(true_pairs)

                        number_noncanonical_pairs = calculate_number_noncanonical_pairs(sequence, true_pairs)
                        Fscore = calculate_F_score(len(true_positives), len(false_negatives), len(false_positives), number_noncanonical_pairs = number_noncanonical_pairs)

                        F_value_sum += Fscore
                F_avg = F_value_sum/len(subopt.structures)
        return F_avg


def AccuracyOf(d1Slice, PolytopeName):
    if PolytopeName.startswith("d.5"):
        rnasuboptfile = "./Data/5S_50/optimal_all_vertices/" + PolytopeName + "/"
    else:
        rnasuboptfile = "./Data/tRNA_50/optimal_all_vertices/" + PolytopeName + "/"
    print(slice_center(d1Slice))
    print(d1Slice.vertices())
    try:
        rnasuboptfile += str(slice_center(d1Slice)[0]).split("/")[0] + ":" + str(slice_center(d1Slice)[0]).split("/")[1] + "_"
    except:
        rnasuboptfile += str(slice_center(d1Slice)[0]) + "_"
    try:
        rnasuboptfile += str(slice_center(d1Slice)[1]).split("/")[0] + ":" + str(slice_center(d1Slice)[1]).split("/")[1] + "_"
    except:
        rnasuboptfile += str(slice_center(d1Slice)[1]) + "_"
    try:
        rnasuboptfile += str(slice_center(d1Slice)[2]).split("/")[0] + ":" + str(slice_center(d1Slice)[2]).split("/")[1]
    except:
        rnasuboptfile += str(slice_center(d1Slice)[2])
    rnasuboptfile += ".rnasubopt"
    print(rnasuboptfile)
    F_avg = compute_slice_accuracy(rnasuboptfile, PolytopeName)
    return F_avg


    
def CreateAccuracyData(AccuracyValuesFileName, PolytopeList, Slices):
    AccuracyFile = open(AccuracyValuesFileName)
    RawData = AccuracyFile.readlines()
    AccuracyFile.close()
    
    AccuracyData = []
    for i in range(len(PolytopeList)):
        AccuracyData.append([])
        for j in Slices[i]:
            AccuracyData[i].append(AccuracyOf(j, RawData[i].split(" = ")[0]))
    return AccuracyData

def CreateNameData(AccuracyValuesFileName):
    AccuracyFile = open(AccuracyValuesFileName)
    RawData = AccuracyFile.readlines()
    AccuracyFile.close()

    NameData = []
    for i in range(len(RawData)):
        NameData.append(RawData[i].split(" = ")[0])
    return(NameData)

def getLValue(filename):
    infile = open(filename)
    rawdata = infile.readlines()
    infile.close()
    return eval(rawdata[0].strip())

