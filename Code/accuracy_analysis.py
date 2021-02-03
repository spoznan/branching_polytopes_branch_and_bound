import csv
import collections
import sys
import os

def calculate_F_score(number_true_positives, number_false_negatives, number_false_positives, number_noncanonical_pairs = 0):
    precision = number_true_positives / float(number_true_positives + number_false_positives)
    recall = number_true_positives / float(number_true_positives + number_false_negatives - number_noncanonical_pairs)
    if precision == 0 or recall == 0:
        return 0
    else:
        return 2 * precision * recall / (precision + recall)           

def dot_bracket_to_base_pairs(dot_bracket):
    # Returns a set of base pars (i,j), where i < j.  Assumes given dot-bracket notation is valid.
    output = set()
    stack = list()
    peek = lambda: stack[-1]
    push = lambda x: stack.append(x)
    pop = lambda: stack.pop()
    empty = lambda: len(stack) == 0
    for i in range(len(dot_bracket)):
        char = dot_bracket[i]
        if char == '(' or char == ')':
            tple = (char, i)
            if empty():
                push(tple)
            else:
                top = peek()
                if tple[0] == top[0]:
                    push(tple)
                else:
                    pop()
                    pair = (top[1] + 1, tple[1] + 1) #.nopct files are 1-indexed, not 0-indexed like python
                    output.add(pair)
    return output

def calculate_number_noncanonical_pairs(sequence, pairs_1_indexed):
    canonical_pairs = ["au", "cg", "gu"]
    output = 0
    for first, second in pairs_1_indexed:
        first -= 1 # Go from 1-indexing to 0-indexing
        second -= 1 # Go from 1-indexing to 0-indexing
        #print(sequence)
        # nucleotides_list = sorted(sequence[first].lower().extend(sequence[second].lower()))
        nucleotides_list = sorted(sequence[first].lower() + sequence[second].lower())
        base_pair = "".join(nucleotides_list)
        if base_pair not in canonical_pairs:
            output += 1
    return output
