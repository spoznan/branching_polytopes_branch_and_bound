import csv
import collections
import sys
import os

#AccuracyResult = collections.namedtuple("AccuracyResult", "accession_number organism sequence structure score_stats")

#ScoreStats = collections.namedtuple("ScoreStats", "predicted_pairs true_pairs true_positives false_negatives false_positives")

##def write_results_to_csv(output_file_name, nopct_files_directory, subopt_files_directory):
##    with open(output_file_name, "ab") as file_handle:
##        # "ab" = append, open-in-binary.  Use append so an accident won't destroy an existing file.  Open-in-binary because that's what csv.writer expects.
##        writer = csv.writer(file_handle)
##        writer.writerow(["Sequence name",
##                         "Accession number", 
##                         "number of true positives", 
##                         "number of false positives", 
##                         "number of false negatives", 
##                         "number of noncanonical pairs", 
##                         "F score (canonical pairs only)", "comments"])
##
##
##        results = make_accuracy_statistics(nopct_files_directory, subopt_files_directory)
##
##        accession_to_multiplicity = collections.defaultdict(int)
##        for accession_number,_,  _, _ , _ in results:
##            accession_to_multiplicity[accession_number] += 1
##
##        accession_to_times_encountered = collections.defaultdict(int)   
##        for accession_number, organism, sequence, structure, score_stats in results:
##            number_true_positives = len(score_stats.true_positives)
##            number_false_positives = len(score_stats.false_positives)
##            number_false_negatives = len(score_stats.false_negatives)
##            number_noncanonical_pairs = calculate_number_noncanonical_pairs(sequence, score_stats.true_pairs)
##            Fscore = calculate_F_score(number_true_positives, number_false_negatives, number_false_positives, number_noncanonical_pairs = number_noncanonical_pairs)
##            row =[]
##            row.append(organism)
##            row.append(accession_number)
##            # row.append(structure)
##            row.append(number_true_positives)
##            row.append(number_false_positives)
##            row.append(number_false_negatives)
##            row.append(number_noncanonical_pairs)
##            row.append(Fscore)
##            writer.writerow(row)
##
##
##   
##    
##
##def make_accuracy_statistics(nopct_files_directory, subopt_files_directory):
##    sequence_to_nopct = create_sequence_to_nopct_dict_from_directory(nopct_files_directory)
##    subopt_files = create_subopt_files_from_directory(subopt_files_directory)
##
##    results = []
##    for subopt in subopt_files:
##        sequence = subopt.sequence
##        nopct = sequence_to_nopct[sequence]
##        for struct in subopt.structures:
##            score_stats = calculate_score_stats(struct, nopct)
##            results.append(build_accuracy_result(nopct_file = nopct, structure = struct, score_stats = score_stats))
##    results.sort(key=lambda accResult: accResult.accession_number)
##    return results
##    
##def create_subopt_files_from_directory(subopt_files_directory):
##    subopt_files = []
##    for basename in os.listdir(subopt_files_directory):
##        filename = os.path.join(subopt_files_directory, basename)
##        with open(filename, 'r') as fileHandle:
##            subopt_files.append(SuboptFile(fileHandle.read()))
##    return subopt_files
##
##def create_sequence_to_nopct_dict_from_directory(nopct_files_directory):
##    sequence_to_nopct = {}
##    for basename in os.listdir(nopct_files_directory):
##        filename = os.path.join(nopct_files_directory, basename)
##        with open(filename, 'r') as fileHandle:
##            parsedFile = NopctFile(fileHandle.read())
##            sequence_to_nopct[parsedFile.sequence] = parsedFile
##    return sequence_to_nopct
##
##
##def calculate_score_stats(predicted_dot_bracket, nopct_file):
##    predicted_pairs = dot_bracket_to_base_pairs(predicted_dot_bracket)
##    true_pairs = nopct_file.canonicalPairs
##    
##    true_positives = true_pairs.intersection(predicted_pairs)
##    false_negatives = true_pairs.difference(predicted_pairs)
##    false_positives = predicted_pairs.difference(true_pairs)
##    
##    return ScoreStats(predicted_pairs = predicted_pairs, 
##                      true_pairs = true_pairs, 
##                      true_positives = true_positives, 
##                      false_negatives = false_negatives, 
##                      false_positives = false_positives)
##
##
##def build_accuracy_result(nopct_file, structure,  score_stats):
##    return AccuracyResult(accession_number = nopct_file.accessionNumber, organism = nopct_file.organism, sequence = nopct_file.sequence, structure = structure, score_stats = score_stats)

## USED IN compute_slice_accuracy < AccuracyOf < CreateAccuracyData
def calculate_F_score(number_true_positives, number_false_negatives, number_false_positives, number_noncanonical_pairs = 0):
    precision = number_true_positives / float(number_true_positives + number_false_positives)
    recall = number_true_positives / float(number_true_positives + number_false_negatives - number_noncanonical_pairs)
    if precision == 0 or recall == 0:
        return 0
    else:
        return 2 * precision * recall / (precision + recall)           

##def nopct_fields_to_base_pairs(nopct_fields):
##    output = set()
##    NUCLEOTIDE_IS_UNPAIRED = 0
##    for index_from_1, _, paired_index_from_1 in nopct_fields:
##        if paired_index_from_1 != NUCLEOTIDE_IS_UNPAIRED:
##            assert(index_from_1 != paired_index_from_1)
##            idx_here = index_from_1 - 1
##            pair_idx = paired_index_from_1 - 1
##            assert(nopct_fields[pair_idx].paired_index_from_1 == index_from_1)
##            if index_from_1 < paired_index_from_1:
##                pair = (index_from_1, paired_index_from_1)
##                output.add(pair)
##    return output
            
## USED IN compute_slice_accuracy < AccuracyOf < CreateAccuracyData
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

## USED IN compute_slice_accuracy < AccuracyOf < CreateAccuracyData
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
