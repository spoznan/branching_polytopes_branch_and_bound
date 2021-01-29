import itertools
import functools

class NopctFile:
	def __init__(self, rawNopct):
        # Remove leading and trailing whitespace
		stripped = rawNopct.strip()
		verify_nopct_format(stripped)

		self.nucleotideIndex = 1
		self.pairIndex = 4


		lines = stripped.split('\n')
		self.fileName = self.parseFileName(lines)
		self.organism = self.parseOrganism(lines)
		self.accessionNumber = self.parseAccessionNumber(lines)
		self.sequence = self.parseSequence(lines)
		canonPairs, noncanonPairs = self.parseCanonAndNoncanonPairs(lines, self.sequence)
		self.canonicalPairs = canonPairs
		self.noncanonPairs = noncanonPairs

	def parseFileName(self, nopctLines):
		fileNameLine = nopctLines[0]
		removePrefix = fileNameLine.replace(fileNamePrefix(), '')
		return removePrefix.strip()

	def parseOrganism(self, nopctLines):
		organismLine = nopctLines[1]
		removePrefix = organismLine.replace(organismPrefix(), '')
		return removePrefix.strip()

	def parseAccessionNumber(self, nopctLines):
		accessionLine = nopctLines[2]
		removePrefix = accessionLine.replace(accessionPrefix(), '')
		return removePrefix.strip()

	def parseSequence(self, nopctLines):
		extractAndAppendNucleotide = lambda accumulator, item: accumulator + item[self.nucleotideIndex]
		data = self.dataLines(nopctLines)
		return functools.reduce(extractAndAppendNucleotide, data, "")


	# returns canonical pairs, noncanonical pairs
	def parseCanonAndNoncanonPairs(self, nopctLines, sequence):
		canonicalPairs = set()
		noncanonPairs = set()

		indexFrom1_index = 0
		for item in self.dataLines(nopctLines):
			if self.isPaired(item):
				# Work with 1-indexed values, because that seems to be how nopct stores where the paired nucleotide is
				indexFrom1 = int(item[indexFrom1_index])
				pairedIndex = int(item[self.pairIndex])
				pair = tuple(sorted([indexFrom1, pairedIndex]))

				nucleotideHere = sequence[indexFrom1 - 1]
				nucleotideThere = sequence[pairedIndex - 1]
				if pairIsCanon(nucleotideHere, nucleotideThere):
					canonicalPairs.add(pair)
				else:
					noncanonPairs.add(pair)
		return canonicalPairs, noncanonPairs

	def isPaired(self, item):
		# item is a list of strings.  Typically, one of the items produced by dataLines
		unpairedSentinelValue = "0"
		return item[self.pairIndex] != unpairedSentinelValue

	def dataLines(self, nopctLines):
		# Deliberately returns a list of lists of strings.
		# Excludes the first four lines, which are metatadata about the file.
		# The data about the sequence seems to start at line 5.
		numberOfMetadataLines = 4
		return map(lambda x: x.split(), nopctLines[numberOfMetadataLines + 1 : ])
		# The format of a nopct data line seems to be:
		# <index from 1> <nucleotide> <index from 0> <(2 + index from 0) mod (sequence length + 1)> <index of paired nucleotide, or 0 if unpaired> <index from 1 again>


	@staticmethod
	def loadFromString(nocptStr):
		return NopctFile(nocptStr)
	@staticmethod
	def loadFromFile(filePath):
		with open(filePath, 'r') as fileHandle:
			return NopctFile(fileHandle.read())


def fileNamePrefix():
	return "Filename:"
def organismPrefix():
	return "Organism:"
def accessionPrefix():
	return "Accession Numbers:"

def pairIsCanon(firstLetter, secondLetter):
	# Validate arguments
	nucleotides = ['a', 'c', 'g', 'u']
	if not firstLetter.lower() in nucleotides:
		raise ValueError("Illegal argument: firstLetter=" + firstLetter)
	if not secondLetter.lower() in nucleotides:
		raise ValueError("Illegal argument: secondLetter=" + secondLetter)
	
	canonical_pairs = ['au', 'cg', 'gu'] # 'au' and 'cg' are Watson-Crick pairs, 'gu' is the so-called 'wobble pair'
	sortedPair = sorted(firstLetter.lower() + secondLetter.lower())
	pairString = ''.join(sortedPair)
	return pairString in canonical_pairs


def verify_nopct_format(nopctString):
	# The nopct files here have a very, very strange structure.
	# I'm not sure the structure is what it should be.
	# This long function lists the assumptions this file makes about the kind of file it is given.
	lines = list(map(lambda s: s.split(), nopctString.split('\n')))

	first_line_header = "Filename:"
	assert(lines[0][0] == first_line_header)
	second_line_header = "Organism:"
	assert(lines[1][0] == second_line_header)
	third_line_header = "Accession"
	assert(lines[2][0] == third_line_header)
	fourth_line_header = "Citation"
	assert(lines[3][0] == fourth_line_header)
	def is_int(x):
	    try:
                int(x)
                return True
	    except ValueError:
                return False
	def is_float_strict(x):
	    try: 
                float(x)
                return not x.isdigit() # isdigit returns false if there is a period
	    except ValueError:
                return False
	assert is_int(lines[4][0])
	assert lines[4][1] == "dG"
	assert lines[4][2] == "="
	assert is_float_strict(lines[4][3])
	assert lines[4][4] == "[initially"
	assert lines[4][5] == "0.0]"
	num_nucleotides = int(lines[4][0])
	assert len(lines) == num_nucleotides + 5, nopctString
	nucleotides = set(['A', 'C', 'G', 'U'])
	for num in range(num_nucleotides):
	    index = num + 5
	    assert lines[index][0] == str(num + 1), str(num)
	    assert lines[index][1] in nucleotides, str(num)
	    assert lines[index][2] == str(num), str(num)
	    assert lines[index][3] == str((num + 2) % (num_nucleotides + 1)), str(num) + "," + str(lines[index][0])
	    bonded_nucleotide = lines[index][4]
	    assert is_int(bonded_nucleotide), str(num)
	    bonded_nucleotide = int(bonded_nucleotide)
	    if bonded_nucleotide > 0:
                other_index = bonded_nucleotide + 5 - 1
                other_bonded = lines[other_index][4]
                assert is_int(other_bonded), str(num)
                other_bonded = int(other_bonded)
                assert other_bonded == num + 1, str(num)
	    assert lines[index][5] == lines[index][0], str(num)
