import re as regular_expression
from fractions import Fraction

def _structure_regex():
    return "[<>().]+"

def _sequence_regex():
    return "[ACGU]+"

def _get_matched_text(regex_match):
    return regex_match.group()

class SuboptFile:
    def __init__(self, raw_text):
        self.original_text = raw_text   
        self.__parseSecondaryStructures()
        self.__parseSequence()
        self.__parseParameters()
        self.__parseSignature()
        self.__parseSignatures()

    def __parseSecondaryStructures(self):
        self.structures = []
        for line in self.original_text.split('\n'):
            if line and line[0].isdigit():
                regex_match = regular_expression.search(pattern = _structure_regex(), string = line)
                self.structures.append(_get_matched_text(regex_match))

    def __parseSequence(self):
        sequence_line_index = 2
        sequence_line = self.original_text.split('\n')[sequence_line_index]
        regex_match = regular_expression.search(pattern = _sequence_regex(), string = sequence_line)
        self.sequence = _get_matched_text(regex_match)


    def __parseParameters(self):
        parameters_line_index = 1
        parameters_line = self.original_text.split('\n')[parameters_line_index]
        # prefix, params = parameters_line.split(":",1)
        # print parameters_line.split()
        afrac = parameters_line.split()[4]
        adec = parameters_line.split()[6].rstrip(',')
        bfrac = parameters_line.split()[9]
        bdec = parameters_line.split()[11].rstrip(',')
        cfrac = parameters_line.split()[14]
        cdec = parameters_line.split()[16].rstrip(',')
        self.parameter_afrac = afrac
        self.parameter_bfrac = bfrac
        self.parameter_cfrac = cfrac
        #self.parameter_adec = float(Fraction(adec))
        #self.parameter_bdec = float(Fraction(bdec))
        #self.parameter_cdec = float(Fraction(cdec))
        self.parameter_adec=adec
        self.parameter_bdec=bdec
        self.parameter_cdec=cdec
  

    def __parseSignature(self):
        signature_line_index = 4
        signature_line = self.original_text.split('\n')[signature_line_index]
        # print signature_line.split()
        m = signature_line.split()[2]
        u = signature_line.split()[3]
        b = signature_line.split()[4]
        w = signature_line.split()[5]
        self.signature = [m,u,b,w]

    def __parseSignatures(self):
        self.signatures = []
        for line in self.original_text.split('\n'):
            if line and line[0].isdigit():
                m = line.split()[2]
                u = line.split()[3]
                b = line.split()[4]
                w = line.split()[5]
                signature = [m,u,b,w]
                self.signatures.append(signature)

