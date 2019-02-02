#!/usr/bin/python
import sys
import argparse

parser = argparse.ArgumentParser()
parser.add_argument("-t", "--threshold", dest= 'th', type=int, default=10,
                  help="Filter variants with less than t samples in any genotype group (default 10)", metavar="INT")
parser.add_argument("-g", "--genotype", dest= 'genotype', type=str, default=10,
                  help="Genotype file in VCF format", metavar="FILE")
args = parser.parse_args()

with open(args.genotype, 'r') as rd, sys.stdout as wr :
        for line in rd:
                if line.startswith('#') :
                        wr.write(line)
                elif line.count('0|0') > 0 or line.count('1|1') > 0 or line.count('0|1') > 0 or line.count('1|0') > 0:
                        if line.count('0|0') < args.th and line.count('0|0') >0:
                                pass
                        elif ( line.count('0|1') + line.count('1|0') ) < args.th and ( line.count('0|1') + line.count('1|0') ) >0:
                                pass
                        elif line.count('1|1') < args.th and line.count('1|1') >0:
                                pass
                        elif line.count('1|1') == 0 and line.count('0|1') == 0 and line.count('1|0') == 0:
                                pass
                        elif line.count('1|1') == 0 and line.count('0|0') == 0:
                                pass
                        elif line.count('0|0') == 0 and line.count('0|1') == 0 and line.count('1|0') == 0:
                                pass
                        else:
                                wr.write(line)

                else:
                        if line.count('0/0') < args.th and line.count('0/0') >0:
                                pass
                        elif line.count('0/1') < args.th and line.count('0/1') >0:
                                pass
                        elif line.count('1/1') < args.th and line.count('1/1') >0:
                                pass
                        elif line.count('1/1') == 0 and line.count('0/1') == 0:
                                pass
                        elif line.count('1/1') == 0 and line.count('0/0') == 0:
                                pass
                        elif line.count('0/0') == 0 and line.count('0/1') == 0:
                                pass
                        else:
                                wr.write(line)
