#!/usr/bin/env python
# 2017.05.01 23:00
# this program was wrote for greping lines from a file based on certain gene names that provided in another file.

import argparse
import csv
import re


def _main():
    greper=argparse.ArgumentParser(description="Greping lines from a file (search file) based on certain gene names that provided in another file (match file). Please set the Nth field you want to match by -f, ex: -f 2.")
    greper.add_argument('searchFile',type=argparse.FileType('r'))
    greper.add_argument('matchFile',type=argparse.FileType('r'))
    greper.add_argument('-f',action="store",dest="f",help="Set which field in match file containing gene names.")
    greper.add_argument('outputFile',type=argparse.FileType('w'))
    options=greper.parse_args()

    match=csv.reader(options.matchFile,delimiter='\t',quotechar=None)
    search=csv.reader(options.searchFile,delimiter='\t',quotechar=None)
    output=csv.writer(options.outputFile,delimiter='\t',quotechar=None)

    dic=[]
    for row in search:
        if "@" not in row[0]:
#            print row
            dic.append(row[0])

    for element in match:
        for gene in dic:
            if element[int(options.f)] == gene:
                output.writerow(element)

if __name__ == '__main__':
    _main()
