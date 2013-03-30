# Read_TargetScan.py
#
#
#

import csv
import os
import xml.etree.ElementTree as ET

print "Begin"

mirlist_file = open("/Users/nithintumma/Desktop/miRNAs.csv", "rU")
genes_file = open("/Users/nithintumma/Desktop/genes_ordered.csv", "rU")
paths_file = open("/Users/nithintumma/Desktop/pathways.csv")

pathways = list(csv.reader(paths_file))    # read unique pathways into list
pathways = [ x[0] for x in pathways ]

# --- READ miRNA MASTER LIST (only those in pathways) ---
mirnas = list(csv.reader(mirlist_file))    # read unique miRNAs into list
mirnas = [ x[0] for x in mirnas ]

for path in pathways:
	