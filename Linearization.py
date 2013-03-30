'''
Created on Jan 10, 2012

@author: nithintumma
'''

import xml.etree.ElementTree as ET
import csv

def getRelation (paths, relations):
    name_rel = []
    dict_rel = {}
    for i in paths:
        if len(i) > 1:
            for j in range(len(i) - 1):
                name_rel.append(relations[( i[j], i[j+1] )])
            dict_rel[tuple(i)] = name_rel
            name_rel = []
        else:
            #name_rel.append("none")
            dict_rel[tuple(i)] = "none"
    return dict_rel

def feedbackLoop (path):
    for i in range(len(path)):
        if path[i] in path[i + 1 : len(path)]: return True
        else: continue
    return False

def AddPath (path, build_relation):
    extend_paths = []
    last_gene = path[len(path) - 1]
    connections = []
    connections.extend(build_relation[last_gene])
    
    if len(connections) == 0 : return[path]
    
    if feedbackLoop(path): 
        #print "Feedback"
        return[path]
    
    temp = []
    
    while (len(connections) > 0):
        temp = []
        temp.extend(path)
        temp.append(connections[0])
        extend_paths.append(temp)
        del connections[0]
    return extend_paths


def CurrPathInhibitors (path):
    
    #Parse the xml file for hte current path
    
    curr_path = path
    file_name = "/Users/nithintumma/Desktop/KEGG/" + curr_path +".xml"
    relations = {}
    
    try:
        tree = ET.parse(file_name)
        elem = tree.getroot()
        relations_dict = elem.getiterator("relation")
        entry_dict = elem.getiterator("entry")
    except:
        print ("Error with path: " + path)
        errors.append(path)
        return relations
    else: pass  
    
    # Figure out what genes interact with each other, and dump this into a dictionary of relations
    
    for i in relations_dict:
        if len(i) > 0:
            relations[(i.get("entry1"), i.get("entry2"))] = i[0].get("name")
        else: continue
    
    
    ent1 =[]
    ent2 = []
    
    for entry1, entry2 in relations.keys():
        ent1.append(entry1)
        ent2.append(entry2)
    
    entries = {}
    for i in entry_dict:
        if i.get("type") == "gene":
            entries[i.get("id")] = i[0].get("name")
    
    lone_genes = []
    
    start_path = []
    
    for i in entries.keys():
        if i in ent1 and i not in ent2:
            start_path.append([i])
        elif i not in ent1 and i not in ent2:
            lone_genes.append([i])
    
    build_relation = {}
    
    temp_list=[]
    
    for i in entries.keys():
        for first, second in relations.keys():
            if first == i and second in entries.keys(): temp_list.append(second)
        build_relation[i] = temp_list
        temp_list = []

    #Run the program
    old_temp = []
    temp_index = []
    paths = []
    
    for i in start_path:
        paths.append(i)
    
    #print paths
    
    while (paths != old_temp):
        old_temp = []
        old_temp.extend(temp_index)
        temp_index = []
        for i in paths:
            c = AddPath(i, build_relation)
            for j in c: temp_index.append(j)
        paths = []
        paths.extend(temp_index)
        #print paths
        
        
    for i in lone_genes:
        paths.append(i)
    
    dict_rel = getRelation (paths, relations)
    
    
    inhibit_list = []
    cnt_inhibit = 0
    
    for path in dict_rel.keys():
        for i in range(len(path)):
            for j in range(i, len(path) - 1):
    #            print dict_rel[path][j]
                #print j
                #print dict_rel[path]
                if dict_rel[path][j] == "inhibition": cnt_inhibit+=1
            if cnt_inhibit == 0: continue
            elif cnt_inhibit % 2 == 0: 
                inhibit_list.append([path[i], 0])
            elif cnt_inhibit % 2 == 1:
                inhibit_list.append([path[i], 1])
                                     
            cnt_inhibit = 0
    
    #print inhibit_list
    #print len(inhibit_list)
    sum_inh = 0
    sum_tot = 0
    
    gene_prop = {}
    
    #print entries.keys()
    #print len(entries.keys())
    
    for gene_id in entries.keys():
        for gene, cnt in inhibit_list:
                if gene == gene_id:
                    sum_inh = sum_inh + cnt
                    sum_tot += 1
                    
        gene_prop[gene_id] = [sum_inh, sum_tot]
        sum_inh = 0
        sum_tot=0
    
    num_path = 0
    for gene_id in entries.keys():
        for path in paths:
                if gene_id in path: num_path += 1
        gene_prop[gene_id].append(num_path)
        num_path = 0
        
    inhibit_dict = {}
    
    # Process gene aliases
    for i in gene_prop.keys():
        if gene_prop[i][0] > 0:
            gene_names = entries[i]
            gene_names = gene_names.replace(".", "")
            gene_names = gene_names.replace(",", "")
            gene_names = gene_names.upper()
            gene_names = gene_names.split()
            for gene in gene_names:
                inhibit_dict[gene] = float(gene_prop[i][0])/float(gene_prop[i][2])
    
    #Check if there are any feedback loops
    #feedback_loops = []
    #for path in paths:
    #    for i in range(len(path)):
            #print entries[path[i]]
            #for j in path[i+1:len(path)]: print entries[j] 
    #        if path[i] in path[i + 1 : len(path)]: 
    #            feedback_loops.append[path]
    
    return inhibit_dict

errors = []
#paths_file = open("/Users/nithintumma/Desktop/pathways.csv")
#pathways = []           # List of all of the pathways

#Reader = csv.reader(paths_file, delimiter = "/")
#for line in Reader:
#    pathways.append(line[0])

print "Begin"
errors = []
mirlist_file = open("/Users/nithintumma/Desktop/miRNAs.csv", "rU")
mirna_file = open("/Users/nithintumma/Desktop/miRNA_ordered.csv", "rU")
scores_file = open("/Users/nithintumma/Desktop/score_ordered.csv", "rU")
pathways_file = open("/Users/nithintumma/Desktop/Pathways_ordered.csv", "rU")
genes_file = open("/Users/nithintumma/Desktop/genes_ordered.csv", "rU")
paths_file = open("/Users/nithintumma/Desktop/pathways.csv")

pathways = list(csv.reader(paths_file))    # read unique pathways into list
pathways = [ x[0] for x in pathways ]

# --- READ miRNA MASTER LIST (only those in pathways) ---
mirnas = list(csv.reader(mirlist_file))    # read unique miRNAs into list
mirnas = [ x[0] for x in mirnas ]

# --- READ RELATIONS FILE INTO LIST ---
# "Pathway","Gene ID","Gene Name","Copies in Path","Transcript","miRNA","Score","Start","End"
#  0         1         2           3                4            5       6       7       8
relation = list(csv.reader(open("/Users/nithintumma/Desktop/relation.csv", "rU")))  # get relations
relation = relation[1:]                                 # discard header

# --- BUILD EMPTY DICTIONARY OF miRNAs and Pathways FROM RELATIONS (there may be some we're not interested in, i.e., not in pathways) ---
allRelMirs = set([])
allRelPaths = set([])
for r in relation:
    #print r
    allRelMirs.add(r[5])
    allRelPaths.add(r[0])
print 'Unique miRNAs in relation file = ',len(allRelMirs)
print 'Unique Pathways in relation file = ',len(allRelPaths)
mir_scores  = dict([ (x,[]) for x in list(allRelMirs) ])        # make dictionary from list of miRNAs (one list for each mir)
mir_targets = dict([ (x,[]) for x in list(allRelMirs) ])
path_genes = dict([ (x,set([])) for x in list(allRelPaths) ])   # make dictionary from list of Paths (one set for each path)
path_mirs  = dict([ (x,set([])) for x in list(allRelPaths) ])

# --- POPULATE DICTIONARIES ---
for r in relation:                                              # go through relations once and populate all dictionaries
    mir_scores[ r[5] ].append( [r[0],r[6]] )
    tmp_name = r[0] + " " + r[2]
    mir_targets[ r[5] ].append( [tmp_name, r[6]] )
    path_genes[ r[0] ].add( r[2] )
    path_mirs[ r[0] ].add( r[5] )
for key in path_genes.keys():                                   # convert sets to lists (just to be nice)
    path_genes[key] = list(path_genes[key])
    path_mirs[key] = list(path_mirs[key])


print("Initialization Complete")




print "Genes"
b = path_genes["hsa05218"]
print path_genes["hsa05218"]

print "Inhibitors"
a = CurrPathInhibitors("hsa05218")
print a

print "Ind. Inhibitors"
for i in a:
    print i
    for j in b:
        if i == j: print "In pathway"

print "mir-126 Targets in Path"
for gene_name, score in mir_targets["hsa-miR-126"]:
    gene_name = gene_name.split()
    if gene_name[0] == "hsa05218": 
        print (gene_name[1] + "," + score)

print "mir-Targets that are inhibitors"
for gene_name, score in mir_targets["hsa-miR-126"]:
    gene_name = gene_name.split()
    if gene_name[0] == "hsa05218": 
        if gene_name[1] in a: print (gene_name[1] + "," + score)

    



#cit_cycle = {'IDH2': 13.758100000000001, 'MDH1': 13.644399999999999, 'IDH3G': 14.394500000000001}
#lung_cancer = {'CDKN2A': 15.351800000000001, 'HRAS': 14.200699999999999, 'FOXO3': 14.2614, 'BAD': 15.3956, 'PIK3R2': 15.176600000000001, 'PRKCB1': 14.2257}

#cnt = 0
#for gene in lung_cancer:
#    if gene in CurrPathInhibitors("hsa05223"):
#        print gene
#        print lung_cancer[gene]
#        cnt +=1
        
#if cnt == 0:
#    print "None"
#cnt = 0
#for path in pathways:
#    print path
#    a = CurrPathInhibitors(path)
#    print a
#    if len(a) > 0: cnt +=1

#print cnt
#print len(pathways)
#print errors

print "Done"
