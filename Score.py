'''
Created on Jan 10, 2012

@author: nithintumma
'''

import csv
import xml.etree.ElementTree as ET
import numpy
import os

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

def Score (mir, path, inhibitor_dict):
    
    end_score = 0.0
  
    for path_gene, gene_score in mir_targets[mir]:
        path_gene = path_gene.split()
        if path_gene[0] == path:
            if path_gene[1] in path_genes[path]: 
                if path_gene[1] in inhibitor_dict:
                    coeff = inhibitor_dict[path_gene[1]]
                    #print coeff
                    # If it inhibits more linear sub-paths than it activates
                    if coeff > 0.5: end_score += float(gene_score) * float(coeff)
                    else: end_score -= float(gene_score) * float(coeff)
                # If it is not an inhibitor of any linear pathways
                else: end_score -= gene_score
                    
    
    end_score = end_score/float(len(path_genes[curr_path]))
    return end_score


def CurrPathInhibitors (curr_path):
    
    #Parse the xml file for hte current path
    relations = {}
    file_name = "/home/nithin/Desktop/Data/KEGG/" + curr_path +".xml"
    
    try:
        tree = ET.parse(file_name)
        elem = tree.getroot()
        relations_dict = elem.getiterator("relation")
        entry_dict = elem.getiterator("entry")
    except:
        print ("Error with path: " + curr_path)
        errors.append(curr_path)
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
    
    
    while (paths != old_temp):
        old_temp = []
        old_temp.extend(temp_index)
        temp_index = []
        for i in paths:
            c = AddPath(i, build_relation)
            for j in c: temp_index.append(j)
        paths = []
        paths.extend(temp_index)
        
    # add single genes with no relationships    
    for i in lone_genes:
        paths.append(i)
    
    dict_rel = getRelation (paths, relations)
    
    
    inhibit_list = []
    cnt_inhibit = 0
    
    for path in dict_rel.keys():
        for i in range(len(path)):
            for j in range(i, len(path) - 1):
                if dict_rel[path][j] == "inhibition": cnt_inhibit+=1
            if cnt_inhibit == 0: continue
            elif cnt_inhibit % 2 == 0: 
                inhibit_list.append([path[i], 0])
            elif cnt_inhibit % 2 == 1:
                inhibit_list.append([path[i], 1])
                                     
            cnt_inhibit = 0
    

    sum_inh = 0
    sum_tot = 0
    
    gene_prop = {}
    
 
    
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
    #        if path[i] in path[i + 1 : len(path)]: 
    #            feedback_loops.append[path]
    
    return inhibit_dict

def writeCSV(score_dict, out_file = "output.csv"):
    writer = csv.writer(open(out_file, 'wb'), delimiter=',')
    
    path_description = {}
    for path in pathways:
        curr_path = path
        file_name = "/home/nithin/Desktop/Data/KEGG/" + curr_path +".xml"
        tree = ET.parse(file_name)
        elem = tree.getroot()
        path_description[curr_path] = elem.get("title")

    tmp = [" "]
    for path in pathways:
        tmp.append(path_description[path])
    
    pathways.insert(0, " ")
    
    writer.writerow(pathways)
    writer.writerow(tmp)
    tmp = []
    
    del pathways[0]
    # Need a try catch when writing for mirs/paths that do not interact - add a 0.0
    #Change how we srite the pathways
    tmp_list = []
    for mir in mirnas:
        tmp_list.append(mir)
        for path in pathways:
            try:
                tmp_list.append(score_dict[(mir, path)])
            except:
                tmp_list.append(0.0)
            else: pass
        writer.writerow(tmp_list)
        tmp_list = []

def Enrichment(score_dict):
    new_score_dict = {}
    
    for mir in mirnas:
        avg = 0.0
        all_scores = []
        for path in pathways:
            all_scores.append(score_dict[mir, path])
        avg = numpy.average(all_scores)
        for path in pathways:
            if avg == 0.0: 
                print "Error divide by 0"
                new_score_dict[(mir, path)] = 0.0
            else: 
                new_score_dict[(mir, path)] = score_dict[mir, path]/avg
    
    print new_score_dict
    writeCSV(new_score_dict, "output_Enrichment")

# ------ getMirsPath is not called currently (used in the old code)-------
def getMirsPath (path):
    temp_mirs = []
    for j in range(len(pathways_ordered)):
        if pathways_ordered[j] == path: temp_mirs.append(mirna_ordered[j])
    temp_mirs = list(set(temp_mirs))
    mir_list = []
    for mir in temp_mirs:
        if mir in mirnas:
            mir_list.append(mir)
    return mir_list

print "Begin"
errors = []
os.chdir("/home/nithin/Desktop/Data")
errors = []
mirlist_file = open("miRNAs.csv", "rU")
paths_file = open("pathways.csv", "rU")
genes_file = open("genes.csv", "rU")
pathways = list(csv.reader(paths_file))    # read unique pathways into list
pathways = [ x[0] for x in pathways ]

# --- READ miRNA MASTER LIST (only those in pathways) ---
mirnas = list(csv.reader(mirlist_file))    # read unique miRNAs into list
mirnas = [ x[0] for x in mirnas ]

# --- READ RELATIONS FILE INTO LIST ---
# "Pathway","Gene ID","Gene Name","Copies in Path","Transcript","miRNA","Score","Start","End"
#  0         1         2           3                4            5       6       7       8
relation = list(csv.reader(open("relation.csv", "rU")))  # get relations
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
    mir_scores[ r[5] ].append( [r[0],float(r[6])] )
    tmp_name = r[0] + " " + r[2]
    mir_targets[ r[5] ].append( [tmp_name, float(r[6])] )
    path_genes[ r[0] ].add( r[2] )
    path_mirs[ r[0] ].add( r[5] )
for key in path_genes.keys():                                   # convert sets to lists (just to be nice)
    path_genes[key] = list(path_genes[key])
    path_mirs[key] = list(path_mirs[key])


print("Initialization Complete")
    
# Remember that mir_scores and mir_targets are different. mir_scores maps each mir to the pathways with scores, mir_targets maps each mir to the genes that it targets and the score

print("Initialization Complete")

score_dict = {} # Dictionary of the form (mir, path) = score

#Need to loop through all of the pathways
for path in pathways:
    curr_path = path
    print curr_path
    mir_list = path_mirs[curr_path]
    curr_path_inhibitors = CurrPathInhibitors(curr_path)
    
    for mir in mirnas:
        if mir in mir_list:
            score_dict[(mir, curr_path)] = Score(mir, curr_path, curr_path_inhibitors)
        else: score_dict[(mir, curr_path)] = 0.0 

# Write Data to file    
print score_dict 
writeCSV (score_dict)
#Enrichment (score_dict)
print (errors)
print ("Done")

