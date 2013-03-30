'''
Created on Jan 10, 2012

@author: nithintumma
'''
# ---- NEED TO MAKE CHANGE: GENE NAMES TO ENTREZ ID FROM THE XML FILE

import xml.etree.ElementTree as ET

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
    
    #Parse the xml file for the current path
    errors = []
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








