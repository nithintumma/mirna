'''
Created on Jan 8, 2012

@author: nithintumma
'''

import xml.etree.ElementTree as ET


def getRelation (paths):
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

def AddPath (path):
    extend_paths = []
    last_gene = path[len(path) - 1]
    #print last_gene
    connections = []
    connections.extend(build_relation[last_gene])
    
    if len(connections) == 0 : return[path]
    
    if feedbackLoop(path): 
        print "Feedback"
        return[path]
    
    temp = []
    
    while (len(connections) > 0):
        temp = []
        temp.extend(path)
        temp.append(connections[0])
        extend_paths.append(temp)
        del connections[0]
    return extend_paths

#Iniitalize the variables

def CheckforFeedback(paths):
    feedbacks = []
    for path in paths:
        for i in range(len(path)):
            if path[i] in path[i + 1 : len(path)]: 
                feedbacks.append(path)
    return feedbacks

curr_path = "hsa05221"
file_name = "/Users/nithintumma/Desktop/KEGG/" + curr_path +".xml"
tree = ET.parse(file_name)
elem = tree.getroot()
relations_dict = elem.getiterator("relation")
entry_dict = elem.getiterator("entry")
relations = {}
print relations_dict

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

gene_entries = {}

for i in entry_dict:
    entries[i.get("id")] = i[0].get("name")

for i in entry_dict:
    if i.get("type") == "gene":
        gene_entries[i.get("id")] = i[0].get("name")
        


    
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
        

#print build_relation
print entries
 


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
        c = AddPath(i)
        for j in c: temp_index.append(j)
    paths = []
    paths.extend(temp_index)
    print paths
    
    
for i in lone_genes:
    paths.append(i)

print "Done"
dict_rel = getRelation (paths)

#for i in dict_rel:
#    print i 
#    print dict_rel[i]
#    for j in i:
#        print entries[j]

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

print gene_prop    
    
#print gene_prop

#def Draw (paths, path_rel_dict):

#print gene_prop["34"]

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

print inhibit_dict
#Check if there are any feedback loops
#for i in paths:////
#    print len(i)






    
