'''
Created on Dec 27, 2011
Things to do: CHANGE GENE NAME FROM NAME TO ENTREZ GENE ID LOCATED IN RELATIONSHIPS FILE AND XML
@author: nithintumma
'''
import random
import numpy
import csv
from Inhibitor_Modules import *
import os
#import Graph as g
import math

def getMirsPath (path):
    tmp_mirs = []
    tmp_mirs.extend(path_mirs[path])
    mir_list = []
    for i in tmp_mirs:
        if i in mirnas:
            mir_list.append(i)
    return mir_list

def getMirs (num):    # This function takes in a number of genes to return as a paramater. Using the list of mir's, it randomly creates a list of mirs and returns this mirlist
    mir_gene = []
    k = []
    i_list = []
    for i in range(num):
        k = random.sample(range(len(mir_list)), num)
        i_list.append(i)
    for j in k:
        mir_gene.append(mir_list[j])
    return mir_gene

def Legalize (pop):            # During initial testing, at some times the population would treat mir_genes of length one as strings instead of lists - which caused errors throughout the rest of the program. This function ensures that each member of the population is packaged as a list.
    for i in range(len(pop)):
        if (isinstance(pop[i], list)): 
            pass
        else:
            pop.append([pop[i]])
            del pop[i]    
    return pop

def Leg_Length_Mir (mir_gene, limit = 4):
    tmp_group = []
    tmp = []
    if len(mir_gene) <= limit: return mir_gene        
    else: 
        divisions = []
        for i in range(20):
            divisions.append(random.sample(range(len(mir_gene)), limit))
        for divide in divisions:
            for i in divide:
                tmp.append(mir_gene[i])
            tmp_group.append(tmp)
            tmp = []
        tmp_sort = Sort_Fitness(tmp_group)
        return tmp_sort[0]
                  
def Sort_Fitness (pop): # Returns the population sorted by fitness
    sort_pop = []
    pop_fit = ScoreFunction (pop, curr_path, mir_scores, inhibitors)
    for mir_gene, score in sorted(zip(pop, pop_fit), key = lambda x:x[1], reverse = 1):
        sort_pop.append(mir_gene)
    return sort_pop
     
def ScoreFunction ( pop , path , mir_scores, inhibitors ):  # Takes in the population, pathway, and the mir_scores dictionary as inputs. Uses the Score function to build the fitness list. 
    curr_fitness = []
    for i in range(len(pop)):
        curr_fitness.append(Score_Covering_No_Inh(pop[i]))  #Need to Change to selectable type, for Score - need to add (path, mir_scores) as paramaters
    return curr_fitness

def Score (mir_gene, path, mir_scores, inhbitors):  # This function takes in a mir_gene and the mir_sore dictionary, as well as the current path as inputs. The mir_score dictionary returns a list of the form [[path_1, score_1], [path_2, score_2], ..., [path_i, score_i]]. This dictionary is unpacked, and the scores are added to create a total_sum for each
                                            #the mir_gene, as well as the pathway_score - the sum of the scores of the mir's in the gene that target the pathway. THe formula for score is computed as score_pathway/len(mir_gene)(score_total - score_pathway). This fitness function favors selective mir_genes, as well as provides an in-build mechanism 
                                            # to prevent indiscriminate growth of the gene. It may be necessary to modify this equation to produce optimal results.    
    if (len(mir_gene) == 0): return 0.0
    
    gene_score = 0.0                            
    score = 0.0
    
    for mir in mir_gene:
        sum_total = 0.0
        sum_pathway = 0.0
        
        for pathway, score in mir_scores[mir]:
            sum_total =  sum_total + score
            if pathway == path: sum_pathway = sum_pathway + score
        
        if (sum_total == sum_pathway):
            gene_score = gene_score + sum_pathway/(sum_total)
        else:
            gene_score = gene_score + sum_pathway/(sum_total - sum_pathway)    
    
    score = 100000.0*gene_score/float(len(mir_gene))
    
    return score    

def Score_Covering (mir_gene, inhbitors):
    # --- Scoring based on the inhibitor value of the gene (coeff). If gene is an inhibitor, subtract/add coeff depending upon value of coeff. Else add 1. Need to ensure genes are distinct, which is what set operator is used for -----
    #---Ensure no length mir_combos are scored zero---
    if len(mir_gene) == 0: return 0.0 
    
    end_score = 0.0
    tmp_genes = []
    for mir in mir_gene:
        for gene_name, score in mir_targets[mir]:
            gene_name = gene_name.split()
            if gene_name[0] == curr_path:
                tmp_genes.append(gene_name[1])
    tmp_genes = list(set(tmp_genes))
    # ---- Normalize score to the length of the curr_path
    for gene in tmp_genes:
        if gene in inhibitors:
            if inhibitors[gene] > 0.5:
                end_score -= inhibitors[gene]
            else: end_score += inhibitors[gene]
        else:
            end_score += 1
    end_score = float(end_score)/float(len(path_genes[curr_path]))
    return end_score

def Score_Covering_No_Inh (mir_gene):
    # --- Ensure that the no-length paths are given a score of zero
    if len(mir_gene) == 0: return 0.0 
    
    end_score = 0.0
    tmp_genes = []
    for mir in mir_gene:
        for gene_name, score in mir_targets[mir]:
            gene_name = gene_name.split()
            if gene_name[0] == curr_path:
                tmp_genes.append(gene_name[1])
    tmp_genes = list(set(tmp_genes))
    # ------- Normalize score to length of the current path ----
    score = float(len(tmp_genes))/float(len(path_genes[curr_path]))
    # --- Select against long combinations of miRNA's
    # ----- Normalizing to the eight root of the length of the combination appears to select for the 3-5 mirs----
    return score/float(math.sqrt((len(mir_gene))))
    
def Mutate (pop, mut_rate):  # Here, we iterate through the population, and create a list of equal length as the current mir_gene, and populate it with random values in the interval [0, 1]. If the value of the list at index i (list[i]) < the mut_rate - the mir at mir_gene[i] is swapped with another mir from the mir_list.
    len_pop = len(pop)
    
    for i in range(len_pop):
        
        pop = Legalize(pop)
        cur_gene = pop[i]
        cur_len = len(cur_gene)
        mut_list = []
        temp_index = []

        for j in range(cur_len):
            mut_list.append(random.random())
        
        for j in range(len(mut_list)): 
            if mut_list[ j ] < mut_rate:
                temp_index.append(cur_gene[j])       
        
        for j in temp_index:
            cur_gene.remove(j)
            cur_gene.append( getMirs( 1 )[0])

        pop[i] = cur_gene        
    
    pop = Legalize(pop)
    return pop

def Insertion (pop, insert_rate):  #Utilizes the same mechanics as the Mutate function - although instead of swapping a gene, if a mir_gene is selected for Insertion, a random mir is appended to the mir_gene.

    len_pop = len(pop)

    for i in range(len_pop):
        cur_gene_x = pop[i]        
        cur_len = len(cur_gene_x)
        insert_list = []        
        temp_index = []

        for j in range(cur_len):
            insert_list.append(random.random())
        
        for j in range( len( insert_list ) ): 
            if insert_list[ j ] < insert_rate:
                temp_index.append(j)    
        
        for j in temp_index:
            cur_gene_x.append( getMirs( 1 )[0])

        pop[i] = cur_gene_x        
    
    #pop = Legalize(pop)
    return pop

def Deletion (pop, del_rate): #Utilizes the same mechanics as the Mutate function - although instead of swapping a gene, if a mir in a mir_gene is selected for Insertion, a mir is deleted from the list. If the length of the mir_gene is one, the function passes - ensuring that blank mir_genes do not accumulate in the population.
    len_pop = len(pop)
    
    for i in range(len_pop):
        cur_gene = pop[i]    
        cur_len = len(cur_gene)
   
        if len(cur_gene) == 1: continue
        
        del_list = []
        temp_index = []

        for j in range(cur_len):
            del_list.append(random.random())
        
        for j in range( len( del_list ) ): 
            if del_list[ j ] < del_rate:
                temp_index.append(cur_gene[j])    
        
        for j in temp_index:
            cur_gene.remove(j)

        pop[i] = cur_gene        
    
    pop = Legalize(pop)
    return pop

def StopCriteria (fitness, history, stall_gens):  #Used to determine the stopping point of the GA. Keeps a history list - containing the average fitness of the last 40 runs of the GA. If the st deviation of history drops below a set value (limit), the GA is stopped.
    
    sum_fit =  numpy.sum(fitness)
    
    avg = sum_fit/len(fitness)

    
    stop_condition = 0        
    
    limit = 0.0  

    if len(history) < stall_gens:
        history.append(avg)
    else:
        del history[0]
        history.append(avg)
    
    if (numpy.std(history) < limit):
        stop_condition = 1

    return[history, stop_condition]
    
def Mate (pop, fitness, curr_path, mir_scores): #This is the most complex function; here we iterate through the population len(pop)/2 times, and randomly select two parents at each iteration to Swap - calling the Swap method to determine the kids. Each mir_gene is assigned a probability of being selected, proportional to the mir_gene's fitness.
                                               # This is achieved through the sum_rand/r_dad loop. After the parents are selected, they are passed to the Swap function, adn the kids from the Swap function are added to a list of kids. After the loop iterates through the population, the kids are added to the population and the parents are deleted. Because it is likely
                                                # that certain parents are chosen multiple times, the program counts the number of parents that are deleted, and only adds that number of kids to the population to ensure that the population does not continue to grow in size.
                                                
                                                # Problems - It might be beneficial to define a max_length of the mir_genes. It might also be beneficial to determine if the kids are actually more fit than the parents, and choose which subset of parents/kids to add to the population and which ot delete. 
    pop_len = len(pop)
    new_pop = []
    
    
    prob = fitness
    prob = prob/numpy.sum(prob)  #Ensure that the sum of prob is = 1
    
    #print ("Prob")
    #print prob
    
    kids = []
    curr_kids = []
    
    parents = []
    
    Num_Genes = len(pop)
    
    #------We are mating 25% of the original population - yielding len(pop)/2 kids ------
    for i in range(Num_Genes/2): #temp change to 50%
        rdad = random.random()    #Random value between 0 and 1
        rmom = random.random()
        
        sum_rand = 0.0
        idad = 0
        while rdad > sum_rand:
            sum_rand += prob[idad]
            idad = idad + 1
        
        
        sum_rand = 0.0
        imom = 0
        while rmom > sum_rand:
            sum_rand += prob[imom]
            imom = imom + 1

        idad = idad - 1
        imom = imom - 1
        dad = pop[idad]
        mom = pop[imom]
        
        # ------ Mom and Dad need to be lists for the rest of the program to proceed -----
        if (isinstance(mom, list)): pass
        else:
            mom = [mom]
        if (isinstance(dad, list)): pass
        else:
            dad = [dad]
        
        #-----Swap the parents based on the fitness----
        curr_kids = Swap( dad, mom, curr_path, mir_scores)
        
        #-----Ensure that Swap function returned the proper values------
        if (len(curr_kids) < 2):                     
            print "Error: Current kids too short"
            continue
        
        kids.append(curr_kids[0])
        kids.append(curr_kids[1])
        
        parents.append(dad)
        parents.append(mom)
    
    # ----- Elitism with 10% of the fittest individuals from the original population carrying over
    elite_prop = 10 #---temp make the elitism 10
    sort_pop = Sort_Fitness (pop)
    for i in range(len(sort_pop)/elite_prop):
        new_pop.append(sort_pop[i])
   
    # ----- Add all of the kids to the new population
    #for kid in kids:
    #    new_pop.append(kid)
    
    # ----- Add the top scoring 50% of kids to population
    sort_kids = Sort_Fitness(kids)
    for index in range(len(sort_kids)/2):
        new_pop.append(sort_kids[index])
        
    # --- Roulette style selection of the remainder of the members of the new population
    num =  pop_len- len(new_pop)    #This is the number that still need to be carried over to maintain an equivalent population size
    
    rest_pop = pop[len(sort_pop)/elite_prop + 1: len(pop)]
    rest_fitness = ScoreFunction ( rest_pop , curr_path , mir_scores, inhibitors )
    rest_prob = rest_fitness/numpy.sum(rest_fitness)
    
    for i in range(num):
        rkeep = random.random()    #Random value between 0 and 1
        
        sum_rand = 0.0
        ikeep = 0
        while rkeep > sum_rand:
            sum_rand += rest_prob[ikeep]
            ikeep = ikeep + 1
        
        ikeep -=1
        tmp = pop[ikeep]
        
        if (isinstance(tmp, list)): pass
        else: tmp = [tmp]
        
        new_pop.append(tmp)

    # ---- Ensure that the genes in the population are lists and they do not exceed the specified length
    for mir_gene_index in range(len(new_pop)):
        new_pop[mir_gene_index] = Leg_Length_Mir(new_pop[mir_gene_index])
    new_pop = Legalize(new_pop)
    return new_pop    

def Swap (dad, mom, curr_path, mir_scores): # This function acts as a supplement to the Mate population, and actually mates two mir_genes. 
    
    if (isinstance(mom, list)): pass
    else:
        mom = [mom]
    
    if (isinstance(dad, list)): pass
    else:
        dad = [dad]
        
    if len(dad) == 0: dad = getMirs(random.randint(1, 5))   # Ensure that parents are not empty lists - should not be called.
    if len(mom) == 0: mom = getMirs(random.randint(1, 5))
    
    temp_dad = []
    temp_mom = [] 
    
    del_list = []
    
    if len(dad) == 1:                                # Ensure that if the gene is of length one
        if random.random() > 0.5: num_dad = 0
        else: num_dad = 1
    else:
        num_dad = random.randint(1, len(dad))
    
    if len(mom) == 1:
        if random.random() > 0.5: num_mom = 0
        else: num_mom = 1
    else:
        num_mom = random.randint(1, len(mom))
    
    while len(dad) < num_dad: num_dad-=1
    
    r_dad = random.sample(range(len(dad)), num_dad)  # Select num_dad numbers randomly from the interval [0, len(dad)]
    
    #print r_dad
    #print "dad"
    #print dad
    
    for i in r_dad:
        temp_dad.append(dad[i])
        del_list.append(dad[i])
    
    for i in del_list:
        dad.remove(i)

    del_list = []
    
    while len(mom) < num_mom: num_mom-=1               # Ensure that there is not a problem with random.sample
    
    r_mom = random.sample(range(len(mom)), num_mom)
    
    #print r_mom
    #print "mom"
    #print mom
   
    for i in r_mom:
        temp_mom.append(mom[i])
        del_list.append(mom[i]) 
        
    for i in del_list:
        mom.remove(i)
    
    kid1 = temp_dad + mom
    kid2 =  temp_mom + dad
    
    temp_mom = []
    temp_dad = []
        
    # -----Ensure no duplication of miRNA's in a single gene ---- 
    kid1 = list(set(kid1))
    kid2 = list(set(kid2))
    
    if len(kid1) == 0:
        #print "Error Kid 1"
        kid1 = getMirs(random.randint(1, 5))
    if len(kid2) == 0:
    	#print "Error kid 2"
        kid2 = getMirs(random.randint(1, 5))
            
    return [kid1, kid2]

def Diversity (pop):
    mirs_tot = []
    for mir_gene in pop:
        for mir in mir_gene:
            mirs_tot.append(mir)
    mirs_tot = list(set(mirs_tot))
    # ---- Compute Hamming distance pairwise between the members of the population --- 
    dist = 0
    for i in range(len(pop)):
        for j in range(len(pop) - i):
            if i == i+j: continue
            else: 
                targets_1 = []
                targets_2 = []
                for mir in pop[i+j]:
                    for gene_name, score in mir_targets[mir]: targets_1.append(gene_name)
                for mir in pop[i]:
                    for gene_name, score in mir_targets[mir]: targets_2.append(gene_name)
                targets_1 = set(targets_1)
                targets_2 = set(targets_2)
                diff = targets_1.symmetric_difference(targets_2)
                dist +=  len(diff) 
    
    return len(mirs_tot), dist

def Diversity_Num (pop):
    mirs_tot = []
    for mir_gene in pop:
        for mir in mir_gene:
            mirs_tot.append(mir)
    mirs_tot = list(set(mirs_tot))
    return len(mirs_tot)

def Calc_Freq (pop, freq):
	# -----  For each mir in the pop, calculate the number of times it appears in a combo and add it to the freq----
	for mir in mir_list:
		for mir_combo in pop:
			if mir in mir_combo: freq[mir] = freq[mir] + 1
	
	return freq
	

curr_path = "hsa04910"    # Just chose the Tryptophan metabolism, should be able to change this to any pathway.

os.chdir("/Users/nithintumma/Downloads/mirfiles")
mirlist_file = open("mirnas.csv", "rU")

# --- Open requisite files ---    
#genes_file = open("genes.csv", "rU")
paths_file = open("pathways.csv", "rU")

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
    mir_scores[ r[5] ].append( [r[0], float(r[6])] )
    tmp_name = r[0] + " " + r[2]
    mir_targets[ r[5] ].append( [tmp_name, float(r[6])] )
    path_genes[ r[0] ].add( r[2] )
    path_mirs[ r[0] ].add( r[5] )
for key in path_genes.keys():                                   # convert sets to lists (just to be nice)
    path_genes[key] = list(path_genes[key])
    path_mirs[key] = list(path_mirs[key])

# --- mir_targets:= mirna ---> ("path gene", score), mir_scores := mirna ----> ("path", score)
# ---- path_mirs := path ----> mirna, path_genes:= path ---> gene
# pathways:= list of pathways (unique human), mirnas:= list of mirnas (unique human)
# -----Initialize output files ----- 
o_file = open("/Users/nithintumma/Desktop/Output/output.txt", "r+")   # The output file
fit_file = open("/Users/nithintumma/Desktop/Output/fitness.csv", "r+") # The fitness output file
maxfit_file = open("/Users/nithintumma/Desktop/Output/maxfit.csv", "r+") #the Maximum fitness output file 
print("Initialization Complete")

# ---- Declare population, fitness, history lists ----
population = []
fitness = []
history = []
    
# ----- Parameters that should be changeable -----
init_size = 60
max_gene_len = 3
stall_gens = 30

# --- Contains all of the human mirs that target the curr_path ----
mir_list = getMirsPath(curr_path)

#--------------------------------
# 			INITIALIZE			|
#--------------------------------

# ----- Populate initial random population -----
#for i in range(max_gene_len):
#    for j in range(init_size):
#        population.append(getMirs(i+1))

for i in range(300):
    population.append(getMirs(max_gene_len))

frequency = {}
for mir in mir_list:
	frequency[mir] = 0

inhibitors = CurrPathInhibitors(curr_path)

# -------- Find fitness for initial population -----
fitness = ScoreFunction(population, curr_path, mir_scores, inhibitors)

ok = 1
stop_condition = 0
cnt = 0

while ok:    
    population = Legalize(population)
    population = Mate(population, fitness, curr_path, mir_scores)    
    population = Mutate(population, 0.01)
    #population = Insertion(population, 0.005)    
    #population = Deletion(population, 0.005)
    population = Legalize(population)
    for mir_index in range(len(population)):
        population[mir_index] = Leg_Length_Mir(population[mir_index])
    
    sort_pop = Sort_Fitness(population)
    fitness = ScoreFunction(population, curr_path, mir_scores, inhibitors)      #Determine the fitness of the population after Mating
    num_mirs = Diversity_Num(population)
    frequency = Calc_Freq(population, frequency)
    
    # ---- Print output after Generation
    print cnt
    print str(numpy.average(fitness))
    print str(numpy.max(fitness))
    print("Current Best Score = " + str(max(fitness)))
    print("By: " + str(population[fitness.index(max(fitness))]))
    print (num_mirs)
    #print (dist)
    
    maxfit_file.write(str(max(fitness)) + "," )
    fit_file.write(str(numpy.average(fitness)) + ", ")
    

	
#    graph_dict["gen"].append(cnt)
#    graph_dict["avg_fit"].append(numpy.average(fitness))
#    G1.plot(graph_dict)
    
    
    
    o_file.write("Cycle Number: " + str(cnt) + ", " + "\n")
    o_file.write("Population: " + str(sort_pop) + ", "  + "\n")
    o_file.write("Avg. Fitness: " + str(fitness) + ", " + "\n")
    o_file.write("Max Fitness: " + str(max(fitness)) + ", " + "\n" )
    o_file.write(("Average: " + str(numpy.average(fitness))) + ", " + "\n")
    o_file.write("Current Best Score (After Mutation): " + str(max(fitness)) + ", " + "\n")
    o_file.write("By: " + str(population[fitness.index(max(fitness))]) + ", " + "\n")
    
    history, stop_condition = StopCriteria(fitness, history, stall_gens)

    if cnt == 1200:      #Currently, just running for 5,000 cycle, will eventually use stop_condition
        ok = 0
    cnt = cnt + 1
    
#print("Best Gene After " + str(cnt) + " = " + str(population[fitness.index(max(fitness))]))
o_file.write("Best Gene After: " + str(cnt) + " = " + str(population[fitness.index(max(fitness))]) + ", " + "\n")        
o_file.close()
print("Finished")

