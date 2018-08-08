#!/usr/bin/env python3

##################################################################
# This script calculates different ethnic origin haplotypes      #
# proportions within IBD clusters. Then, it finds polyethinc and #
# monoethnic IBD clusters enriched with affected haplotypes.     #
# This script requires python-3.6.3                              #
# To start type: ./ibdmig.py 22 ibdmig.list 9 6                  #
# Read more in ibdmig.readme                                     #
# By Ural Yunusbaev 2018                                         #
##################################################################

# module load python-3.6.3

import sys
import csv
import os
import random
from decimal import *
from datetime import datetime

file_hcl  = "Missing file_hcl"
file_list = "Missing file_list"
poly_threshold = "Missing polyethnic cluster threshold"
mono_threshold = "Missing monoethnic cluster threshold"
file_map = "Missing file_map"


if   sys.argv[1:]: file_hcl  = sys.argv[1]
if   sys.argv[2:]: file_list = sys.argv[2]
if   sys.argv[3:]: poly_threshold  = sys.argv[3]
if   sys.argv[4:]: mono_threshold  = sys.argv[4]
if   sys.argv[5:]: file_map  = sys.argv[5]
if file_hcl == "Missing file_hcl" or file_list == "Missing file_list" or poly_threshold == "Missing polyethnic cluster threshold" or mono_threshold == "Missing monoethnic cluster threshold":
	print ("Missing command line argument.")
	sys.exit("To start type: ./ibdmig.py 22 ibdmig.list 9 6\nwhere:\n22 - the number of chromosomes in acording to DUSH output files (clust_1.hcl ... clust_22.hcl)\nibdmig.list - the file containing a list of individuals with following columns: ind_id, population, phenotype\n9  - the size threshold for affected polyethnic cluster\n6  - the size threshold for affected monoethnic cluster\nRead more in ibdmig.readme")

if file_map == "Missing file_map":
	print ("Missing file_map argument in command line. Genetic distances in output will be 0.")

work_dir = os.getcwd()
log_file = 'Program           :  ibdmig.py\nAuthor            :  Copyright (C) 2018 Ural Yunusbaev\nStart Time        :  ' + datetime.now().strftime('%Y-%m-%d %H:%M:%S') + '\nWork Directory    :  ' + work_dir + '\n\nParameters\n  chromosomes     :  ' + file_hcl       + '\n' + '  individuals list:  ' + file_list      + '\n' + '  poly  threshold :  ' + poly_threshold + '\n' + '  mono  threshold :  ' + mono_threshold + '\n' + '  genetic distances loaded from :  ' + file_map + '\n\n'
with open("ibdmig.out.log", "w") as f: f.write(log_file)
print (log_file)

if int(poly_threshold) < 6: sys.exit("The poly_threshold can not be less than 6!")
if int(mono_threshold) < 5: sys.exit("The mono_threshold can not be less than 5!")


#####################################################################
#                                                                   #
# 1.Reading the list of individuals with population & phenotype data#
#                                                                   #
#####################################################################

IndPopList = []
with open(file_list) as file:
        for row in file:
                line = (row.strip().split('\t'))
                IndPopList.append(line)
IndPopListLen = len(IndPopList)

log_file =  (str(len(IndPopList)) + " individuals loaded from " + file_list + "\n")

with open("ibdmig.out.log", "a") as f: f.write(log_file)

UniqPopList=[]
for i in range(len(IndPopList)): UniqPopList.append(IndPopList[i][1])
UniqPopList=set(UniqPopList)
UniqPopList = sorted(list(UniqPopList))

log_file = (str(len(UniqPopList)) + " populations detected in " + file_list + ":\n")
UniqPopListNumber = []
for k in  range(len(UniqPopList)): UniqPopListNumber.append([list(UniqPopList)[k], 0])

for i in range(len(IndPopList)):
	for k in  range(len(UniqPopList)):
		if IndPopList[i][1] == list(UniqPopList)[k]:
			UniqPopListNumber[k][1] += 1
#print (UniqPopListNumber)
n=0
for i in range(len(UniqPopListNumber)): 
	n=n+1
	log_file =  log_file + (" " + str(n) + ". " + UniqPopListNumber[i][0]  + " = " + str(UniqPopListNumber[i][1]) + " individuals\n")
with open("ibdmig.out.log", "a") as f: f.write(log_file)
print (log_file)

if len(UniqPopList) > 9: sys.exit("The number of populations can not be more than 9!")


#####################################################################
#                                                                   #
# 1.1.Reading the genetic dastances from map/bim file               #
#                                                                   #
#####################################################################

GenDists = []

if (os.path.exists(file_map)):
	with open(file_map) as file:
        	for row in file:
                	line = (row.strip().split('\t'))
                	GenDists.append(line)
else:
	print(file_map)


log_file =  (str(len(GenDists)) + " position's genetic dastances loaded from " + file_map + "\n")
print (log_file)
with open("ibdmig.out.log", "a") as f: f.write(log_file)


GenDists1col = []
for i in GenDists:
#	print ( str(i[0]) + ":" + str(i[3]) )
	i1 =  ( str(i[0]) + ":" + str(i[3]) )
	GenDists1col.append(i1)


if (len(GenDists)) == (len(GenDists1col)):
	log_file =  ('Uniq genetic dastances positions check - passed.')
else:
	log_file =  ('Uniq genetic dastances positions check - NOT passed!!!')
print (log_file)
with open("ibdmig.out.log", "a") as f: f.write(log_file)


####################################################################
#                                                                  #
# 2. Reading & analysing the IBD clusters                          #
#                                                                  #
####################################################################

####################################################################
#                                                                  #
# Creating the ibdmig.out.total and ibdmig.out.affected files      #
#                                                                  #
####################################################################


clusters = []
ChromosomesNumber = 0
ClustersNumberTotal = 0
log_file = ("\n")
for ch in range(1, int(file_hcl)+1):
	hcl = ("clust_" + str(ch) + ".hcl")
	chrom = [str(ch)]
	ClustersNumber = 0
	with open(hcl) as file:
		for row in file:
			line = (row.strip().replace('.0', '').replace('.1', '').replace('\t', ' ').split(' '))
			line1 = chrom + line[0:3] + line[3::4]
			clusters.append(line1)
			ClustersNumber += 1
			ClustersNumberTotal += 1
	log_file = log_file + (str(ClustersNumber) + " clusters loaded from " + hcl + "\n")
	ChromosomesNumber += 1


log_file = log_file + ("Totaly " + str(ClustersNumberTotal) + " clusters loaded from " + str(ChromosomesNumber) + " chromosomes. \n")
with open("ibdmig.out.log", "a") as f: f.write(log_file)
print (log_file)


IndPopList0 = [ IndPopList[i][0] for i in range(len(IndPopList)) ]

PopsNumber0 = []
for k in  range(len(UniqPopList)): PopsNumber0.append( 0 )

PopsNumberTotal = PopsNumber0
ClstSizeTotal = 0
OutputTable   = []
OutputTable1  = []
OutputTable2  = []
OutputRowGdist   = []
OutputTableGdist = []

for clst in range(len(clusters)):
	PopsNumber0 = []
	for k in  range(len(UniqPopList)): PopsNumber0.append( '0' )
	PopsNumber = PopsNumber0
	ClstSize = 0
	ClstSizePatients = 0
	OutputRow = []
	for element in range(len(clusters[clst])):
		if element < 4:
			OutputRow.append(clusters[clst][element])
		else:
			elementsPop = IndPopList[IndPopList0.index(clusters[clst][element])][1]
			elementsStatus = int(IndPopList[IndPopList0.index(clusters[clst][element])][2])
			for k in  range(len(UniqPopList)):
				if elementsPop == list(UniqPopList)[k]:
					PopsNumber[k] = str(int(PopsNumber[k]) + 1)
					PopsNumberTotal[k] = str(int(PopsNumberTotal[k]) + 1)
					ClstSize += 1
					ClstSizeTotal += 1
			if elementsStatus == 2: ClstSizePatients += 1
	OutputRow.append(str(ClstSize))
	OutputRow.append(str(ClstSizePatients))
	for x in PopsNumber: 
		OutputRow.append(x)
		if int(x) == ClstSize and int(x) >= int(mono_threshold) and int(x) == ClstSizePatients: OutputTable2.append(OutputRow)
	OutputTable.append(OutputRow)

	OutputRowGdist = list(OutputRow)

	if (os.path.exists(file_map)):
		clust = clusters[clst]
		pos_indx1 = ( GenDists1col.index( str(clust[0]) + ":" + str(clust[2]) ) )
		pos_indx2 = ( GenDists1col.index( str(clust[0]) + ":" + str(clust[3]) ) )
		Length_cM = Decimal(GenDists[pos_indx2][2]) - Decimal(GenDists[pos_indx1][2])
		Start_cM  = Decimal(GenDists[pos_indx1][2])
		End_cM    = Decimal(GenDists[pos_indx2][2])
	else:
		Length_cM = 0.0
		Start_cM  = 0.0
		End_cM    = 0.0

	OutputRowGdist.append( Length_cM )
	OutputRowGdist.append( Start_cM )
	OutputRowGdist.append( End_cM )

	OutputTableGdist.append(OutputRowGdist)

	if ClstSizePatients >= int(poly_threshold) and ClstSizePatients == ClstSize: OutputTable2.append(OutputRow)

TableHeader = ['CHR', 'CLASTER', 'START', 'END', 'SIZE', 'AFFECT']
for k in range(len(UniqPopList)): TableHeader.append(str(list(UniqPopList)[k]))

####################################################################
#                                                                  #
# Saving the ibdmig.out.total file                                 #
#                                                                  #
####################################################################

with open("ibdmig.out.total","w") as f:
	f.write( '\t'.join(TableHeader) )
	f.write( '\n' )
	wr = csv.writer(f, delimiter="\t")
	wr.writerows(OutputTable)

log_file = ("\n" + str(len(OutputTable)) + " clusters analysed and saved to ibdmig.out.total\n")
with open("ibdmig.out.log", "a") as f: f.write(log_file)
print (log_file)

####################################################################
#                                                                  #
# Saving the ibdmig.out.total.gdist file                           #
#                                                                  #
####################################################################

TableHeaderGdist = list(TableHeader)
TableHeaderGdist.append('LENGTH_cM')
TableHeaderGdist.append('START_cM')
TableHeaderGdist.append('END_cM')
#print (TableHeaderGdist)

with open("ibdmig.out.total.gdist","w") as f:
	f.write( '\t'.join(TableHeaderGdist) )
	f.write( '\n' )
	wr = csv.writer(f, delimiter="\t")
	wr.writerows(OutputTableGdist)

log_file = ("\n" + str(len(OutputTableGdist)) + " clusters analysed and saved to ibdmig.out.total.gdist\n")
with open("ibdmig.out.log", "a") as f: f.write(log_file)
print (log_file)


####################################################################
#                                                                  #
# The ibdmig.out.total.end file                                    #
#                                                                  #
####################################################################

MaxClusterSize  = max([int(x[4]) for x in OutputTable])
MinClusterSize  = min([int(x[4]) for x in OutputTable])
MeanClusterSize = round((sum([int(x[4]) for x in OutputTable]) / len([int(x[4]) for x in OutputTable])), 1)

TotalEnd = []
TotalEndmax = []
TotalEndmin = []
TotalEndmean = []

for col in range(len(TableHeader)):
	if col > 3: 
		TotalEndmax.append(max([int(x[col]) for x in OutputTable]))
		TotalEndmin.append(min([int(x[col]) for x in OutputTable]))
		TotalEndmean.append(round((sum([int(x[col]) for x in OutputTable]) / len([int(x[col]) for x in OutputTable])), 1))
	elif col == 0: 
		TotalEndmax.append('max')
		TotalEndmin.append('min')
		TotalEndmean.append('mean')
	else: 
		TotalEndmax.append('-')
		TotalEndmin.append('-')
		TotalEndmean.append('-')
TotalEnd.append(TotalEndmax)
TotalEnd.append(TotalEndmin)
TotalEnd.append(TotalEndmean)


####################################################################
#                                                                  #
# Saving the ibdmig.out.total.end file                             #
#                                                                  #
####################################################################

with open("ibdmig.out.total.end","w") as f:
	f.write( '\t'.join(TableHeader) )
	f.write( '\n' )
	wr = csv.writer(f, delimiter="\t")
	wr.writerows(TotalEnd)

log_file = ("Common statistics of clusters analysis saved to ibdmig.out.total.end\n")
with open("ibdmig.out.log", "a") as f: f.write(log_file)
print (log_file)




####################################################################
#                                                                  #
# The ibdmig.out.total.gdist.end file                              #
#                                                                  #
####################################################################

MaxClusterSize  = max([int(x[4]) for x in OutputTableGdist])
MinClusterSize  = min([int(x[4]) for x in OutputTableGdist])
MeanClusterSize = round((sum([int(x[4]) for x in OutputTableGdist]) / len([int(x[4]) for x in OutputTableGdist])), 1)

TotalEnd = []
TotalEndmax = []
TotalEndmin = []
TotalEndmean = []

for col in range(len(TableHeaderGdist)):
	if col > 3: 
		TotalEndmax.append(max([int(x[col]) for x in OutputTableGdist]))
		TotalEndmin.append(min([int(x[col]) for x in OutputTableGdist]))
		TotalEndmean.append(round((sum([int(x[col]) for x in OutputTableGdist]) / len([int(x[col]) for x in OutputTableGdist])), 1))
	elif col == 0: 
		TotalEndmax.append('max')
		TotalEndmin.append('min')
		TotalEndmean.append('mean')
	else: 
		TotalEndmax.append('-')
		TotalEndmin.append('-')
		TotalEndmean.append('-')
TotalEnd.append(TotalEndmax)
TotalEnd.append(TotalEndmin)
TotalEnd.append(TotalEndmean)


####################################################################
#                                                                  #
# Saving the ibdmig.out.total.gdist.end file                       #
#                                                                  #
####################################################################

with open("ibdmig.out.total.gdist.end","w") as f:
	f.write( '\t'.join(TableHeaderGdist) )
	f.write( '\n' )
	wr = csv.writer(f, delimiter="\t")
	wr.writerows(TotalEnd)

log_file = ("Common statistics of clusters analysis saved to ibdmig.out.total.gdist.end\n")
with open("ibdmig.out.log", "a") as f: f.write(log_file)
print (log_file)





####################################################################
#                                                                  #
# Saving the ibdmig.out.affected file                              #
#                                                                  #
####################################################################



with open("ibdmig.out.affected","w") as f:
	f.write( '\t'.join(TableHeader) )
	f.write( '\n' )
	wr = csv.writer(f, delimiter="\t")
	wr.writerows(OutputTable2)

log_file = (str(len(OutputTable2)) + " clusters passed the threshold saved to ibdmig.out.affected.\n")
with open("ibdmig.out.log", "a") as f: f.write(log_file)
print (log_file)


####################################################################
#                                                                  #
# Creating the ibdmig.out.proportion file                          #
#                                                                  #
# Making the cluster size_list                                     #
#                                                                  #
####################################################################

OutputTable3=[]
for i in range(len(OutputTable)): OutputTable3.append(int(OutputTable[i][4]))
OutputTable3 = set(OutputTable3)
OutputTable3 = sorted(list(set(OutputTable3)))


####################################################################
#                                                                  #
# Getting the combinations of populations list[PopCombinations]    #
#                                                                  #
####################################################################

from itertools import combinations

UniqPopListSize = len(UniqPopList)

UniqPopListSizeString = ''
UniqPopListSizeString000 = ''

for s in range(UniqPopListSize):
	UniqPopListSizeString = UniqPopListSizeString + str(s+1)

UniqPopListNameString = []
UniqPopListNameString000 = []
PopCombinations = []
PopCombinations1 = []
CombPopNameTable = []
CombPopNameTableId = []

for i in range(len(UniqPopList)): UniqPopListNameString.append('0')

for x in range(UniqPopListSize):
	UniqPopListNameString000 = list(UniqPopListNameString)
	c = list(combinations(UniqPopListSizeString, x+1))
	for x in range(len(c)):
		UniqPopListNameString000 = list(UniqPopListNameString)
		for y in range(len(c[x])):
			UniqPopListNameString000[int(c[x][y])-1] = '1'
		PopCombinations.append(UniqPopListNameString000)

	for x in c:
		CombPopNameRow = []
		CombPopNameRowId = []
		for y in x:
			CombPopNameRowId.append(y)
			CombPopNameRow.append(UniqPopList[int(y)-1])

		CombPopNameTableId.append(CombPopNameRowId)
		CombPopNameTable.append(CombPopNameRow)


for x in PopCombinations:
	PopCombinations1.append(''.join(x))
PopCombinations = list(PopCombinations1)


####################################################################
#                                                                  #
# Getting the ibdmig.out.proportion [PopCombinationsTable]         #
#                                                                  #
####################################################################

PopCombinationsTable = []
PopCombinationsTableInd = []

for x in PopCombinations:
	for y in OutputTable3:
		TableRow = []
		TableRow.append(x+str(y))
		TableRow.append(x)
		TableRow.append(y)
		TableRow.append(0)
		PopCombinationsTable.append(TableRow)
for x in PopCombinationsTable:
	PopCombinationsTableInd.append(x[0])

PopCombinationsTableGdist = list(PopCombinationsTable)
#print ("PopCombinationsTable")
#print (PopCombinationsTable)

####################################################################
#                                                                  #
# Fullihg the ibdmig.out.proportion [PopCombinationsTable]         #
#                                                                  #
####################################################################

#print ("OutputTable")
#print (OutputTable)
#print ("OutputTableGdist")
#print (OutputTableGdist)


for x in OutputTable:
	ClstSize = x[4]
	string1 = ''
	for y in range(len(UniqPopList)):
		string0 = '0'
		if int(x[6+y]) > 0:
			string0 = '1'
		string1 = string1 + string0
	string2 = string1 + ClstSize
	index1 = (PopCombinationsTableInd.index(string2))
	PopCombinationsTable[index1][3] = PopCombinationsTable[index1][3] + 1

#print("PopCombinations")
#print(PopCombinations)
#print("OutputTable3")
#print(OutputTable3)
#print("PopCombinationsTableInd")
#print(PopCombinationsTableInd)
#print("PopCombinationsTable")
#print(PopCombinationsTable)


	
OutputTable40 = []
for x in PopCombinations:
	OutputTable4 = []
	string2total_for_row = 0
	OutputTable4.append(x)
	for y in OutputTable3:
		string2 = str(x) + str(y)
		index1 = (PopCombinationsTableInd.index(string2))
		string2 = PopCombinationsTable[index1][3]
		#print("string2 = " + str(string2))
		string2total_for_row = string2total_for_row + int(string2)
		OutputTable4.append(string2)
	OutputTable4.append(string2total_for_row)
	#print("string2total_for_row")
	#print(string2total_for_row)
	OutputTable40.append(OutputTable4)
	#print("OutputTable4")
	#print(OutputTable4)

#print("OutputTable40")
#print(OutputTable40)

OutputTable40head = ['POPS']
for x in OutputTable3: OutputTable40head.append(str(x))
OutputTable40head.append('TOTAL')


####################################################################
#                                                                  #
# Saving the ibdmig.out.proportion file                            #
#                                                                  #
####################################################################

with open("ibdmig.out.proportion","w") as f:
	f.write( '\t'.join(OutputTable40head) )
	f.write( '\n' )
	wr = csv.writer(f, delimiter="\t")
	wr.writerows(OutputTable40)

log_file = (str(len(OutputTable40)) + " variants of population combinations saved to ibdmig.out.proportion.\n")
with open("ibdmig.out.log", "a") as f: f.write(log_file)
print (log_file)

####################################################################
#                                                                  #
# Getting the ibdmig.out.proportion.header [CombPopNameTableOut]   #
#                                                                  #
####################################################################


CombPopNameTableOut = []
for x in range(len(PopCombinations)):
	CombPopNameTableOutRow = []
	CombPopNameTableOutRow.append((PopCombinations[x]))
	CombPopNameTableOutRow.append(('_'.join(CombPopNameTable[x])))
	CombPopNameTableOut.append(CombPopNameTableOutRow)
#print ('CombPopNameTableOut')
#print (CombPopNameTableOut)

####################################################################
#                                                                  #
# Saving the ibdmig.out.proportion.header file                     #
#                                                                  #
####################################################################

with open("ibdmig.out.proportion.header","w") as f:
	wr = csv.writer(f, delimiter="\t")
	wr.writerows(CombPopNameTableOut)

log_file = ("Variants of population combinations with population names saved to ibdmig.out.proportion.header.\n")
with open("ibdmig.out.log", "a") as f: f.write(log_file)
print (log_file)





######################################################################
#                                                                    #
# Fullihg the ibdmig.out.proportion.gdist [PopCombinationsTableGdist]#
#                                                                    #
######################################################################

#print ("OutputTableGdist")
#print (OutputTableGdist)

row_len = 6 + len(UniqPopList)
#print ("row_len" + str(row_len))

for x in OutputTableGdist:
	ClstSize = x[4]
	ClstLength = x[row_len]
	string1 = ''
	for y in range(len(UniqPopList)):
		string0 = '0'
		if int(x[6+y]) > 0:
			string0 = '1'
		string1 = string1 + string0
	string2 = string1 + ClstSize
	index1 = (PopCombinationsTableInd.index(string2))
	PopCombinationsTableGdist[index1][3] = PopCombinationsTableGdist[index1][3] + ClstLength

#print("PopCombinationsTableGdist")
#print(PopCombinationsTableGdist)


OutputTable4  = []	
OutputTable40Gdist = []
for x in PopCombinations:
	OutputTable4 = []
	string2total_for_row = 0
	OutputTable4.append(x)
	for y in OutputTable3:
		string2 = str(x) + str(y)
		index1 = (PopCombinationsTableInd.index(string2))
		string2 = PopCombinationsTableGdist[index1][3]
		#print("string2 = " + str(string2))
		string2total_for_row = string2total_for_row + string2
		OutputTable4.append(str(string2))
	OutputTable4.append(str(string2total_for_row))
	#print("string2total_for_row")
	#print(str(string2total_for_row))
	OutputTable40Gdist.append(OutputTable4)
	#print("OutputTable4")
	#print(OutputTable4)

#print("OutputTable40Gdist")
#print(OutputTable40Gdist)


OutputTable50Gdist = []

for x in range(len(OutputTable40Gdist)):
	OutputRow50Gdist = []
	#print (OutputTable40Gdist[x])
	#print (OutputTable40[x])
	for y in range(len(OutputTable40Gdist[x])):
		#print ("y=" + str(y))
		#print (OutputTable40Gdist[x][y])
		#print (OutputTable40[x][y])
		if y > 0 and int(OutputTable40[x][y]) > 0 :
			z = Decimal(OutputTable40Gdist[x][y]) / int(OutputTable40[x][y])
		else:
			z = OutputTable40Gdist[x][y]
		#print (z)
		OutputRow50Gdist.append(str(z))
	OutputTable50Gdist.append(OutputRow50Gdist)

#print("OutputTable50Gdist")
#print(OutputTable50Gdist)


####################################################################
#                                                                  #
# Saving the ibdmig.out.proportion.gdist file                      #
#                                                                  #
####################################################################

with open("ibdmig.out.proportion.gdist","w") as f:
	f.write( '\t'.join(OutputTable40head) )
	f.write( '\n' )
	wr = csv.writer(f, delimiter="\t")
	wr.writerows(OutputTable50Gdist)

log_file = (str(len(OutputTable50Gdist)) + " variants of population combinations saved to ibdmig.out.proportion.gdist.\n")
with open("ibdmig.out.log", "a") as f: f.write(log_file)
print (log_file)






####################################################################
#                                                                  #
# Getting the p-val for ibdmig.out.affected.significant            #
#                                                                  #
####################################################################

##############################################
#                                            #
# 3 dimensions (1-case, 2-control, 3-ethnos) #
#                                            #
##############################################

#print ( "OutputTable2 =" )
#for aa in OutputTable2: print ( aa )

#print ( "clusters =" )
#for bb in clusters: print ( bb )

#print ( "OutputTable =" )
#for cc in OutputTable: print ( cc )

if len(OutputTable2) > 0:
	Dataset = []
	for x in UniqPopList:
		for y in IndPopList:
			if y[1] == x:
				Dataset.append(int(str(UniqPopList.index(x)+1) + y[2]))

	#print ( "Dataset = " + str(Dataset) )
	NumberOfPermutations = 100000
	P_perm = 1
	P_corr = 1
	OutputTable2RowIndex = 0

	for a in OutputTable2:
		ChromNum = a[0]
		ClstrNum = a[1]

		for x in clusters:
			if x[0] == ChromNum and x[1] == ClstrNum:
				ClstStr = []
				for y in range(len(x)):
					if y > 3: ClstStr.append(x[y])
		ClstStr.sort()
		#print ( "ClstStr = " + str(ClstStr) )

		ClstStrRow = []
		for x in OutputTable:
			if x[0] == ChromNum and x[1] == ClstrNum:
				ClstStrRow = x
		#print ( "ClstStrRow = " + str(ClstStrRow) )


		ClstStrRowsTotal = []
		for x in OutputTable:
			if x[4] == ClstStrRow[4]:
				ClstStrRowsTotal.append(x)
		#print ( "ClstStrRowsTotal = " + str(ClstStrRowsTotal) )
		#print ( "len(ClstStrRowsTotal) = " + str(len(ClstStrRowsTotal)) )

		ClstStrRowsExactMatch45Columns = []
		for x in OutputTable:
			if x[4] == ClstStrRow[4] and x[5] == ClstStrRow[5]:
				ClstStrRowsExactMatch45Columns.append(x)
		#print ( "ClstStrRowsExactMatch45Columns = " + str(ClstStrRowsExactMatch45Columns) )
		#print ( "len(ClstStrRowsExactMatch45Columns) = " + str(len(ClstStrRowsExactMatch45Columns)) )

		ClstStrRowInds = []
		for x in range(len(ClstStrRow)):
			if int(x) > 3: 
				ClstStrRowInds.append(x)
		#print ( "ClstStrRowInds = " + str(ClstStrRowInds) )


		ClstStrRowsExactMatchAllColumns = []
		for x in ClstStrRowsTotal:
			x_match = 1
			for Ind in ClstStrRowInds:
				if x[Ind] != ClstStrRow[Ind]:
					x_match = 0
			if x_match == 1:
				ClstStrRowsExactMatchAllColumns.append(x)

		#print ( "ClstStrRowsExactMatchAllColumns = " + str(ClstStrRowsExactMatchAllColumns) )
		#print ( "len(ClstStrRowsExactMatchAllColumns) = " + str(len(ClstStrRowsExactMatchAllColumns)) )


		ClstStrRowPopsStr = ""
		xxx = 0
		for x in ClstStrRow:
			if xxx > 5: 
				Str = "0"
				if int(x) > 0: Str = "1"
				ClstStrRowPopsStr = ClstStrRowPopsStr + Str
			xxx += 1
		#print ( "ClstStrRowPopsStr = " + str(ClstStrRowPopsStr) )


		ClstStrRowsPopCompositionMatch = []
		for yy in ClstStrRowsTotal:
			ClstStrRowPopsStr1 = ""
			xxx = 0
			for x in yy:
				if xxx > 5: 
					Str = "0"
					if int(x) > 0: Str = "1"
					ClstStrRowPopsStr1 = ClstStrRowPopsStr1 + Str
				xxx += 1
			if ClstStrRowPopsStr1 == ClstStrRowPopsStr: 
				ClstStrRowsPopCompositionMatch.append(yy)

		#print ( "ClstStrRowsPopCompositionMatch = " + str(ClstStrRowsPopCompositionMatch) )
		#print ( "len(ClstStrRowsPopCompositionMatch) = " + str(len(ClstStrRowsPopCompositionMatch)) )


		Clst = []
		for x in UniqPopList:
			for y in IndPopList:
				if y[1] == x:
					for z in ClstStr:
						if y[0] == z: Clst.append(int(str(UniqPopList.index(x)+1) + y[2]))
		Clst.sort()
		#print ( "Clst = " + str(Clst) )
		ClstSize = len(Clst)
		#print ( "ClstSize = " + str(ClstSize) )
		#print ( "NumberOfPermutations = " + str(NumberOfPermutations) )


		Clst1lettersList = []
		Clst1lettersStr  = ""
		for x in Clst: Clst1lettersList.append(str(x)[0])
		ClstUniqSet = sorted(set(Clst1lettersList))
		Clst1lettersStr = ''.join(sorted(set(Clst1lettersList)))
		#print ( "Clst1lettersStr = " + str(Clst1lettersStr) )


		CounterForClst = []
		CounterForTotalColumn = []
		for perm in range(NumberOfPermutations):
			CounterTotalRow = []
			for x in range(ClstSize): CounterTotalRow.append(random.choice(Dataset))
			CounterTotalRow.sort()
			if CounterTotalRow == Clst: CounterForClst.append(CounterTotalRow)

			CounterPermuted1lettersList = []
			CounterPermuted1lettersStr  = ""
			for x in CounterTotalRow: CounterPermuted1lettersList.append(str(x)[0])
			CounterPermutedUniqSet = sorted(set(CounterPermuted1lettersList))
			CounterPermuted1lettersStr = ''.join(sorted(set(CounterPermuted1lettersList)))
			#print ( "CounterPermuted1lettersStr = " + str(CounterPermuted1lettersStr) )
			#print ( "Clst1lettersStr = " + str(Clst1lettersStr) )
			#print ( str(CounterPermuted1lettersStr) + "=" + str(Clst1lettersStr) )
			if str(CounterPermuted1lettersStr) == str(Clst1lettersStr): CounterForTotalColumn.append(CounterTotalRow)

		#print ( "CounterForTotalColumn = " + str(CounterForTotalColumn) )
		#print ( "len(CounterForTotalColumn) = " + str(len(CounterForTotalColumn)) )
		#print ( "CounterForClst = " + str(CounterForClst) )
		#print ( "len(CounterForClst) = " + str(len(CounterForClst)) )

		Cell      = len(ClstStrRowsPopCompositionMatch )
		ColumnSum = len(ClstStrRowsTotal)
		ProportionEmpirical = Cell / ColumnSum

		CounterPermutedAllMatch  = len(CounterForClst)
		CounterPermutedSizeMatch = len(CounterForTotalColumn)
		ProportionPermuted = CounterPermutedSizeMatch / NumberOfPermutations

		#print ( "ProportionPermuted = " + str(ProportionPermuted) )

		if ProportionPermuted == 0:
			Alpha = ProportionEmpirical
		else:
			Alpha = ProportionEmpirical / ProportionPermuted

		P_perm =   CounterPermutedAllMatch / NumberOfPermutations
		P_corr = ( CounterPermutedAllMatch * Alpha ) / NumberOfPermutations

		OutputTable2[OutputTable2RowIndex].append(P_perm)
		OutputTable2[OutputTable2RowIndex].append(P_corr)
		OutputTable2[OutputTable2RowIndex].append(Cell)
		OutputTable2[OutputTable2RowIndex].append(ColumnSum)
		OutputTable2[OutputTable2RowIndex].append(ProportionEmpirical)
		OutputTable2[OutputTable2RowIndex].append(ProportionPermuted)
		OutputTable2[OutputTable2RowIndex].append(CounterPermutedAllMatch)
		OutputTable2[OutputTable2RowIndex].append(CounterPermutedSizeMatch)
		OutputTable2[OutputTable2RowIndex].append(Alpha)

		OutputTable2RowIndex = OutputTable2RowIndex + 1

####################################################################
#                                                                  #
# Saving the ibdmig.out.affected file with p-val                   #
#                                                                  #
####################################################################
	TableHeader.append('P_perm')
	TableHeader.append('P_corr')
	TableHeader.append('Cell')
	TableHeader.append('ColumnSum')
	TableHeader.append('ProportionEmpirical')
	TableHeader.append('ProportionPermuted')
	TableHeader.append('CounterPermutedAllMatch')
	TableHeader.append('CounterPermutedSizeMatch')
	TableHeader.append('Alpha')

	with open("ibdmig.out.affected","w") as f:
		f.write( '\t'.join(TableHeader) )
		f.write( '\n' )
		wr = csv.writer(f, delimiter="\t")
		wr.writerows(OutputTable2)
	log_file = ("p-values added to ibdmig.out.affected.\n")
	with open("ibdmig.out.log", "a") as f: f.write(log_file)
	print (log_file)


####################################################################

log_file = ("\nEnd Time          :  " + datetime.now().strftime('%Y-%m-%d %H:%M:%S') + "\n")
with open("ibdmig.out.log", "a") as f: f.write(log_file)
print (log_file)

