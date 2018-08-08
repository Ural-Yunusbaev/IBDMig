# IBDMig
Python tool for IBD clusters admixture estimation

By Ural Yunusbaev 2018

The IBDMig uses DUSH output IBD clusters to calculate the number of clusters including individuals of one ethnic origin and of different ethnic origin. It creates a table where rows are populations combinations and columns are cluster sizes. The table shows the counts of monoethnic and polyethnic clusters. This reflects admixture process in studying population. IBDMig calculates multi-IBD sharing between individuals of different ethnic origin. Thus it shows the haplotype contribution from one population to other. Furthermore the IBDMig estimates average length of haplotypes in each cluster size bins.

To start type: 
./ibdmig.py 22 ibdmig.list mapfile

where:
22 - the number of chromosomes in according to DUSH output files (clust_1.hcl ... clust_22.hcl);
ibdmig.list - the file containing a list of individuals with following columns: ind_id, population, phenotype;
mapfile - the map/bim file with genetic distances (not mandatory).

Input files examples:
ibdmig.list - the list of individuals with following tab delimeted columns: ind_id, population, phenotype
head ibdmig.list
10BO	pop1	1
103B	pop1	1
9i	pop1	2
88N	pop1	2
9RE	pop2	1
98RE	pop2	1
103A	pop3	2
102N	pop3	2
101N	pop3	2
100N	pop3	2

* - columns are folowing: individual ID, population, phenotype.
**- maximum number of different populations is 7.


head -n 3 clust_1.hcl
c1	16504399	17593685	19N 19N.0	19N 19N.0	182A 182A.0	182A 182A.0	66i 66i.1	66i 66i.1	153A 153A.1	153A 153A.1
c2	16504399	17799529	62BB 62BB.0	62BB 62BB.0	55k 55k.0	55k 55k.0	190k 190k.0	190k 190k.0	51A 51A.1	51A 51A.1
c3	16504399	17823261	164B 164B.0	164B 164B.0	38BO 38BO.1	38BO 38BO.1	36i 36i.1	36i 36i.1	100k 100k.1	100k 100k.1

* - for details see http://www1.cs.columbia.edu/~gusev/dash/


head -n 3 mapfile.bim
1       rs3094315       0.48877594      752566  G       A
1       rs12562034      0.49571378      768448  A       G
1       rs12124819      0.49944228      776546  G       A

* - for details see http://zzz.bwh.harvard.edu/plink/data.shtml#map


Output files are following:
ibdmig.out.affected
ibdmig.out.proportion
ibdmig.out.total

Output files examples:

cat ibdmig.out.proportion
POPS	4	5	6	7	8	9	10	11	15	16	TOTAL
100	1	1	0	0	0	0	0	0	0	0	2
010	1	1	1	0	0	0	0	0	0	0	3
001	1	0	0	0	0	0	0	0	0	0	1
110	1	1	1	0	0	0	0	0	0	0	3
101	1	1	1	0	0	0	0	0	0	0	3
011	1	1	1	0	1	0	0	0	0	0	4
111	1	1	1	1	1	1	1	1	1	1	10

* The header of ibdmig_proportion.out is following: 
POPS - variants of population combinations;
4-16 - sizes of IBD clusters;
TOTAL - total number for row
** populations combinations in column 1 presented in the file ibdmig.out.proportion.header.

cat ibdmig.out.proportion.header
100	pop1
010	pop2
001	pop3
110	pop1_pop2
101	pop1_pop3
011	pop2_pop3
111	pop1_pop2_pop3


head -n 4 ibdmig.out.total
CHR	CLASTER	START	END	SIZE	AFFECT	Pop1	Pop2	Pop3
1	c1	165043	175936	4	4	1	3	0
1	c2	165043	177995	6	2	4	0	2
1	c3	165043	178232	4	1	2	0	2

* - columns are folowing:
Chromosome number
Cluster identifier
Cluster start position
Cluster end position
Cluster size
The number of affected individuals
The number of individuals from Pop1
The number of individuals from Pop2
The number of individuals from Pop3

cat ibdmig.out.total.end
CHR	CLASTER	START	END	SIZE	AFFECT	Bas	Rus	Tat
max	-	-	-	6	4	4	3	2
min	-	-	-	4	1	1	0	0
mean	-	-	-	5	2	2	1	1

