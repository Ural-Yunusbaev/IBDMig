# IBDMig
Python tool for IBD clusters admixture estimation

By Ural Yunusbaev 2018

The IBDMig uses DUSH output IBD clusters to calculate the number of clusters including individuals of one ethnic origin and of different ethnic origin. It creates a table where rows are populations combinations and columns are cluster sizes. The table shows the counts of monoethnic and polyethnic clusters. This reflects admixture process in studying population. IBDMig calculates multi-IBD sharing between individuals of different ethnic origin. Thus it shows the haplotype contribution from one population to other. Furthermore the IBDMig estimates average length of haplotypes in each cluster size bins.

To start type: <br>
./ibdmig.py 22 ibdmig.list mapfile

where:<br>
22 - the number of chromosomes in according to DUSH output files (clust_1.hcl ... clust_22.hcl);<br>
ibdmig.list - the file containing a list of individuals with following columns: ind_id, population, phenotype;<br>
mapfile - the map/bim file with genetic distances (not mandatory).<br>

<pre>
Table 1.
POPS	4	5	6	7	8	9	10	11	12	13	14	15	16	17	20	TOTAL
Bas	2958	870	227	91	38	13	7	3	0	0	1	0	0	0	0	4208
Rus	3697	829	224	66	18	12	3	0	0	0	0	0	0	0	0	4849
Tat	1593	282	54	14	2	1	0	0	0	0	0	0	0	0	0	1946
Bas_Rus	6627	2002	767	357	128	58	17	5	1	0	0	0	0	0	0	9962
Bas_Tat	11207	4096	1646	853	396	212	104	33	25	14	4	0	1	0	0	18591
Rus_Tat	15982	5757	2344	1190	472	273	108	29	11	6	2	1	0	0	0	26175
Bas_Rus_Tat	15367	8784	5327	3640	2042	1342	745	346	158	74	39	11	2	1	1	37879
TOTAL	57431	22620	10589	6211	3096	1911	984	416	195	94	46	12	3	1	1	103610
</pre>
Input files examples:<br>
ibdmig.list - the list of individuals with following tab delimeted columns: ind_id, population, phenotype<br>
<pre>
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
* - columns are folowing: individual ID, population, phenotype.<br>
**- maximum number of different populations is 7.
</pre>

head -n 3 clust_1.hcl<br>
c1	16504399	17593685	19N 19N.0	19N 19N.0	182A 182A.0	182A 182A.0	66i 66i.1	66i 66i.1	153A 153A.1	153A 153A.1<br>
c2	16504399	17799529	62BB 62BB.0	62BB 62BB.0	55k 55k.0	55k 55k.0	190k 190k.0	190k 190k.0	51A 51A.1	51A 51A.1<br>
c3	16504399	17823261	164B 164B.0	164B 164B.0	38BO 38BO.1	38BO 38BO.1	36i 36i.1	36i 36i.1	100k 100k.1	100k 100k.1<br>

* - for details see http://www1.cs.columbia.edu/~gusev/dash/


head -n 3 mapfile.bim<br>
1       rs3094315       0.48877594      752566  G       A<br>
1       rs12562034      0.49571378      768448  A       G<br>
1       rs12124819      0.49944228      776546  G       A<br>
* - for details see http://zzz.bwh.harvard.edu/plink/data.shtml#map


Output files are following:<br>
ibdmig.out.affected<br>
ibdmig.out.proportion<br>
ibdmig.out.total<br>

Output files examples:

cat ibdmig.out.proportion<br>
POPS	4	5	6	7	8	9	10	11	15	16	TOTAL<br>
100	1	1	0	0	0	0	0	0	0	0	2<br>
010	1	1	1	0	0	0	0	0	0	0	3<br>
001	1	0	0	0	0	0	0	0	0	0	1<br>
110	1	1	1	0	0	0	0	0	0	0	3<br>
101	1	1	1	0	0	0	0	0	0	0	3<br>
011	1	1	1	0	1	0	0	0	0	0	4<br>
111	1	1	1	1	1	1	1	1	1	1	10<br>

* The header of ibdmig_proportion.out is following: <br>
POPS - variants of population combinations;<br>
4-16 - sizes of IBD clusters;<br>
TOTAL - total number for row<br>
** populations combinations in column 1 presented in the file ibdmig.out.proportion.header.<br>

cat ibdmig.out.proportion.header<br>
100	pop1<br>
010	pop2<br>
001	pop3<br>
110	pop1_pop2<br>
101	pop1_pop3<br>
011	pop2_pop3<br>
111	pop1_pop2_pop3<br>


head -n 4 ibdmig.out.total<br>
CHR	CLASTER	START	END	SIZE	AFFECT	Pop1	Pop2	Pop3<br>
1	c1	165043	175936	4	4	1	3	0<br>
1	c2	165043	177995	6	2	4	0	2<br>
1	c3	165043	178232	4	1	2	0	2<br>

* - columns are folowing:<br>
Chromosome number<br>
Cluster identifier<br>
Cluster start position<br>
Cluster end position<br>
Cluster size<br>
The number of affected individuals<br>
The number of individuals from Pop1<br>
The number of individuals from Pop2<br>
The number of individuals from Pop3<br>

cat ibdmig.out.total.end<br>
CHR	CLASTER	START	END	SIZE	AFFECT	Bas	Rus	Tat<br>
max	-	-	-	6	4	4	3	2<br>
min	-	-	-	4	1	1	0	0<br>
mean	-	-	-	5	2	2	1	1<br>

