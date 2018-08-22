# IBDMig

IBDMig is a Python3 tool to assess the admixture process in mixed cohort via IBD sharing. IBDMig assesses IBD sharing for individuals of different ethnic origin in DASH (Gusev et al., 2011) generated IBD clusters. Thus it shows the haplotype contribution from one population to other. Furthermore, IBDMig detects IBD clusters enriched with patients of one/different ethnic origin.<br>
Dowmload IBDMig files from https://github.com/Ural-Yunusbaev/IBDMig/archive/master.zip

### Usage

<pre>./ibdmig.py 22 ibdmig.list mapfile.bim</pre>
where:<br>
22 - the number of chromosomes in according to number of DASH output files (clust_1.hcl ... clust_22.hcl);<br>
ibdmig.list - the file containing a list of individuals;<br>
mapfile.bim - the map/bim file with genetic distances (not mandatory).<br>
9  - the size threshold for affected polyethnic cluster (not mandatory, 9 if not defined)<br>
6  - the size threshold for affected monoethnic cluster (not mandatory, 6 if not defined)<br>

IBDMig generates the following output files:<br>
ibdmig.out.cluster_counts - counts of clusters for each populations combinations and cluster size category (see Output files examples);<br>
ibdmig.out.cluster_length - average length of haplotypes for each populations combinations and cluster size category (see Output files examples).<br>

### Output files examples

<pre>
cat ibdmig.out.cluster_counts
POPS	3	4	5	6	7	8	9	10	11	12	13	14	15	16	17	20	TOTAL
100	0	2958	870	227	91	38	13	7	3	0	0	1	0	0	0	0	4208
010	0	3697	829	224	66	18	12	3	0	0	0	0	0	0	0	0	4849
001	0	1593	282	54	14	2	1	0	0	0	0	0	0	0	0	0	1946
110	0	6627	2002	767	357	128	58	17	5	1	0	0	0	0	0	0	9962
101	1	11207	4096	1646	853	396	212	104	33	25	14	4	0	1	0	0	18592
011	0	15982	5757	2344	1190	472	273	108	29	11	6	2	1	0	0	0	26175
111	0	15367	8784	5327	3640	2042	1342	745	346	158	74	39	11	2	1	1	37879
</pre>
Counts of clusters for populations combination and cluster size categories.<br>
Rows are populations combinations, and columns are clusters sizes.<br>
The header of ibdmig.out.cluster_counts is following: <br>
POPS - populations combinations;<br>
4-20 - sizes of clusters;<br>
TOTAL - total number for the row.<br>
Populations combinations in column 1 presented in the file ibdmig.out.cluster_header.<br>

<pre>
cat ibdmig.out.cluster_length
POPS	4	5	6	7	8	9	10	11	12	13	14	15	16	17	20	TOTAL
100	3.8	3.3	3.0	3.2	2.8	3.1	2.7	3.5			3.3					3.6
010	2.5	2.4	2.3	2.3	2.2	2.0	2.0									2.5
001	3.1	2.7	2.7	2.3	2.3	3.6										3.0
110	2.6	2.5	2.4	2.3	2.3	2.1	2.0	2.3	2.0							2.5
101	3.1	2.9	2.7	2.6	2.6	2.5	2.6	2.3	2.4	2.3	2.6		2.3			3.0
011	2.6	2.4	2.3	2.2	2.2	2.1	2.0	2.1	2.0	2.0	2.3	2.0				2.5
111	2.6	2.5	2.4	2.3	2.3	2.1	2.2	2.1	2.0	2.1	2.0	1.8	2.1	2.0	2.1	2.5
</pre>
The average length of haplotypes for populations combination and cluster size categories.<br>
Rows are populations combinations, and columns are clusters sizes.<br>
The header of ibdmig.out.cluster_counts is following: <br>
POPS - populations combinations;<br>
4-20 - sizes of clusters;<br>
TOTAL - average for the row.<br>
Populations combinations in column 1 presented in the file ibdmig.out.cluster_header.<br>

<pre>
cat ibdmig.out.cluster_header
100    pop1
010    pop2
001    pop3
110    pop1_pop2
101    pop1_pop3
011    pop2_pop3
111    pop1_pop2_pop3
</pre>

### Input files examples

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
</pre>
Columns: individual ID, source population, phenotype.<br>
The maximum number of source populations is 7.

<pre>
head -n 3 clust_1.hcl
c1    16504399    17593685    19N 19N.0    19N 19N.0    182A 182A.0    182A 182A.0    66i 66i.1    66i 66i.1    153A 153A.1    153A 153A.1
c2    16504399    17799529    62BB 62BB.0    62BB 62BB.0    55k 55k.0    55k 55k.0    190k 190k.0    190k 190k.0    51A 51A.1    51A 51A.1
c3    16504399    17823261    164B 164B.0    164B 164B.0    38BO 38BO.1    38BO 38BO.1    36i 36i.1    36i 36i.1    100k 100k.1    100k 100k.1
</pre>
For details see http://www1.cs.columbia.edu/~gusev/dash/

<pre>
head -n 3 mapfile.bim
1       rs3094315       0.48877594      752566  G       A
1       rs12562034      0.49571378      768448  A       G
1       rs12124819      0.49944228      776546  G       A
</pre>
For details see http://zzz.bwh.harvard.edu/plink/data.shtml#map

### Additional output files

<pre>
head -n 4 ibdmig.out.cluster_list
CHR	CLUSTER	START	END	SIZE	AFFECT	Pop1	Pop2	Pop3	LENGTH_cM	START_cM	END_cM
1	c1	1152631	2996602	5	3	0	3	2	0.0	0.0	0.0
1	c2	1310924	3147030	4	0	3	0	1	0.0	0.0	0.0
1	c3	1493727	2754512	4	1	0	2	2	0.0	0.0	0.0
</pre>
Columns are folowing:<br>
Chromosome number;<br>
Cluster identifier;<br>
Cluster start position;<br>
Cluster end position;<br>
Cluster size;<br>
The number of affected individuals (patients);<br>
The number of individuals from Pop1;<br>
The number of individuals from Pop2;<br>
The number of individuals from Pop3;<br>
Genetic length in centimorgans.<br>
Genetic distanse for start position in centimorgans.<br>
Genetic distanse for end position in centimorgans.<br>
<pre>
cat ibdmig.out.cluster_list.end
CHR	CLUSTER	START	END	SIZE	AFFECT	Pop1	Pop2	Pop3	LENGTH_cM	START_cM	END_cM
max	-	-	-	15	9	6	8	6	0	0	0
min	-	-	-	4	0	0	0	0	0	0	0
mean	-	-	-	6.4	3.6	2.0	2.7	1.6	0.0	0.0	0.0
</pre>

### Contact
Ural Yunusbaev<br>
uralub@gmail.com<br>

### Citing
this tool was developed for<br>
Ural Yunusbayev, Mait Metspalu, Bayazit Yunusbayev. (2018). Reconstructing recent population history while mapping rare variants.<br>

### References
Gusev, A., Kenny, E. E., Lowe, J. K., Salit, J., Saxena, R., Kathiresan, S., Altshuler, D., Friedman, J., Breslow, J., Pe’er, I. (2011). DASH: a method for identical-by-descent haplotype mapping uncovers association with recent variation. American Journal of Human Genetics, 88(6), 706–717.
