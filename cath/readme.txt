--------------------------------------------------------------------------------

OK, from the CATH files, what we want?

We want domains, the 3D coordinates of the atoms in the domains. The coordinates
will be grouped. For each residue, the back-bone and side-chain centres will be 
caculated. 

The residue contacts maps will then be obtained from the coordinates of the 
centres. 

We want all residues in the domains to be associated with the full sequences of 
the protein. More precisely, we need the breaks in the PDB chains to be filled
and mapped.  

It means we need to get related sequences from PsiBlast output and then generate
better multiple alignments using Muscle and Kalign.

And from PSIBlast we can obtain the PSSMs, position specific scoring matrices. 

We might need to calculate phylogenetic trees from the multiple alignments.

--------------------------------------------------------------------------------
1. Get the domain files

#download the Cath domains from 
 wget http://release.cathdb.info/v3.5.0/CathDomainPdb.v3_5_0.tgz
  
# get the index from 
 wget http://release.cathdb.info/v3.5.0/CathDomainList.S100.v3.5.0
 
# remove the transmembrane domains from http://www.cathdb.info/sfam/membrane/
# (69 H transmem families)
awk '{print $2 "." $3 "." $4 "." $5, $0 }' index/CathDomainList.S100.v3.5.0 \
 > t.in
grep -w -v -f  transmem.ls t.in | awk '{$1="";print $0}' \
 > index/CathDomainList.S100.v3.5.0
#44376 domains left, 1492 domain removed for being transmembrane

tar xvzf  CathDomainPdb.v3_5_0.tgz
mkdir pdb

# we have all domains in the index S100
awk '{print "mv dompdb/" $1 " pdb/"}' index/CathDomainList.S100.v3.5.0 > t.ls
source t.ls

--------------------------------------------------------------------------------

2. check residue names are of the standard types (some UNK, PCA, ASX etc.)
 --- chkResName.cpp
 
 Done. 
 
 Some (terminal) UNK segments are removed. 
 3fi2A02 3fv8A02 1gw5B00 1gw5A00 1vkvB01 1vkvA01 3ed7A00 3bcnA00 3bcnB00 2pscA00
 
 
 in 3e2oA02, 235 X to D
 in 156b000, 21,22,23 XXX to DNA 
 in 2fmdA00, 204,205,206 XPX to LPD
 in 2v5iA00, 206 X to G
 in 3jxvA02, 284 X to G
 in 2atcA02, 234 X to N
 in 3atcA02, 234 X to N
 in 3atcC02, 234 X to N
 in 1o0eA00, 117 X to N
 in 1o0eB00, 117 X to N
 in 2pscA00, 199 X to A
 
 Broken domains 2rxn000 3bcl000 1kp0A01 1kp0B01 1kp0A02 1kp0B02 are deleted.
 
--------------------------------------------------------------------------------
3. check strange atoms in residues, removing AE1, AD3, OXT and H atoms.
mkdir pdb2

 --- rmHAtoms.cpp
 
 We don't want the H atoms... 2788 domains with H and other characters in the 
 leading position of the atom name, has those atoms removed. 
 
 Then only the regular heavy atoms are allowed. 
 
mv pdb2/* pdb

--------------------------------------------------------------------------------
4. check the negtively occuring atoms. 
 --- chkNegOccu.cpp
 
 From 1mwcA00, 147 NZ and CE deleted. Clearly guessed coordinates.
 From 1m6mA00, 147 NC and CE deleted. Same reason.
--------------------------------------------------------------------------------

5. check the alternative location of atoms
 --- chkMultiLoc.py
 Now check for different AltLoc marks across the whole Domains. 
 4721 domains index/multiLoc.ls, 
 3562 PDBs, downloaded from rcsb to rcsb
 
 In the 4721 there are 4573 with AB or ABC altLocs
 
 --- replaceAltLocAB.cpp

 First remove the 'B' and 'C' atoms when there is a 'A' there. It should fix 
 most cases. (4483 out of the 4573)
 
 Now we have 238 domains to be fixed
 
 Use B and nothing else, fixed 155, 83 to be fixed.
 
 Use 1 instead of 2, fixed 6, 77 left
 
 remove 2v93A01 2v93A02 
 
 manually fix the 77
 
 now many atoms have multiple entries (from AB -> AA etc)
 --- rmSameAtoms.cpp
 
 Check if there are multiple ATOM entries from the same residue with the same 
 atom name.
 
 --- chkAltlocRes.cpp
 
 Nothing found there. 

--------------------------------------------------------------------------------

5. check atoms are not too close to each other. 
 --- chkAtomDis.cpp
 Manually fixed some domains.
 
--------------------------------------------------------------------------------

6. check atoms in the same residue should be of appropoiately close
 --- chkResDis.cpp
 Manually fixed some domains.
--------------------------------------------------------------------------------

7. check residues (missing CA, N or C atoms)
 --- chkCANC.cpp
 deleted 1w5cC02 1grlA03
 there are small programs to help: remove ophan residue with only one atom, 
 remove the first or the last residues if they are broken etc.

--------------------------------------------------------------------------------

8. check for segments and breaks of the mainchain. 
 --- chkBreak.cpp 
    
    The distance between the C and the N atoms of the next residue is mostly
    between 0 and 2.5A (99%) (jumped from 1.67 to 2.68). 
    That's why I use 2.5A as a threshold here. 
    
    out of 44366 domains. 
    22413 domains with no break,
    9839 with one, 4387 with two, 
    
    This is bad. More than half of the domains have gaps in them. We need to 
    use the gaps. 
    
--------------------------------------------------------------------------------

9. basic blast search of the atom sequences. (need biopython, or not)
 --- get nr protein sequence database
 --- filter it with Davic Jone's tiny program pfilt. 
 --- search the sequences with the filtered database
--------------------------------------------------------------------------------
header file in the brc cluster:

#!/bin/sh
#$ -S /bin/sh
#$ -o  /scratch/home/klinbrc/dev/cath
#$ -e  /scratch/home/klinbrc/dev/cath

-------------------------------------------------------------------------------- 
~/bin/blastp -db ~/scratch/dev/nr/nrfilt \
 -outfmt "7 sseqid sstart send slen pident" \
 -query ~/scratch/dev/cath/seq/1waaC00 \
 -out ~/scratch/dev/cath/1waaC00_blast.hits

 I used the most basic option. Anyway, this is to get the reference sequences
 only. 
 
 blastp -db ~/dev/nr/nrfilt -outfmt "7 sseqid sstart send slen pident" 
--------------------------------------------------------------------------------

10. get the blast sequences. 

 ls blast_hits/ > dom.ls
 mkdir blast_seqs
 --- parseBlastHits.py
 
 use only the hits with over 50% sequence identity. 
 
 For 40056 domain sequences, 5980684 hits. That's 149 hits per domain. 
 Without the sequence identity threshold, there are 227 hits per domain.
 
 The old entry_batch way: (way faster)
 awk '{print "~/bin/blastdbcmd -db ~/scratch/dev/nr/nr -entry_batch \
~/scratch/dev/cath/blast_gids/" $1 \
" > ~/scratch/dev/cath/blast_seqs/" $1 }' dom.ls > t.sh

 t25.sh twice and t63.sh once, 
 BLAST engine error: Failed to parse sequence range 
 (start cannot be larger than stop)
 the hits are 
1larA02 gi|229442495|gb|AAI72930.1|     1136    1431    517     79.73
1larB01 gi|229442495|gb|AAI72930.1|     889     1134    517     71.14
2fh7A02 gi|229442495|gb|AAI72930.1|     1136    1431    517     79.05
 that's blast's problem. The record gi|229442495 is wrong
 remove them .
--------------------------------------------------------------------------------

11. get the sequences to be aligned. 

 The blast hit sequences may have long long long names, shrinking the names to 
 GI ids etc helps. 
 
 Remember we had difficulty aligning some (2000?) of them last time. Sometime 
 kalign complains about some 'J' in the sequences. Checking it. 

#########################################################################
#!/usr/bin/python

indexF = open('/home/klinbrc/scratch/dev/cath/dom.ls')

for domName in indexF:
    domName=domName.strip();
    
    inF=open('/home/klinbrc/scratch/dev/cath/blast_seqs/'+domName);
    outF=open('/home/klinbrc/scratch/dev/cath/blast_seqs2/'+domName,'w');
    outF.write(inF.readline())
    for line in inF:
        if line.startswith('>'):
            line=line.split()[0]+'\n'
            
        outF.write(line)
    inF.close()
    outF.close()
    
indexF.close()
###########################################################################
 
--------------------------------------------------------------------------------

12. get sequence multiple alignments 

 kalign blast_seqs2/102mA00 > kaln/102mA00
 
 muscle -in blast_seqs2/102mA00 -out maln/102mA00
 
 No error messages.
 
 In case there were the domain sequence alone, no blast hits, the output will 
 be a empty file for kalign, and a file with the domain sequence for muscle.

 Use muscle alignment anyway. They are both reasonably fast. (No we don't have 
 the extremly long complete sequences.)
 
--------------------------------------------------------------------------------

13. get the pairwise alignment of the domain and the consensus sequences. 
 -- getConAlign.cpp
 
 Refined the parameters. Realised we had very bad kalign alignment. Mainly 
 because we used the huge, complete sequence from NR. Went back to grab Blast 
 subject sequences with different options.
 
 Anyway, changed the initial blast searching output format and the 
 parseBlastHits.py script, we have the inputs for getting multiple alignments. 
 
 Out of the 40056 Sequences, 37736 are having different consensus sequences.
 That's 94.2%. 
 
 Check: 
 1. in the pairwise alignment, the domain sequences, after removing the gaps, 
 should be identical to the input domain sequences.
 2. the consensus sequences, should be identical to the domain sequence when 
 it's gaps are removed. 
 3. the consensus sequences have their extra gaps. These gaps should correspond
 to the backbone-breaks in the domains... Not always true. Should be often true.
 
 Did some checking and for the consensus alignment with many gaps, the domains 
 do have quite some backbone breaks. 
 
 On average, there are 1.22 gap breaks in an alignment. It sums to 48928 breaks.
 48% of domains have no break, 25% with one, 11% with two, 6% with three 
 and 4% with four.
  
 This is less than the backbone break we have... but very similar numbers. 
 
 After removing the 23 domains which doesn't produce PSSM in the next step, 
 there are 40033 domains left. 
 
 In these domains, we found 51379 backbone breaks and 48926 consensus sequence
 alignment breaks. 36885 backbone breaks are assigned to sequence breaks. 
 
 Venn diagram:
 
 (backbone breaks 14494 ( 36885 ) 12041 sequence breaks)
 
 If we had a backbone break, 72% of chance it will be assigned to a sequence
 break. 
 If we had a sequence break, 75% of chance it will correspond to a backbone
 break.
 
 We do over-predict some sequence breaks, but I think that's acceptable.
 
--------------------------------------------------------------------------------
14. get the PSSMs of the consensus sequences

 ~/bin/psiblast -db ~/scratch/dev/nr/nrfilt -num_iterations 2 \
  -query ~/scratch/dev/cath/cseq/1d66B02  \
  -out_ascii_pssm ~/scratch/dev/cath/cpssm/1d66B02
  
 Some 23 PSSMs are not generated. 
 We have five error messages like: Error: basic_string::substr
 For the other 18, no idea. 
 
 The 23 are 
  1d66B02, 1ex4B02, 1g3iF02, 1g3iW02, 1gmnB01, 1hf9A00, 1iodG00, 1kyiS02
  1larA02, 1larB01, 1mkmA02, 1nl0G00, 1r30A00, 2fh7A02, 2jwaA00, 2k1aA00
  2pp6A01, 3bxjA01, 3bxjB01, 3cj8A01, 3cj8B01, 3cj8C01, 4hb1A00

 Removed the 23 from domain index.
 
--------------------------------------------------------------------------------
15. give residues new residue numbers, according to the sequence alignment.

 --- newResidueNum.cpp
 
 Remove the chain insert id and alternative location id, change the old residue
 numbers to the positions at the consensus alignments.
 
 (It seems the program 'stride' doesn't like alternative locations.)
 
 Obtained 40033 domains in the cdom directory. (The directories are all 'c' 
 something after I have the consensus sequences.)
 
 The file domRes.map in the index directory have all the mapping between the old
 CATH (PDB) domain residue number (and Chain insert id) and new residue number. 
 
 Get the CA only files into the directory cdom_cas as well.
 
--------------------------------------------------------------------------------
16. STRIDE assignment of secondary structures.
 
 stride cdom/1waaC00 > ssdef/1waaC00
 
 Didn't use DSSP because DSSP assigns the coil state to the terminal residues. 

 Error message: 
 Residue SER 109 of chain cdom/1c3cA02A is involved in 6 hydrogen bonds (5 are allowed)
 Residue ASP 195 of chain cdom/1gegA00A is involved in 6 hydrogen bonds (5 are allowed)
 Residue SER 109 of chain cdom/2pfmA02A is involved in 6 hydrogen bonds (5 are allowed)
 Residue ASP 40 of chain cdom/2q3sF00F is involved in 6 hydrogen bonds (5 are allowed)
 Residue SER 63 of chain cdom/2r8oA01A is involved in 6 hydrogen bonds (5 are allowed)

 Manual check. The residues were rightly assigned anyway.
 
--------------------------------------------------------------------------------
17. get residue side-chain and backbone centers. 
 --- getS35Ls.cpp
 
 Get the sequence family representatives (11065 domains). 11330 is the Cath S35
 number. We lost 265 in the QC. Not too bad. 
 
 Or the 2477 H families should be enough? 
 
 Using the S35 families, 0.8% of all residues miss some atoms.
 
 --- getCbProject.cpp
 
 guess the location of side-chain centres using the CA, C and N atoms.
 
 --- getCb.cpp
 write the real and guessed cb with c atoms to the new director ccbc/ 
 
--------------------------------------------------------------------------------
18. get the contact maps

 --- getTet.cpp
 put the contact maps to the directory conDef
 --- mapToCseq.cpp
 Change the residue contact maps into the consensuse sequence coordinates.
 
--------------------------------------------------------------------------------

7. start nn training of YASPIN and CAO exposure prediction

In 5090 domains (681152 residues), we got 3051987 edges
 from Delaunay tetrahedralizations. That's about 9 edges for each residue. 
 
However, the edges are not necessary residue-residue contacts. 
Some are too long. They are on the out-side boundries of domains. 

Anyway, in the 3051987 edges, we have every posible contacts. 

Now we need to make the decision on where to cut the long edges. 

Group edges according to their length, look at the similarities of compositions. 

../res/edgeCompTree.fig

Now set the distance threshold to 8 A. 

2028480 contacts. 6 edges for each residues. 

look at the composition of contacts of local and global contacts. 

../res/disCompTree.fig

It seems dis 1, 2, 3, 4 and beyond are 5 different groups. 

So, we have CAO1, CAO2, CAO3, CAO4 and CAO5 five different series of matrices 
to estimate now. 

Now the contact numbers are 
CAO1 : 324058, 0.95 per residue 
CAO2 : 186660, 0.55 per residue 
CAO3 : 234152, 0.69 per residue 
CAO4 : 163284, 0.48 per residue 
CAO5 :1120326, 3.29 per residue 

from 9 million psiblast and 12,000 sap alignments, we have 
955722752 aligned contacts for c1
572853766 aligned contacts for c2
690025576 aligned contacts for c3
473712990 aligned contacts for c4
3.3621e+09 aligned contacts for c5

for sequence identity 0-5%, 5-10%, 10-15%, ... , 95-100%

2.2536e+05, 1.2339e+06, 8.9735e+06, 1.3095e+08, 4.9954e+08, 
9.5952e+08, 1.1759e+09, 9.7143e+08, 6.778e+08, 4.2939e+08, 
3.2877e+08, 2.0485e+08, 1.4574e+08, 1.1389e+08, 8.1308e+07, 
6.9563e+07, 5.8878e+07, 4.5562e+07, 5.0813e+07, 1.0013e+08

mainly on 25-50%

regression of number of contacts against number of residues:
number_residue = 17.931 + 0.29086 * number_contact
number_contact = -57.113 + 3.4043 * number_residue 
Correlation coefficient			 = 0.9950772

10. 5416 segments, 727775 residues, 2168190 contacts, 6 edges on every residues,
9636185 psiblast alignment. 54948*2=109896 sap alignments. 

11. Remove total identical aligned sequences (90769).
 Removed unaligned sequences (6). 
 9655306 alignments left
 
13. from the dssp files, get secondary structure definition. 

14. changed the first 'CHH' coil definition to 'HHH' and 'CEE' to 'EEE'
    changed the last 'HHC' to 'HHH' and 'EEC' to 'EEE'
    
15. sumE = 166965 (22.9%) E: 6862  EE: 3825*2=7650
    sumH = 257833 (35.4%)
    sumC = 302977 (41.6%)
    after changing the tails, about 27% of tails are now 'E' or 'H'. Anyway, 
    better than 0.

