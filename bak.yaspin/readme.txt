--------------------------------------------------------------------------------
                                YASPIN 
                            
    Yet Another Secondary structure Predictor usIng Neual network
--------------------------------------------------------------------------------

The method was published in Bioinformatics. 21(2):152-9. (Epub 2004 Sep 17). 

It is a hidden neural network model for protein secondary structure prediction. 

Three-states secondary structures are re-classified to 7 states. Helix begins 
and ends, strand begins and ends are seperated from the helix and stand classes.

A back-propagationly supervised feed-forward neural network was trained for the
prediction of emission probabilities of the 7 states. 

A hidden Markov model was then applied on the p-values to obtain the likely 
transimissions between states. 

--------------------------------------------------------------------------------

I think it's time to retrain the model with current data. It has been 7 years
scince the method was published. PDB, CATH and SCOP coverages should have 
improved a lot. We can have some more complete training set now. 

I have some unimplemented ideas. One is the building of consensus sequences for 
query sequences. The query sequence might have some unknown gaps. The consensus
alignment might be used to find and fill them. The sequences beyond the 
terminals of the query sequence is unknown to us. By applying the consensuse, we
might be able to get good PSSMs (Position Specific Scoring Matrices) for some of
the wing amino acids. 

The idea was tried in the cathsrc scripts. I am reasonably happy about the 
outcome. Now I need to check if the PSSMs are better. So far I'm very confident
about that. 

Another idea is to apply the second order hidden Markov model and include burial
states. 

This might be difficult to implement. That's not too bad. 

This will lead to much more states. That might be very bad. Too many states 
means overfitting? 

--------------------------------------------------------------------------------
1. copy the index, caln, cpssm, condef and ssdef directories from cath
 
 Use the 11065 domains of the S35 families, 
 or the 2477 domains of the H familes.
 
 the sizes of domains: 
 min: 16; max: 1146; mean: 146.849; sdv: 83.2491
 the sizes of consensus sequences:
 min: 33; max: 1223; mean: 172.367; sdv: 83.6271

 In the H families: 
 min: 16; max: 759; mean: 131.966; sdv: 80.0986
 min: 38; max: 760; mean: 156.761; sdv: 80.9108

--------------------------------------------------------------------------------
2. check the distribution of values in the PSSM

 --- chkPSSM.cpp

 We have 11065 domains, 1907241 residues in the consensus sequences. 
 In the PSSMs: 
 min: -16; max: 15; mean: -1.60253; sdv: 2.82846
 
 However, instead of 2.828, I still tend to use 5 for scaler. 
 And the mean is -1.6 doesn't matter here either. 

--------------------------------------------------------------------------------
3. get secondary structure definition. 

 I was using STRIDE. However, in addition to the general 8 states 
 (H, G, E, B, I, S, T, -), there are a few residues assigned to 'b' by STRIDE.
 I manually checked what they are. Some of them are assigned to B by DSSP, but 
 most of them can be very different. 
 
 I realised that STRIDE is very different from DSSP. And since everybody else
 was using DSSP, I was asking for trouble there. 
 
 Switched to dssp2 downloaded from 
 ftp://ftp.cmbi.ru.nl/pub/software/dssp/dssp-2-linux-i386.gz
 
 It worked for most of the domains, but complained for 3gylB01, saying
 "empty protein, or no valid complete residues". 
 
 Used the old DSSP cmbi version on that domain, no problem. 
 
 Checked the domain with rasmol. The first residues do have small atom numbers.
 But still can't see anything specially bad there.  
 
 --- chkDSSP.cpp
 
 Found 91 AA in 89 domains are missing DSSP assignments. Need to fix the domain
 file, the caln file and the PSSM file and the bloody condef files? 
 
 Found the 91 AA often have incomplete residue backbones. It means the O atom
 was missing. Two of the domains (1o4xA02 and 2pw8I00) are wrong or broken from 
 the start. Anyway removed the 89 domains. 
 
 lost 22 in H families. We now have 10976 domains in 2455 H families.
 
 From 8 states DSSP definitions to 3 states definitions:
 H,G -> H
 E,B -> E
 I,S,T,- -> C
 
 For the 1611638 residues with SS definitions in the S35 domains:
 669247 are C 41.5%
 584040 are H 36.2%
 358351 are E 22.2%
 
 For the 324117 residues in the CATH H domains:
 132229 are C 40.8%
 121914 are H 37.6%
 69974  are E 21.6%
 
 From 3 states to 7 states:
 CHH -> Hb -> G
 EHH -> Hb -> G
 HHC -> He -> I
 HHE -> He -> I
 HEE -> Eb -> D
 CEE -> Eb -> D
 EEC -> Ee -> F
 EEH -> Ee -> F
 
 It should be like CHHHHCCCEEEEC to CGHHICCCDEEFC. 
 
 In the S35 familes:
 
 669247  C 41.5%
 62042   G 3.85%
 459986  H 28.5%
 62012   I 3.85%
 64324   D 3.99%
 229692  E 14.3%
 64335   F 3.99%
 
 In the H families: 
 
 132229 C 40.8%
 12297  G 3.79%
 97327  H 30.0%
 12290  I 3.79%
 12075  D 3.73%
 45815  E 14.1%
 12084  F 3.73%

 It is pretty unbalanced. I think I need all GIDF residues in S35 families and 
 HCE residues in H families for training? That's means the E state will be the 
 least trained state...
 
--------------------------------------------------------------------------------
4. the transition table 
 --- getSSTransTable.cpp
 
 in the H families, using the Laplace's rule of succession,
 with the 7 states definitions: C G H I D E F
 The transition table:
 
 101780	11481	16  	1   	11762	2852	1	
 1   	1   	12298	1   	1   	1   	1	
 14  	1   	84927	12291	1   	1   	1	
 11819	1   	1   	1   	315 	159 	1	
 1   	1   	1   	1   	1   	10333	1744	
 2788	224 	1   	1   	1   	32399	10342	
 11491	595 	1   	1   	1   	1   	1	

--------------------------------------------------------------------------------
5. the relative solvent accessibility 

110.2 (Ala), 144.1 (Asp), 140.4 (Cys), 174.7 (Glu), 200.7 (Phe), 
 78.7 (Gly), 181.9 (His), 185.0 (Ile), 205.7 (Lys), 183.1 (Leu), 
200.1 (Met), 146.4 (Asn), 141.9 (Pro), 178.6 (Gln), 229.0 (Arg), 
117.2 (Ser), 138.7 (Thr), 153.7 (Val), 240.5 (Trp), 213.7 (Tyr)

ARNDC 106,248,157,163,135,
QEGHI 198,194, 84,184,196,
LKMFP 164,205,188,197,136,
STWYV 130,142,227,222,142

from Rost and Sander 1994

--------------------------------------------------------------------------------

5. get numbers of contact for each residue
 --- getContNum.cpp
 
 with CATH H domains,
 min: 0; max: 15; mean: 5.88596; sdv: 2.23433

 roughly normally distributed, with mode 5.
  
 
 save contact number distributions like
    0           1           2         ...
 G	0.00161438	0.0221214	0.0809808 ...
 K	0.00155342	0.0167368	0.0610343 ...
 E	0.000903148	0.0104507	0.0510064 ...
 ...
 in R:
 contDist=read.table('t.out', header=TRUE)
 hc <- hclust(dist(contDist))
 plot(hc)
 
 Two types of amino acids, 
 small and exposed: G,(ATH),(SNQPRKED)
 large and buried:  C,(LIWF),(VMY)
 
 The amino acids are classified the same way with KL-distances in the average
 tree. 
 
 Look at the cluster the other way: 
 in R:
 contDist=read.table('t.out', header=TRUE)
 hc <- hclust(dist(contDist))
 hc <- hclust(dist(contDist),method='average')
 plot(hc)
 
 Two groups, 4 sub groups. ((0,1,2,3,4),(5,6,7)),((8,9),(10,11,12,13))
 If we had the residues in four states... So what...
 
 However, if we look at the cluster tree generated using the  Kullback-
 Leibler divergence instead of the Euclidean distance, we might have a very 
 different conclusion. 
 
 mb=read.table('t.out'); # the 14 by 14, symmetric, zero-diagonal KL dis matrix
 b=as.dist(mb)
 plot(hclust(b,method='average'))
 plot(hclust(b))
 
 It should be a grouping of (0,1,2,3,4,5,6),(7,8),(9,10,11,12,13,14).
 
 7 and 8 contacts residues are 13 and 10%, together 23%. Not too few...
 
--------------------------------------------------------------------------------

6. To build a table for each domain. 
 ---getTable.cpp
 
domainName, length;  // the length of consensus sequence
AA ifIsInDomain, SS, #contacts, PSSM // for each consensus residue
AA ifIsInDomain, SS, #contacts, PSSM
...
--------------------------------------------------------------------------------
7 the distributions of secondary structure and burial states

 from the H familes, it can be see:
 
    exposed neutral buried
  C 101552  22221  8456
  G 7354    3108   1835
  H 50786   25684  20857
  I 8321    2480   1489
  D 7672    3080   1323
  E 24336   13217  8262
  F 6887    3535   1662
  
  const float STATE_DIS[21]={
    0.313319 ,0.0226893 ,0.15669  ,0.0256728 ,0.0236705 ,0.075084 ,0.0212485, 
    0.0685586,0.00958913,0.079243 ,0.00765156,0.00950274,0.0407785,0.0109066, 
    0.0260893,0.00566154,0.0643502,0.00459402,0.00408186,0.0254908,0.00512778};

  const float LOG_STATE_DIS[21]={
    -1.16053, -3.78586, -1.85348, -3.66232, -3.74353, -2.58915, -3.85147, 
    -2.68007, -4.64713, -2.53524, -4.87285, -4.65618, -3.1996 , -4.51839, 
    -3.64623, -5.17406, -2.74341, -5.383  , -5.5012 , -3.66944, -5.27308};
  
 Coil  (C) is way more likely to be exposed.
 Helix (H) is more often buried than other states.
 
 Burial states and secondary structures are not independent.
 
 The buried Strand begings (Db) is about 76 times less than exposed coil (Ce).
 That's a very unbalanced training set. 
 
 Even in the S35 familes, we have only 7735 Db and 7963 Ib. 
 Prior duplication is needed.  
 
 const size_t PRIOR_STATE_DUP[21]=
    {1,14,2,12,13,4,15,4,32,4,40,33,8,28,12,55,5,68,76,12,61};

--------------------------------------------------------------------------------
8. simple split of cath H in to three parts.
  // found 2467 H families from the S35 domains, not 2455...
  
    srandom(514);
    random();random();random();random();random();
    
    for (size_t i=0;i<2467;i++){
        cout<<random()%3<<'\n';
    }
 0: 3696 test.ls   the training set
 1: 3959 testtr.ls the test training set
 2: 3321 train.ls  the testing set

--------------------------------------------------------------------------------
9. train the neural networks
 --- trainNn.cpp

 redo the training program. Now it accept only wing_length, learning rate, 
 hidden units and random reeds. 
 
 1. WING_LENGTH      5    6    7    8
 2. learning rate    1e-3 1e-4 1e-5 1e-6
 3. hidden units     30   40   50   60
 4. random seeds     1    2    3    4
 
 so far with the online training, the best two are: 
 
 ==> 6_3_30_1.out <==
 5       9.15713e+06     1.50971e+06
 6       9.15734e+06     1.50971e+06
 
 ==> 7_3_50_3.out <==
 12      9.235e+06       1.52032e+06
 13      9.23528e+06     1.52033e+06
 
 in the second round trainning, the best two are:
 
 ==> 6_4_30_2.out <==
 21      9.58896e+06     1.54289e+06
 22      9.58792e+06     1.54289e+06
 
 ==> 6_4_60_2.out <==
 46      9.56094e+06     1.53746e+06
 47      9.56053e+06     1.53746e+06

 use the 6_4_60_2.Nn
 
--------------------------------------------------------------------------------
10. PCA of the input vectors
 --- pcaTest.cpp
 
 use the R nFactor package to determine the number of factors/components to
 retrain: 
 Acceleration Factor (n=22);
 Optimal Coordinates (n=22);
 
 OK, n is 22. But I still don't know the right number of clusters here...
 
--------------------------------------------------------------------------------
11. training of the network (again)
 --- trainNn.cpp
 --- testNn.cpp
 
 Found the simple prior duplication approach get us into trouble. 
 The most popular states, Ce, He and Cn are always under-predicted. 
 The log(Ppredict)-log(Freq) is alway negtive for the three.
 
 Try the snow ball approach. 
 
 within 200 iterations, all are under-trained...
 
 ==> 9_4_50_1.out <==
 198     1.06765e+06     1.17363e+06
 199     1.06764e+06     1.17362e+06
 ==> 9_4_30_4.out <==
 198     1.06772e+06     1.17335e+06
 199     1.06772e+06     1.17335e+06

 changed to 500 training iterations and re-run
 
 ==> 9_4_30_4.out <==
 498     1.06763e+06     1.17314e+06
 499     1.06763e+06     1.17314e+06
 ==> 9_4_50_2.out <==
 448     1.06745e+06     1.17333e+06
 449     1.06745e+06     1.17333e+06 
 ==> 9_4_40_4.out <==
 264     1.06764e+06     1.17336e+06
 265     1.06764e+06     1.17336e+06
 
 rerun with 
 1. WING_LENGTH      7    8    9    10
 2. learning rate    1e-2 1e-3 1e-4 1e-5
 3. hidden units     30   40   50   60
 4. random seeds     1    2    3    4

 ==> 10_4_40_2.out <==
 498     1.06743e+06     1.17319e+06
 499     1.06743e+06     1.17319e+06
 ==> 10_4_50_3.out <==
 498     1.06738e+06     1.1733e+06
 499     1.06738e+06     1.1733e+06

--------------------------------------------------------------------------------
12. the YASPIN algorithm
