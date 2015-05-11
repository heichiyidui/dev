################################################################################
# CGH NimbleGen Comparative genomic hybridization data                         #
################################################################################

mkdir cgh 

################################################################################
# Compared to the plot in the manual, the plot of our data is not very 
# encouraging at all...

# cgh/qc.png 

# get the file OID39087_Experimental_Report.txt from the dvd
# get files *_2800bp_avg_segMNT.txt from the dvd 
cp /var/run/media/kuang/39087_2of2/OID39087_Experimental_Report.txt cgh/
cp /var/run/media/kuang/39087_2of2/Processed_Data_TXT_Files/*_2800bp*.zip cgh/

cd cgh/ 
for zfile in *.zip ; do unzip $zfile ; done
cd ..
rm cgh/*.zip 


tail -n +2 cgh/OID39087_Experimental_Report.txt | awk '{print $4,$8}' | \
  awk '{print "mv cgh/" $1 "_2800bp_avg_segMNT.txt cgh/" $2 "_seg.txt"} ' > t.sh
source t.sh

for ifile in cgh/MF*.txt ; do 
    tail -n +2 $ifile | grep -v "chrY" | grep -v "chrX" | \
    awk '{print $1,$2,$5}' | sed 's/chr//' > t.out 
    mv t.out $ifile
done 

####################################
# log2 ratio of expression levels comparison 
# cgh/lr_comp.png 
# can find deletions in one or two persons.

# log2 ratio :
# size:      6529040
# min:       -4.4507
# max:       3.8796
# mean:      0.0070
# std:       0.1203
# quartiles: -4.4507 -0.0589 0.0018 0.0658 3.8796
# 8 tiles:   -4.4507 -0.1050 -0.0589 -0.0266 0.0018 0.0308 0.0658 0.1213 3.8796

# mostly between -+1

for ifile in cgh/MF_* ; do awk '{if ($3>0.8) print $0}' $ifile | wc ; done 
#    156     468    3048
#    120     360    2352
#     91     273    1788
#    120     360    2346
#    116     348    2264
#    164     492    3222
#    190     570    3761
#    201     603    3957

for ifile in cgh/MF_* ; do awk '{if ($3<-0.8) print $0}' $ifile | wc ; done 
#    285     855    5853
#    307     921    6380
#    545    1635   11226
#    566    1698   11671
#    328     984    6772
#    276     828    5704
#    366    1098    7551
#    375    1125    7704

# more deletions than duplications 
################################################################################

