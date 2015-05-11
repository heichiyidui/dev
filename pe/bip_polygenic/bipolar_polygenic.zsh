# 1. Bipolar data set Nat Genet. 2011 Sep 18;43(10):977-83

wget http://www.med.unc.edu/pgc/files/resultfiles/pgc.bip.2012-04.zip

unzip pgc.bip.2012-04.zip

pgc.bip.full.2012-04.txt
pgc.bip.clump.2012-04.txt
# 2427220 SNPs in the full set
# 108834  SNPs in the clumped set 

# 2. pe set 
pe19 

awk '{print $2}' pe19.bim > t.ls
grab -f t.ls pgc.bip.clump.2012-04.txt | wc
# 48878/108834 
# many SNPs are missing ... only about 45% are available.

pe20_poly
# 108353/108834 
# The imputed set is way more complete. 
# Most SNPs are available.

# 3. scores
# P thresholds: 1, 0.5, 0.4, 0.3, 0.2, 0.1, 0.05, 0.01

tail -n +2 pgc.bip.clump.2012-04.txt \
    | awk '{if ($8<1.00) print $1,$4,log($6)}' > t.profile;
plink --noweb --bfile pe20_poly --score t.profile; 
cp plink.profile pe20_bip_100.profile

tail -n +2 pgc.bip.clump.2012-04.txt |\
 awk '{if ($8<0.50) print $1,$4,log($6)}'> t.profile;\
 plink --noweb --bfile pe20_poly --score t.profile; \
 cp plink.profile pe20_bip_050.profile
tail -n +2 pgc.bip.clump.2012-04.txt | awk '{if ($8<0.40) print $1,$4,log($6)}'> t.profile;plink --noweb --bfile pe20_poly --score t.profile; cp plink.profile pe20_bip_040.profile
tail -n +2 pgc.bip.clump.2012-04.txt | awk '{if ($8<0.30) print $1,$4,log($6)}'> t.profile;plink --noweb --bfile pe20_poly --score t.profile; cp plink.profile pe20_bip_030.profile
tail -n +2 pgc.bip.clump.2012-04.txt | awk '{if ($8<0.20) print $1,$4,log($6)}'> t.profile;plink --noweb --bfile pe20_poly --score t.profile; cp plink.profile pe20_bip_020.profile
tail -n +2 pgc.bip.clump.2012-04.txt | awk '{if ($8<0.10) print $1,$4,log($6)}'> t.profile;plink --noweb --bfile pe20_poly --score t.profile; cp plink.profile pe20_bip_010.profile
tail -n +2 pgc.bip.clump.2012-04.txt | awk '{if ($8<0.05) print $1,$4,log($6)}'> t.profile;plink --noweb --bfile pe20_poly --score t.profile; cp plink.profile pe20_bip_005.profile
tail -n +2 pgc.bip.clump.2012-04.txt | awk '{if ($8<0.01) print $1,$4,log($6)}'> t.profile;plink --noweb --bfile pe20_poly --score t.profile; cp plink.profile pe20_bip_001.profile



