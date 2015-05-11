#!/bin/sh
#$-S /bin/sh
#$ -M kuang.lin@kcl.ac.uk
#$ -m e
#$ -N cv_train_expo
#$ -cwd
#$ -q long.q,bignode.q

hu=$1
rs=$2

./cv_expos_train.py -hu $hu -rs $rs  > cv_exp_tra_${hu}_${rs}.out


