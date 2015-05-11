################################################################################
#!/bin/sh
#$-S /bin/sh
# Deliver an email when the job ends
#$ -M kuang.lin@kcl.ac.uk
#$ -m beas
#$ -N uz_fastq
#use current directory as working directory
#$ -cwd
# -o /home/klinbrc/scratch/dev/wgs
# -e /home/klinbrc/scratch/dev/wgs

# Choose the queues 
# -q 78gb.q,54gb.q
# Declare how much memory is required PER slot - default is 2Gbytes
# -l h_vmem=8G
# Declare how many slots/jobs are required - if not busy 100
# -pe mpi 100
################################################################################
