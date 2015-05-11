#!/home/klinbrc/bin/python3
import os
import sys
import subprocess

################################################################################
# global parameters                                                            #
################################################################################

blast='~/bin/blastp ' 
dbcmd='~/bin/blastdbcmd '

nrdb=' -db ../nr/nr '
nrfilt_db=' -db ../nr/nrfilt '

dbtype=' -dbtype prot '
ofmt=' -outfmt \"7 sseqid sstart send slen pident\" '
# table: subject seq id, subject seq align start, subject seq align end, 
# subject seq length, percent identity of align 

WIN_LENGTH=9 
# try to add extra AAs to the termini in extracting subject sequences

################################################################################
# 1. make the output directions                                                #
################################################################################

# the hits of the initial PsiBlast search
if not os.path.exists('./bhits'):
    subprocess.call('mkdir ./bhits',shell=True)

# the query and subject sequences to be aligned
if not os.path.exists('./c_align_in'):
    subprocess.call('mkdir ./c_align_in',shell=True)

################################################################################
# 2. read the list of query sequence ids                                       #
################################################################################

ifile=open(sys.argv[1])
domain_ids=[]
for line in ifile:
    domain_ids.append(line.split()[0])

ifile.close()

################################################################################
# 3. PsiBlast and align query sequences                                        #
################################################################################

for id in domain_ids:
    if not os.path.isfile('./seq/'+id):
        sys.stderr.write('Failed to open the query sequence file seq/'+id+'\n')
        exit(1)
    
    ###################################
    # 3.1 blast it 
    print(id,'PsiBlast searching...')
    subprocess.call(blast+nrfilt_db+ofmt+\
                    ' -query seq/'+id+\
                    ' -out bhits/'+id,shell=True)
    
    try:
        ifile=open('bhits/'+id)
        hits_lines=ifile.read().split('\n')[:-1]
        ifile.close()
    except:
        sys.stderr.write('Failed to open PsiBlast output file bhits/'+id+'\n')
        sys.stderr.write('Exiting\n')
        exit(2)
    
    ###################################
    # 3.2 get subject sequences
    print('Retrieving the Blast hits...')
    hits=''
    for line in hits_lines:
        if line.startswith('#'):
            continue;
        
        if float(line.split('\t')[4])<50 : # alignment sequence identity
            continue
        
        s_id=line.split('|')[1]
        
        s_begin=int(line.split()[1])-WIN_LENGTH
        s_begin=max( (1, s_begin) )
        
        s_end=int(line.split()[2])+WIN_LENGTH
        s_end=min( (s_end, int(line.split()[3]) ) )
        
        cmd_line=dbcmd+nrdb+'-entry '+s_id \
                 +' -range '+str(s_begin)+'-'+str(s_end)
        # after formatting nrfilt.fasta, got trouble getting sequences of some 
        # ids such as 449265696. 
        # using the downloaded nr database instead. 
        
        s_seq=subprocess.getoutput(cmd_line).split('\n')
        if not s_seq[0].startswith('>'):
            sys.stderr.write('Error in getting subject sequence\n')
            sys.stderr.write(cmd_line)
        
        s_seq[0]=s_seq[0][:80]
        # cut the title line because very long title might be trouble to align 
        hits+='\n'.join(s_seq)+'\n'
    
    # write the sequences to be aligned.
    ofile=open('c_align_in/'+id,'w')
    
    ifile=open('seq/'+id)
    ofile.write(ifile.read())
    ifile.close()
    
    if hits=='': # no Blast hits at all
        ifile=open('seq/'+id)
        hits+=ifile.read()
        ifile.close()
    
    ofile.write(hits)
    ofile.close()
    

