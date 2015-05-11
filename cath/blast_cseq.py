#!/home/klinbrc/bin/python3
import os
import sys
import subprocess

################################################################################
# global parameters                                                            #
################################################################################

blast='~/bin/blastp ' 
dbcmd='~/bin/blastdbcmd '

db=' -db ../nr/nr '

dbtype=' -dbtype prot '
ofmt=' -outfmt \"7 sseqid sstart send slen pident\" '
# table: subject seq id, subject seq align start, subject seq align end, 
# subject seq length, percent identity of align 

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

for dom_id in domain_ids:
    if not os.path.isfile('./cseq/'+dom_id):
        sys.stderr.write('Failed to open the query seq file cseq/'+dom_id+'\n')
        exit(1)
    
    ###################################
    # 3.1 blast it 
    print(dom_id,'Blast searching...')
    subprocess.call(blast+db+ofmt+\
                    ' -query cseq/'+dom_id+\
                    ' -out bhits/'+dom_id,shell=True)
    
    try:
        ifile=open('bhits/'+dom_id)
        hits_lines=ifile.read().split('\n')[:-1]
        ifile.close()
    except:
        sys.stderr.write('Failed to open BLAST hits file bhits/'+dom_id+'\n')
        sys.stderr.write('Exiting\n')
        exit(2)
    
    ###################################
    # 3.2 get subject sequences
    print('Retrieving the Blast hits...')
    
    # read the cseq first 
    cseq=''.join(open('cseq/'+dom_id).read().split('\n')[1:])
    
    hits=''
    hits_seq=set()
    for line in hits_lines:
        if line.startswith('#'):
            continue;
        
        # if float(line.split('\t')[4])<50 : # alignment sequence identity
        #     continue
        
        s_id    = line.split('|')[1]
        s_begin = int(line.split()[1])
        s_end   = int(line.split()[2])
        
        if s_end - s_begin + 1 < len(cseq) * 0.5 :
            continue

        cmd_line=dbcmd+db+'-entry '+s_id \
                 +' -range '+str(s_begin)+'-'+str(s_end)
        
        s_seq=subprocess.getoutput(cmd_line).split('\n')
        if not s_seq[0].startswith('>'):
            sys.stderr.write('Error in getting subject sequence\n')
            sys.stderr.write(cmd_line)
        
        s_seq[0]=s_seq[0][:80]
        # cut the title line because very long title might be trouble to align 

        if ''.join(s_seq[1:]) in hits_seq: # the sequence is included already 
            continue

        hits+='\n'.join(s_seq)+'\n'
        hits_seq.add(''.join(s_seq[1:]))
    
    # write the sequences to be aligned.
    ofile=open('c_align_in/'+dom_id,'w')
    ofile.write(open('cseq/'+dom_id).read())
    if hits=='': # no Blast hits at all
        ofile.write(open('cseq/'+dom_id).read())
    
    ofile.write(hits)
    ofile.close()
