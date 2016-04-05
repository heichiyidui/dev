#include <cstdlib>
#include <cassert>
#include <cstring>
#include <string>
#include <vector>
#include <map>
#include <set>
#include <algorithm>
#include <iostream>
#include <fstream>
#include <sstream>
#include "../util/mathutil.hpp"

int main(int argc,char**argv){
    using namespace std;
    srandom(514);
    
    if (argc<2){cerr<<"Give me some domain id file!\n"; exit(1); }
    
    ifstream id_file(argv[1]);
    if (!id_file.is_open()){
        cerr<<"Failed to open the domain id file: "<<argv[1]<<'\n';
        exit(1);
    }
    string ts;
    vector <string> dom_ids;
    for (;id_file>>ts;) dom_ids.push_back(ts);
    
    id_file.close();
    set <string> dom_id_set(dom_ids.begin(),dom_ids.end());
    
    ////////////////////////////////////////////////////////////////////////////
    // 1. read the VTML and CAO matrices                                      //
    ////////////////////////////////////////////////////////////////////////////
    
    char aa_set[]={'A','R','N','D','C','Q','E','G','H','I',\
                   'L','K','M','F','P','S','T','W','Y','V'};
    const set <char> AA_SET(aa_set,aa_set+sizeof(aa_set)/sizeof(aa_set[0]));
    
    map <char,size_t> AA_CODE_TO_INT;
    AA_CODE_TO_INT['A']=0 ; AA_CODE_TO_INT['R']=1 ; AA_CODE_TO_INT['N']=2 ; 
    AA_CODE_TO_INT['D']=3 ; AA_CODE_TO_INT['C']=4 ; AA_CODE_TO_INT['Q']=5 ; 
    AA_CODE_TO_INT['E']=6 ; AA_CODE_TO_INT['G']=7 ; AA_CODE_TO_INT['H']=8 ; 
    AA_CODE_TO_INT['I']=9 ; AA_CODE_TO_INT['L']=10; AA_CODE_TO_INT['K']=11; 
    AA_CODE_TO_INT['M']=12; AA_CODE_TO_INT['F']=13; AA_CODE_TO_INT['P']=14; 
    AA_CODE_TO_INT['S']=15; AA_CODE_TO_INT['T']=16; AA_CODE_TO_INT['W']=17;
    AA_CODE_TO_INT['Y']=18; AA_CODE_TO_INT['V']=19; AA_CODE_TO_INT['-']=20;
    AA_CODE_TO_INT['u']=21;
    
    const float GAP_VTML =  -1.5 ;
    const float END_VTML =  -0.3 ;
    
    float vtml[22][22];
    ifstream vfile("../cath/vtml/vtml100");
    for (size_t i=0;i<20;i++)
        for (size_t j=0;j<20;j++)
            vfile>>vtml[i][j];
    vfile.close();
    
    for (size_t i=0;i<22;i++)
        for (size_t j=0;j<22;j++){
            if (i==20 || j==20) // gaps
                vtml[i][j] = GAP_VTML;
            if (i==21 || j==21) // termini gaps
                vtml[i][j] = END_VTML;
        }
    
    // const float GAP_CAO  = atof(argv[2]);
    // const float END_CAO  = atof(argv[3]);
    const float GAP_CAO  = -2.5;
    const float END_CAO  = -0.3;
    
    float cao[484][484];
    for (size_t i=0;i<484;i++)
        for (size_t j=0;j<484;j++)
            cao[i][j]=-10.0;
    
    ifstream cfile("../cath/cao/cao150");
    for (size_t i=0;i<20;i++)
        for (size_t j=0;j<20;j++)
            for (size_t m=0;m<20;m++)
                for (size_t n=0;n<20;n++)
                    cfile>>cao[i*22+j][m*22+n];
    cfile.close();
    
    for (size_t i=0;i<22;i++)
        for (size_t j=0;j<22;j++)
            for (size_t m=0;m<22;m++)
                for (size_t n=0; n<22; n++){
                    if (i==20 || j==20 || m==20 || n==20) // gaps
                        cao[i*22+j][m*22+n] = GAP_CAO ;
                    if (i==21 || j==21 || m==21 || n==21) // termini gaps
                        cao[i*22+j][m*22+n] = END_CAO ;
                }
    
    
    ////////////////////////////////////////////////////////////////////////////
    // 2. read alignment and contact map                                      //
    ////////////////////////////////////////////////////////////////////////////
    map <string , vector <vector <int> > > dom_alns;
    for (size_t dom_index=0;dom_index<dom_ids.size();dom_index++){
        string dom_id = dom_ids[dom_index];
        
        vector <string> aln_lines;
        ifstream aln_file(("../cath/seq_aln/"+dom_id).c_str());
        assert(aln_file);
        
        while (!aln_file.eof()){
            getline(aln_file,ts);
            if (ts.size()>0)
                aln_lines.push_back(ts);
        }
        aln_file.close();
        
        // well, gaps on the edges are different
        for (size_t i=0; i<aln_lines.size();i++){
            ts=aln_lines[i];
            ts.erase(ts.find_last_not_of("-")+1);
            ts.append(aln_lines[i].length()-ts.length(),'u');
            
            size_t leading_gap_length=ts.find_first_not_of("-");
            ts.erase(0,leading_gap_length);
            ts.insert(0,leading_gap_length,'u');
            
            aln_lines[i]=ts;
        }
        
        // transpose the alignment 
        vector <vector <int> > aln_cols;
        for (size_t i=0;i<aln_lines[0].length();i++){
            if (AA_SET.count(aln_lines[0][i])==0) continue;
            
            vector <int> aln_col;
            for (size_t j=0;j<aln_lines.size();j++)
                aln_col.push_back(AA_CODE_TO_INT[aln_lines[j][i]]);
            
            aln_cols.push_back(aln_col);
        }
        dom_alns[dom_id]=aln_cols;
    }
    
    map <string, vector <vector<bool> >  > dom_cont_maps;
    ifstream cont_file("index/cath_s35.condef");
    assert(cont_file);
    for (;getline(cont_file,ts);){
        string cont_line ;
        getline(cont_file,cont_line);
        string dom_id=ts.substr(1,7);
        
        if (dom_id_set.count(dom_id)==0)
            continue ;
        
        size_t dom_len , cont_num;
        stringstream tss(ts);
        string ts2;
        tss>>ts2>>dom_len>>cont_num;
        
        vector <vector<bool> > cont_map;
        for (size_t i=0;i<dom_len;i++){
            vector <bool> map_row(dom_len,false);
            cont_map.push_back(map_row);
        }
        
        stringstream tss2(cont_line);
        size_t c_i,c_j;
        for (size_t i=0;i<cont_num;i++){
            tss2>>c_i>>c_j;
            cont_map[c_i][c_j]=true;
            cont_map[c_j][c_i]=true;
        }
        dom_cont_maps[dom_id]=cont_map; 
    }
    cont_file.close();
    
    ////////////////////////////////////////////////////////////////////////////
    // 3 calculate domain alignment feature scores                            //
    ////////////////////////////////////////////////////////////////////////////
    for (size_t l =0; l<dom_ids.size();l++){
        string dom_id = dom_ids[l];
        
        vector <vector <int> > aln=dom_alns[dom_id];
        size_t dom_len=aln.size();
        size_t aln_len=aln[0].size();
        size_t mat_len=(dom_len-8)*(dom_len-9)/2;
        
        vector <int> cont_map;
        for (size_t i=4;i<dom_len;i++)
            for (size_t j=i+5;j<dom_len;j++)
                cont_map.push_back(dom_cont_maps[dom_id][i][j]);
        
        ///////////////////////////////
        // 3.1 calculate the original all-against-all CAO scores
        vector <float> org_cao_scores;
        for (size_t i=4;i<dom_len;i++){
            for (size_t j=i+5;j<dom_len;j++){
                float cao_score=0;
                for (size_t m=0;m<aln_len;m++){
                    for (size_t n=m+1;n<aln_len;n++){
                        cao_score+=cao[aln[i][m]*22+aln[j][m]] \
                                      [aln[i][n]*22+aln[j][n]];
                    }
                }
                org_cao_scores.push_back(cao_score);
            }
        }
    
        ///////////////////////////////
        // 3.2 permute the alignment 100 times 
        // get the z-score of org_cao_score and the std_dev 
        vector <float> z_scores(mat_len) ;
        vector <float> std_devs(mat_len) ;
        
        size_t mat_index=0;
        for (size_t i=4;i<dom_len;i++){
            for (size_t j=i+5;j<dom_len;j++){
                float org_cao_score = org_cao_scores[mat_index];
                
                vector <float> per_cao_scores(100);
                for (size_t iter=0;iter<100;iter++){
                    float per_cao_score=0;
                    
                    random_shuffle(aln[j].begin(),aln[j].end());
                    
                    for (size_t m=0;m<aln_len;m++){
                        for (size_t n=m+1;n<aln_len;n++){
                            per_cao_score+=cao[aln[i][m]*22+aln[j][m]] \
                                              [aln[i][n]*22+aln[j][n]];
                        }
                    }
                    per_cao_scores[iter]=per_cao_score;
                }
                
                double mean = accumulate(per_cao_scores.begin(),   \
                                         per_cao_scores.end(),0.0) \
                            / per_cao_scores.size();
                
                vector<double> diff(per_cao_scores.size());
                transform(per_cao_scores.begin(), per_cao_scores.end(), \
                          diff.begin(), bind2nd(minus<double>(), mean));
                
                double sq_sum = inner_product(diff.begin(), diff.end(), \
                                              diff.begin(), 0.0);
                
                double stdev = sqrt(sq_sum / per_cao_scores.size());
                
                if (stdev>0.0){
                    z_scores[mat_index] = (org_cao_score - mean) / stdev;
                }else{
                    z_scores[mat_index] = 0.0;
                }
                std_devs[mat_index] = stdev ;
                
                mat_index+=1;
            }
        }
        ///////////////////////////////
        // 3.3 sort and output 
        
        vector <float> sorted_z_scores(z_scores);
        sort(sorted_z_scores.begin(),sorted_z_scores.end());
        
        vector <float> sorted_std_devs(std_devs);
        sort(sorted_std_devs.begin(),sorted_std_devs.end());
        
        for (size_t i=0;i<mat_len ; i++){
            if (cont_map[i]==0 && rand()%50 != 0 )
                continue;
            
            float norm_score3,norm_score4;
            
            float score ; int rank;
            
            score=z_scores[i];
            rank = find(sorted_z_scores.begin(),\
                        sorted_z_scores.end(),score) \
                      - sorted_z_scores.begin();
            norm_score3= -norm_isf((rank+1.0)/(mat_len+2.0));
            
            score=std_devs[i];
            rank = find(sorted_std_devs.begin(),\
                        sorted_std_devs.end(),score) \
                      - sorted_std_devs.begin();
            norm_score4= -norm_isf((rank+1.0)/(mat_len+2.0));
            
            cout<<norm_score3<<' '<<norm_score4<<' '<<cont_map[i]<<'\n';
        }
        // break;             
    }
    return(0);
}
        
