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
    
    ifstream id_file("index/h.ls");
    string ts;
    vector <string> dom_ids;
    for (;id_file>>ts;) dom_ids.push_back(ts);
    id_file.close();
    
    set <string> dom_id_set(dom_ids.begin(),dom_ids.end());
    
    ////////////////////////////////////////////////////////////////////////////
    // 1. read the VTML and CAO matrices                                      //
    ////////////////////////////////////////////////////////////////////////////
    
    map <char,size_t> AA_CODE_TO_INT;
    AA_CODE_TO_INT['A']=0 ; AA_CODE_TO_INT['R']=1 ; AA_CODE_TO_INT['N']=2 ; 
    AA_CODE_TO_INT['D']=3 ; AA_CODE_TO_INT['C']=4 ; AA_CODE_TO_INT['Q']=5 ; 
    AA_CODE_TO_INT['E']=6 ; AA_CODE_TO_INT['G']=7 ; AA_CODE_TO_INT['H']=8 ; 
    AA_CODE_TO_INT['I']=9 ; AA_CODE_TO_INT['L']=10; AA_CODE_TO_INT['K']=11; 
    AA_CODE_TO_INT['M']=12; AA_CODE_TO_INT['F']=13; AA_CODE_TO_INT['P']=14; 
    AA_CODE_TO_INT['S']=15; AA_CODE_TO_INT['T']=16; AA_CODE_TO_INT['W']=17;
    AA_CODE_TO_INT['Y']=18; AA_CODE_TO_INT['V']=19; AA_CODE_TO_INT['-']=20;
    AA_CODE_TO_INT['u']=21;
    
    // const float GAP_VTML = -1.3;
    // const float END_VTML = -0.1;
    const float GAP_VTML = atof(argv[1]);
    const float END_VTML = atof(argv[2]);
    
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
    
    ////////////////////////////////////////////////////////////////////////////
    // 2. read alignment and contact map                                      //
    ////////////////////////////////////////////////////////////////////////////
    map <string , vector <vector <int> > > dom_alns;
    for (size_t dom_index=0;dom_index<dom_ids.size();dom_index++){
        string dom_id = dom_ids[dom_index];
        
        size_t aln_depth, dom_len;
        ifstream aln_file(("seq_aln/"+dom_id).c_str());
        assert(aln_file);
        aln_file>>dom_len>>aln_depth;
        vector < vector < int> > dom_aln(dom_len);
        for (size_t i=0;i<dom_len;i++){
            aln_file>>ts;
            for (size_t j=0;j<aln_depth;j++)
                dom_aln[i].push_back(AA_CODE_TO_INT[ts[j]]);
        }
        
        dom_alns[dom_id]=dom_aln;
    }
    
    map <string, vector <vector<bool> >  > dom_cont_maps;
    ifstream cont_file("index/cath_s35.condef");
    assert(cont_file);
    
    size_t dom_len,cont_num;
    string cont_line ;
    for (;cont_file>>ts>>dom_len>>cont_num;){
        getline(cont_file,cont_line);
        getline(cont_file,cont_line);
        string dom_id=ts.substr(1,7);
        
        if (dom_id_set.count(dom_id)==0)
            continue ;
        
        vector <vector<bool> > cont_map;
        for (size_t i=0;i<dom_len;i++)
            cont_map.push_back(vector <bool> (dom_len,false));
        
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
    vector <float> all_cont_scores;
    
    for (size_t dom_index=0;dom_index<dom_ids.size();dom_index++){
        string dom_id = dom_ids[dom_index];
        
        vector <vector <int> > aln=dom_alns[dom_id];
        size_t dom_len = aln.size();
        size_t aln_len = aln[0].size();
        size_t mat_len = (dom_len-4)*(dom_len-5)/2;
        
        vector <bool> cont_map;
        for (size_t i=0;i<dom_len;i++)
            for (size_t j=i+5;j<dom_len;j++)
                cont_map.push_back(dom_cont_maps[dom_id][i][j]);
        
        ///////////////////////////////
        // 3.1 calculate the residue conservation score products
        
        vector <float> res_cons_scores(dom_len);
        for (size_t i=0;i<dom_len;i++){
            float res_cons_score=0;
            for (size_t m=0;m<aln_len;m++){
                for (size_t n=m+1;n<aln_len;n++){
                    res_cons_score += vtml[aln[i][m]][aln[i][n]];
                }
            }
            res_cons_scores[i]=res_cons_score;
        }
        vector <float> res_score_prods;
        for (size_t i=0;i<dom_len;i++)
            for (size_t j=i+5;j<dom_len;j++)
                res_score_prods.push_back( res_cons_scores[i] \
                                         * res_cons_scores[j]);
        
        ///////////////////////////////
        // 3.3 sort and normalize score 
        
        vector <float> sorted_res_score_prods (res_score_prods);
        sort(sorted_res_score_prods.begin(),sorted_res_score_prods.end());
        
        for (size_t i=0;i<mat_len ; i++){
            if (cont_map[i]==0 )
                continue;
            
            float norm_score1;
            
            float score ; int rank;
            
            score=res_score_prods[i];
            rank = find(sorted_res_score_prods.begin(),\
                        sorted_res_score_prods.end(),score) \
                      - sorted_res_score_prods.begin();
            norm_score1= -norm_isf((rank+1.0)/(mat_len+2.0));
            
            //  cout<<norm_score1<<' '<<cont_map[i]<<'\n';
            all_cont_scores.push_back(norm_score1);
        }
    }
    
    double mean = accumulate(all_cont_scores.begin(), \
                             all_cont_scores.end(), 0.0) \
                           / all_cont_scores.size();
    cout<<mean<<'\n';
    
    return(0);
}

