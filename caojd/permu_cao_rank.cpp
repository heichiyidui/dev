#include <cstdlib>
#include <cassert>
#include <cstring>
#include <string>
#include <vector>
#include <map>
#include <set>
#include <cmath>
#include <algorithm>
#include <iterator>
#include <iostream>
#include <fstream>
#include <sstream>
#include "../util/mathutil.hpp"

int main(int argc,char**argv){
    using namespace std;
    
    if (argc<2){
        cerr<<"Give me some domain id list!\n";
        exit(1);
    }
    ifstream id_file(argv[1]);
    if (!id_file.is_open()){
        cerr<<"Failed to open the domain id file: "<<argv[1]<<'\n';
        exit(1);
    }
    string ts;
    vector <string> dom_ids;
    for (;id_file>>ts;)
        dom_ids.push_back(ts);
    id_file.close();
    
    srand(514);
    
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
    
    float cao[400][400];
    ifstream cfile("../cath/cao/cao150");
    for (size_t i=0;i<400;i++)
        for (size_t j=0;j<400;j++)
            cfile>>cao[i][j];
    cfile.close();
    
    const float GAP_CAO  = -3.5;
    
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
        
        // well, gaps on the edges don't count
        for (size_t i=0; i<aln_lines.size();i++){
            ts=aln_lines[i];
            ts.erase(ts.find_last_not_of("-")+1);
            ts.append(aln_lines[i].length()-ts.length(),'u');
            
            size_t leading_gap_length=ts.find_first_not_of("-");
            ts.erase(0,leading_gap_length);
            ts.insert(0,leading_gap_length,'u');
            
            aln_lines[i]=ts;
        }
        
        // transpose alignment 
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
    ifstream cont_file("../cath/index/cath_s35.condef");
    assert(cont_file);
    for (;getline(cont_file,ts);){
        string cont_line ;
        getline(cont_file,cont_line);
        string dom_id=ts.substr(1,7);
        
        if (find(dom_ids.begin(),dom_ids.end(),dom_id)==dom_ids.end())
            continue ;
        
        vector <vector<bool> > cont_map;
        size_t dom_len=dom_alns[dom_id].size();
        for (size_t i=0;i<dom_len;i++){
            vector <bool> map_row(dom_len,false);
            cont_map.push_back(map_row);
        }
        
        replace(cont_line.begin(),cont_line.end(),'-',' ');
        stringstream tss(cont_line);
        size_t c_i,c_j;
        for (;tss>>c_i>>c_j;){
            if (c_j - c_i < 5)
                continue;
            cont_map[c_i][c_j]=true;
            cont_map[c_j][c_i]=true;
        }
        dom_cont_maps[dom_id]=cont_map; 
    }
    cont_file.close();
    
    ////////////////////////////////////////////////////////////////////////////
    // 3 calculate permutation contact scores                                 //
    ////////////////////////////////////////////////////////////////////////////
    for (size_t l =0; l<dom_ids.size();l++){
        string dom_id = dom_ids[l];
        
        vector <vector <int> > aln=dom_alns[dom_id];
        size_t dom_len=aln.size();
        size_t aln_len=aln[0].size();
        size_t mat_len=(dom_len-8)*(dom_len-9)/2;
        
        // 3.1 calculate the original all-against-all CAO score
        vector <float> org_cao_scores;
        vector <int> cont_map;
        for (size_t i=4;i<dom_len;i++){
            for (size_t j=i+5;j<dom_len;j++){
                float cao_score=0;
                for (size_t m=0;m<aln_len;m++){
                    for (size_t n=m+1;n<aln_len;n++){
                        // GAP 
                        if (aln[i][m]==20 || aln[j][m]==20)
                            {cao_score += GAP_CAO; continue;}
                        if (aln[i][n]==20 || aln[j][n]==20)
                            {cao_score += GAP_CAO; continue;}
                        // termini GAP 'u'
                        if (aln[i][m]==21 || aln[j][m]==21)
                            continue;
                        if (aln[i][n]==21 || aln[j][n]==21)
                            continue;
                        int c_m=aln[i][m]*20+aln[j][m];
                        int c_n=aln[i][n]*20+aln[j][n];
                        cao_score+=cao[c_m][c_n];
                    }
                }
                org_cao_scores.push_back(cao_score);
                cont_map.push_back(dom_cont_maps[dom_id][i][j]);
            }
        }
        
        // 3.2 permute the alignment 100 times 
        vector <float> per_cao_scores[mat_len];
        for (size_t iter=0;iter<100;iter++){
            for (size_t i=0;i<dom_len;i++)
                random_shuffle(aln[i].begin(),aln[i].end());
            
            size_t mat_index=0;
            for (size_t i=4;i<dom_len;i++){
                for (size_t j=i+5;j<dom_len;j++){
                    float per_cao_score=0;
                    
                    for (size_t m=0;m<aln_len;m++){
                        for (size_t n=m+1;n<aln_len;n++){
                            // GAP 
                            if (aln[i][m]==20 || aln[j][m]==20)
                                {per_cao_score += GAP_CAO; continue;}
                            if (aln[i][n]==20 || aln[j][n]==20)
                                {per_cao_score += GAP_CAO; continue;}
                            // termini GAP 'u'
                            if (aln[i][m]==21 || aln[j][m]==21)
                                continue;
                            if (aln[i][n]==21 || aln[j][n]==21)
                                continue;
                            int c_m=aln[i][m]*20+aln[j][m];
                            int c_n=aln[i][n]*20+aln[j][n];
                            per_cao_score+=cao[c_m][c_n];
                        }
                    }
                    per_cao_scores[mat_index].push_back(per_cao_score);
                    mat_index++;
                }
            }
        }
        
        // 3.3 get the rank-normal of org_cao_score
        for (size_t i=0;i<mat_len ; i++){
            float org_score=org_cao_scores[i];
            vector <float> per_scores=per_cao_scores[i];
            per_scores.push_back(org_score);
            sort(per_scores.begin(),per_scores.end());
            
            int rank = find(per_scores.begin(),\
                            per_scores.end(),org_score) \
                          - per_scores.begin();
            cout<<-norm_isf((rank+1.0)/(101+2.0))<<' '<<cont_map[i]<<'\n';
        }
        
        // break;             
    }
    return(0);
}

