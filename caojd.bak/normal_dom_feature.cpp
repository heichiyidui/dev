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
        cerr<<"Give me some domain feature file!\n";
        exit(1);
    }
    ifstream df_file(argv[1]);
    if (!df_file.is_open()){
        cerr<<"Failed to open the domain feature file: "<<argv[1]<<'\n';
        exit(1);
    }
    ////////////////////////////////////////////////////////////////////////////
    // read domain id and length
    
    string ts;
    for (;df_file>>ts;){
        string dom_id=ts.substr(1,7);
        size_t dom_len; 
        df_file>>dom_len;
        size_t mat_len=(dom_len-8)*(dom_len-9)/2;
        vector <double> res_cons_prds(mat_len);
        vector <double> cao_scores(mat_len);
        vector <double> cao_z_scores(mat_len);
        vector <double> std_devs(mat_len);
        vector <int>    cont_map(mat_len);
        
        double td1,td2,td3,td4;
        int tm1;
        for (size_t i=0;i<mat_len;i++){
            df_file>>td1>>td2>>td3>>td4>>tm1;
            res_cons_prds[i] = td1;
            cao_scores   [i] = td2;
            cao_z_scores [i] = td3;
            std_devs     [i] = td4;
            cont_map     [i] = tm1;
        }
        vector <double> sorted_res_cons_prds = res_cons_prds ;
        vector <double> sorted_cao_scores    = cao_scores    ;
        vector <double> sorted_cao_z_scores  = cao_z_scores  ;
        vector <double> sorted_std_devs      = std_devs      ;
        
        sort(sorted_res_cons_prds.begin(), sorted_res_cons_prds.end());
        sort(sorted_cao_scores   .begin(), sorted_cao_scores   .end());
        sort(sorted_cao_z_scores .begin(), sorted_cao_z_scores .end());
        sort(sorted_std_devs     .begin(), sorted_std_devs     .end());
        
        cout<<'>'<<dom_id<<' '<<dom_len<<'\n';
        for (size_t i=0;i<mat_len;i++){
            float rank1 = find(sorted_res_cons_prds.begin(),\
                               sorted_res_cons_prds.end(),res_cons_prds[i]) \
                       - sorted_res_cons_prds.begin();
            rank1=-norm_isf((rank1+1.0)/(mat_len+2.0));
            
            float rank2 = find(sorted_cao_scores.begin(),\
                               sorted_cao_scores.end(),cao_scores[i]) \
                       - sorted_cao_scores.begin();
            rank2=-norm_isf((rank2+1.0)/(mat_len+2.0));
            
            float rank3 = find(sorted_cao_z_scores.begin(),\
                               sorted_cao_z_scores.end(),cao_z_scores[i]) \
                       - sorted_cao_z_scores.begin();
            rank3=-norm_isf((rank3+1.0)/(mat_len+2.0));
            
            float rank4 = find(sorted_std_devs.begin(),\
                               sorted_std_devs.end(),std_devs[i]) \
                       - sorted_std_devs.begin();
            rank4=-norm_isf((rank4+1.0)/(mat_len+2.0));
            
            cout<<rank1<<' '<<rank2<<' '<<rank3<<' '<<rank4<<' ' \
                <<cont_map[i]<<'\n';
            
        }
    }
    
    df_file.close();
    return 0;
}