#include <iostream>
#include <fstream>
#include <cassert>
#include <algorithm>
#include <vector>
#include <set>
#include <cstdlib>
#include "util.hpp"
#include "seq.hpp"

int main(int argc, char** argv)
{
    using namespace std;
    using namespace cao;
    
    ifstream index;index.open("index/CathDomainList.H");
    assert(index);
    string ts;
    vector <string> ind;
    for (;index>>ts;){
        ind.push_back(ts);
        getline(index,ts);
    } 
    index.close();
    
    for (size_t l=0;l<ind.size();l++){
        string domName=ind[l];
        
        ifstream ifile;ifile.open(("caln/"+domName).c_str());
        assert(ifile); 
        Seq dseq,sseq;
        ifile>>dseq>>sseq;
        ifile.close();
        
        bool isAAInDomain[dseq.length()];
        for (size_t i=0;i<dseq.length();i++)
            isAAInDomain[i]=!abc().isGap(dseq[i]);
        
        string states;
        for (size_t i=0;i<dseq.length();i++)
            states=states+"U";
        
        ifstream dfile;dfile.open(("ssdef/"+domName).c_str());
        assert(dfile); 
        for (;dfile>>ts;){
            if (ts=="#") break;
            getline(dfile,ts);
        }
        getline(dfile,ts);
        for (;getline(dfile,ts);){
            if (ts[13]=='!')
                continue;
            size_t loc;
            fromString(ts.substr(5,5),loc);
            loc--;
            char AA=ts[13];
            char SS=ts[16];
            if (dseq[loc] != abc().state(AA))
                cout<<domName<<'\n';
            states[loc]=SS;
        }
        dfile.close();
        for (size_t i=0;i<dseq.length();i++){
            if (isAAInDomain[i] && states[i]=='U')
                cout<<"missing SS in "<<domName<<' '<<i<<'\n';
            if (!isAAInDomain[i] && states[i] != 'U')
                cout<<"over assign SS in "<<domName<<'\n';
        }
        
        ////////////////////////////////////////////////////////////////////////
        // from 8 states of DSSP definition to 3 states SS definition. 
        ////////////////////////////////////////////////////////////////////////
        
        string states3;
        for (size_t i=0;i<states.length();i++){
            if (states[i]=='U'){
                states3+="U";
            }else if (states[i]=='H' || states[i]=='G'){
                states3+="H";
            }else if (states[i]=='E' || states[i]=='B'){
                states3+="E";
            }else if (states[i]=='S' || states[i]=='T' 
                   || states[i]=='I' || states[i]==' '){
                states3+="C";
            }else{
                cerr<<"What? unrecognized states from DSSP: "
                    <<states[i]<<" in domain "<<domName<<'\n';
                states3+="C";
                break;
            }
        }
        ////////////////////////////////////////////////////////////////////////
        // from 3 states to 7 states
        ////////////////////////////////////////////////////////////////////////
        string states7=states3;
        for (size_t i=0;i<states3.length()-2;i++){
            ts=states3.substr(i,3);
            if (ts=="CHH" || ts=="EHH"){
                states7[i+1]='G';
            }else if (ts=="HHC" || ts=="HHE"){
                states7[i+1]='I';
            }else if (ts=="HEE" || ts=="CEE"){
                states7[i+1]='D';
            }else if (ts=="EEC" || ts=="EEH"){
                states7[i+1]='F';
            }
        }
        for (size_t i=0;i<states7.length();i++){
            cout<<states7[i]<<'\n';
        }
    }
    return 0;
}
