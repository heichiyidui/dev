#include <fstream>
#include <cassert>
#include <cstdlib>
#include <set>
#include <vector>
#include <algorithm>
#include "util.hpp"


int main(int argc, char** argv)
{
    using namespace std; 
    
    ifstream index;index.open("index/CathDomainList.S100.v3.4.0");
    assert(index); 
    set <string> classes;
    set <string> domNames;
    vector <string> tlist;
    string ts;
    for (;getline(index,ts);){
        strSplit(ts,tlist);
        ts=tlist[1]+'.'+tlist[2]+'.'+tlist[3]+'.'+tlist[4]+'.'+tlist[5];
        if (!classes.count(ts)){
            classes.insert(ts);
            domNames.insert(tlist[0]);
        }
    }
    index.close();
    
    ifstream ifile;ifile.open("index/CathDomainList.S100.v3.4.0");
    assert(ifile); 
    for (;getline(ifile,ts);){
        strSplit(ts,tlist);
        if (domNames.count(tlist[0]))
            cout<<ts<<'\n';
    }
    ifile.close();
    return 0;
}

