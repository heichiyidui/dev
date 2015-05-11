#include <iostream>
#include <fstream>
#include <cassert>
#include <algorithm>
#include <vector>
#include <set>
#include <cstdlib>
#include "util.hpp"
#include "alphabet.hpp"
#include "xyz.hpp"
#include "pdbatom.hpp"
#include "aacid.hpp"
#include "seq.hpp"
#include "aacid.hpp"
#include "comptetr.hpp"

// g++ -c predicates.cxx -O2
// g++ -c tetgen.cxx -O2
// g++ getTet.cpp predicates.o tetgen.o -O2

int main(int argc, char** argv)
{
    using namespace std;
    using namespace cao;
    
    ifstream index;index.open("t.ls");
    assert(index); 
    vector <string> ind;
    string ts;
    for (;index>>ts;){
        ind.push_back(ts);
        getline(index,ts);
    }
    index.close();
    ////////////////////////////////////////////////////////////////////////////
    // finished reading the index of domain names                             //
    ////////////////////////////////////////////////////////////////////////////
    
    for (size_t l=0;l<ind.size();l++){
        string domName=ind[l];
        
        ifstream ifile;ifile.open(("ccbc/"+domName).c_str());
        assert(ifile); 
        PdbAtom atom;
        vector <XYZ> coords; // coordinates of side-chain center and backbone C
        for (;ifile>>atom;){
            coords.push_back(atom.toXYZ());
        }
        // for each residue, the CB, then the C atom coordinates
        ifile.close();
        ////////////////////////////////////////////////////////////////////////
        // finished reading pdb residues                                      //
        ////////////////////////////////////////////////////////////////////////
        
        vector <pair<size_t, size_t> > conMap=computeTetra(coords);
        
        ofstream ofile;ofile.open(("conDef/"+domName).c_str());
        assert(ofile); 
        for (size_t i=0;i<conMap.size();i++){
            if (conMap[i].first%2 || conMap[i].second%2)
                continue; // we don't care about main-chain contacts
            if (conMap[i].first > conMap[i].second)
                continue;
            if (abs(coords[conMap[i].first]-coords[conMap[i].second]) > 8)
                continue;
            
            ofile<<conMap[i].first/2<<' '<<conMap[i].second/2<<'\n';
        }
        ofile.close();
        
        cout<<domName<<'\n';
//        break;
    }
    return 0;
}
