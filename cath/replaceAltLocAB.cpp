#include <fstream>
#include <cassert>
#include <cstdlib>
#include <set>
#include <algorithm>
#include "pdbatom.hpp"

int main(int argc, char** argv)
{
    using namespace std; 
    using namespace cao;
    
    ifstream index;index.open("ab_abc.ls");
    assert(index); 
    string ts,domName;
    PdbAtom atom;
    for (;index>>domName;){
        string domFileName="pdb/"+domName;
        getline(index,ts);
        
        vector <PdbAtom> atoms;
        ifstream ifile;ifile.open(domFileName.c_str());
        assert(ifile); 
        for (;ifile>>atom;){
            atoms.push_back(atom);
        }
        ifile.close();
        
        //// finished reading atoms ////////////////////////////////////////////
        
        ifstream pfile;pfile.open(("rcsb/"+domName.substr(0,4)+".pdb").c_str());
        
        assert(pfile); 
        vector <PdbAtom> pAtoms;
        for (;pfile>>atom;)
            pAtoms.push_back(atom);
        pfile.close();
        
        //use 'A' atoms instead of the 'B' and 'C' atoms;
        for (size_t i=0;i<atoms.size();i++){
            if (atoms[i].altLoc_==' ' || atoms[i].altLoc_=='A')
                continue;
            
            bool isCathReplacementFound=false;
            for (size_t j=0;j<atoms.size();j++){
                if (atoms[i].altLoc_!=' ' && atoms[i].altLoc_!='A')
                    continue;
                if (  isFromSameResidue(atoms[j],atoms[i]) 
                    &&atoms[j].atomName_==atoms[i].atomName_){
                    isCathReplacementFound=true;
                    atoms[i].occu_=-199.0;
                }
            }
            if (isCathReplacementFound)
                continue;
            
            bool isPdbReplacementFound=false;
            for (size_t j=0;j<pAtoms.size();j++){
                if (!isFromSameResidue(atoms[i],pAtoms[j]))
                    continue;
                if (pAtoms[j].altLoc_!=' ' && pAtoms[j].altLoc_!='A')
                    continue;
                if (pAtoms[j].atomName_!= atoms[i].atomName_)
                    continue;
                isPdbReplacementFound=true;
                atoms[i]=pAtoms[j];
                break;
            }
            if (!isPdbReplacementFound){
                cout<<atoms[i]<<domName<<'\n';
            }
        }
        
        ofstream ofile;ofile.open(("pdb2/"+domName).c_str());
        assert(ofile); 
        for (size_t i=0;i<atoms.size();i++)
            if (atoms[i].occu_>=0)
                ofile<<atoms[i];
        ofile.close();
    }
    index.close();
    return 0;
}


