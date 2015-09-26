#include <cassert>
#include <fstream>
#include <iostream>
#include "pdbatom.hpp"

int main(int argc, char** argv)
{
    
    using namespace std;
    using namespace cao;
    
    ifstream index;index.open(argv[1]);
    assert(index); 
    string ts,domName;
    PdbAtom atom;
    for (;index>>domName;){
        string domFileName="dompdb/"+domName;
        getline(index,ts);
        
        vector <PdbAtom> atoms;
        vector <XYZ> atom_xyzs;
        ifstream ifile;ifile.open(domFileName.c_str());
        assert(ifile); 
        for (;ifile>>atom;){
            atoms.push_back(atom);
            atom_xyzs.push_back(atom.toXYZ());
        }
        ifile.close();
        /// finished reading atoms 
        
        bool has_close_atoms = false;
        for (size_t i=0;i<atoms.size();i++){
            for (size_t j=i+1;j<atoms.size();j++){
                if (abs(atom_xyzs[i] - atom_xyzs[j]) < 0.5){
                    has_close_atoms = true;
                    cout<<atoms[i]<<'\n'<<atoms[j]<<'\n';
                    cout<<"distance: "<<abs(atom_xyzs[i] - atom_xyzs[j])<<'\n';
                }
            }
        }
        if (has_close_atoms){
            cout<<"file: "<<domName<<'\n'<<'\n';
        }
        
        
    }
    index.close(); 
    return 0;
}
