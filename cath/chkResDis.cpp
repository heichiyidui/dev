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
        
        vector <vector <PdbAtom> > residues;
        for (size_t i=0;i<atoms.size();i++){
            bool isResFound=false;
            for (size_t j=0;j<residues.size();j++){
                if (isFromSameResidue(atoms[i],residues[j][0])){
                    isResFound=true;
                    residues[j].push_back(atoms[i]);
                    continue;
                }
            }
            if (!isResFound){
                vector <PdbAtom> tRes;
                tRes.push_back(atoms[i]);
                residues.push_back(tRes);
            }
        }
        
        bool hasStrangeResidue=false;
        for (size_t i=0;i<residues.size();i++){
            for (size_t j=0;j<residues[i].size();j++){
                for (size_t k=j+1;k<residues[i].size();k++){
                    XYZ v=residues[i][j].toXYZ() - residues[i][k].toXYZ();
                    if (abs(v)>10.0 || abs(v)<0.8){
                        cout<<residues[i][j]
                            <<residues[i][k]
                            <<"dis: "<<abs(v)<<'\n';
                        hasStrangeResidue=true;
                    }
                }
            }
        }
        
        if (!hasStrangeResidue)
            continue;
        cout<<"e "<<domFileName<<" & "<<'\n'
            <<"r "<<domFileName<<" & \n"
            <<"www.cathdb.info/domain/"<<domName.substr(0,7)<<'\n'
            <<"www.rcsb.org/pdb/files/"<<domName.substr(0,4)<<".pdb\n"
            <<'\n';
        
    }
    index.close(); 
    return 0;
}
