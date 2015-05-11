#include <fstream>
#include <cassert>
#include <cstdlib>
#include <algorithm>
#include "util.hpp"
#include "aacid.hpp"
#include "alphabet.hpp"
#include "pdbatom.hpp"

int main(int argc, char** argv)
{
    using namespace std; 
    using namespace cao;
    using namespace Aacid;
    
    ifstream index;index.open("t.ls");
    assert(index); 
    string ts;
    vector <string> ind;
    for (;index>>ts;){
        ind.push_back(ts);
        getline(index,ts);
    }
    index.close();
    ////////////////////////////////////////////////////////////////////////////

    size_t completeAaSum=0;
    size_t incompleteAaSum=0;
    
    for (size_t l=0;l<ind.size();l++){
        string domName=ind[l];
        
        ifstream ifile;ifile.open(("dompdb/"+domName).c_str());
        assert(ifile); 
        vector <PdbAtom> atoms;
        PdbAtom atom;
        for (;ifile>>atom;){
            size_t aaType = abc().state(atom.resName_);
            
            if (!Aacid::isAtomInAA(atom.atomName_, aaType)){
                cout<<atom;
                cout<<domName<<'\n';
                continue;
            }
                
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
        ////////////////////////////////////////////////////////////////////////
        // finished reading pdb residues                                      //
        ////////////////////////////////////////////////////////////////////////
        
        ofstream ofile;ofile.open(("ccbc/"+domName).c_str());
        assert(ofile); 
        for (size_t i=0;i<residues.size();i++){
            size_t type = abc().state(residues[i][0].resName_);
            
            XYZ c,ca,n,cb;
            PdbAtom caAtom,cAtom;
            for (size_t j=0;j<residues[i].size();j++){
                const string atomName=residues[i][j].atomName_;
                if (atomName==" CA ") {
                    ca=residues[i][j].toXYZ();
                    caAtom=residues[i][j];
                }
                if (atomName==" C  ") {
                    c =residues[i][j].toXYZ();
                    cAtom=residues[i][j];
                }
                if (atomName==" N  ") n =residues[i][j].toXYZ();
            } // this three atoms are guarrented to be there. 
            
            ////////////////////////////////////////////////////////////////////
            // when the side-chain is NOT complete                            //
            // we need to guess the coordination of the center                //
            ////////////////////////////////////////////////////////////////////
            if (residues[i].size()<Aacid::AA_ATOM_NUM[type]){ 
                XYZ v1=n-ca;
                XYZ v2=c-ca;
                XYZ v3=v1*v2;
                v2=v1*v3;
                v1/=abs(v1); v2/=abs(v2); v3/=abs(v3);
                cb=ca+v1*Aacid::AA_CB[type][0]
                     +v2*Aacid::AA_CB[type][1]
                     +v3*Aacid::AA_CB[type][2];
                
                caAtom=cb;
                caAtom.atomName_=" CB ";
                ofile<<caAtom<<cAtom;
                incompleteAaSum++;
                continue;
            }
            
            ////////////////////////////////////////////////////////////////////
            // when the side-chain is complete                                //
            // we need to calculate the coordination of the center            //
            ////////////////////////////////////////////////////////////////////
            
            for (size_t j=0;j<residues[i].size();j++){
                const string atomName=residues[i][j].atomName_;
                if (atomName==" C  " || atomName==" N  " || atomName==" O  ")
                    continue;
                
                // C 12, N 14, S 32, O 16;
                if (atomName[1]=='C')
                    cb+=(residues[i][j].toXYZ()-ca) * 12;
                if (atomName[1]=='N')
                    cb+=(residues[i][j].toXYZ()-ca) * 14;
                if (atomName[1]=='S')
                    cb+=(residues[i][j].toXYZ()-ca) * 32;
                if (atomName[1]=='O')
                    cb+=(residues[i][j].toXYZ()-ca) * 16;
            }
            
            cb/=Aacid::SIDECHAIN_MASS[type];
            cb+=ca;
            caAtom=cb;
            caAtom.atomName_=" CB ";
            ofile<<caAtom<<cAtom;
            completeAaSum++;
            ////////////////////////////////////////////////////////////////////
            // OK, now we have all the atoms we need                          //
            ////////////////////////////////////////////////////////////////////
        }
        ofile.close();
//        break;
    }
    cout<<incompleteAaSum<<' '<<completeAaSum<<'\n';
    return 0;
}

