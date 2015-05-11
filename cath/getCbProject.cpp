#include <fstream>
#include <cassert>
#include <cstdlib>
#include <algorithm>
#include "aacid.hpp"
#include "alphabet.hpp"
#include "pdbatom.hpp"
#include "xyz.hpp"

int main(int argc, char** argv)
{
    using namespace std; 
    using namespace cao;
    using namespace Aacid;
    
    ifstream index;index.open("index/CathDomainList.S35");
    assert(index); 
    string ts;
    vector <string> ind;
    for (;index>>ts;){
        ind.push_back(ts);
        getline(index,ts);
    }
    index.close();
    ////////////////////////////////////////////////////////////////////////////

    size_t AAsum[20]; 
    XYZ AAprojects[20];
    for (size_t i=0;i<20;i++){
        AAsum[i]=0;
        AAprojects[i].x=0;AAprojects[i].y=0;AAprojects[i].z=0;
    }
    
    for (size_t l=0;l<ind.size();l++){
        string domName=ind[l];
        
        ifstream ifile;ifile.open(("cdom/"+domName).c_str());
        assert(ifile); 
        vector <PdbAtom> atoms;
        PdbAtom atom;
        for (;ifile>>atom;)
            atoms.push_back(atom);
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
        
        for (size_t i=0;i<residues.size();i++){
            size_t type = abc().state(residues[i][0].resName_);
            
            if (residues[i].size()<Aacid::AA_ATOM_NUM[type]) // NOT complete
                continue;                              
            
            XYZ c,ca,n,cb;
            for (size_t j=0;j<residues[i].size();j++){
                const string atomName=residues[i][j].atomName_;
                if (atomName==" CA ") ca=residues[i][j].toXYZ();
                if (atomName==" C  ") c =residues[i][j].toXYZ();
                if (atomName==" N  ") n =residues[i][j].toXYZ();
            }
            
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
            
            ////////////////////////////////////////////////////////////////////
            // OK, now we have all the atoms we need                          //
            ////////////////////////////////////////////////////////////////////
            XYZ v1=n-ca;
            XYZ v2=c-ca;
            XYZ v3=v1*v2;
            v2=v1*v3;
            
            v1/=abs(v1); v2/=abs(v2); v3/=abs(v3);
            
            AAprojects[type].x+=cb.dot(v1);
            AAprojects[type].y+=cb.dot(v2);
            AAprojects[type].z+=cb.dot(v3);
            AAsum[type]++;
        }
//        break;
    }
    
    for (size_t i=0;i<20;i++){
        AAprojects[i]/=AAsum[i];
        cout<<AAprojects[i].x<<','
            <<AAprojects[i].y<<','
            <<AAprojects[i].z<<'\n';
    }
    return 0;
}

