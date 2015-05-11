////////////////////////////////////////////////////////////////////////////////
// test the PdbAtom class                                                     //
////////////////////////////////////////////////////////////////////////////////
#include <cassert>
#include <fstream>
#include <iostream>
#include "pdbatom.hpp"

int main(int argc, char** argv)
{
    
    using namespace std;
    cao::PdbAtom atom;
    
    ifstream ifile;ifile.open("t.in");
    assert(ifile); 
    
    for (;ifile>>atom;){
        cout<<atom<<'\n';
    }
    
    for (;!ifile.eof();){
        ifile.clear();
    
        cout<<"restart reading"<<'\n';
        for (;ifile>>atom;){
            cout<<atom<<'\n';
        }
    }
    
    ifile.close();
    return 0;
}
