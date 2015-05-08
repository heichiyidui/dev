///////////////////////////////////////////////////////////////////////////////
// Copyright (C) 2005-11  Kuang Lin kuanggong@gmail.com                      //
///////////////////////////////////////////////////////////////////////////////

///////////////////////////////////////////////////////////////////////////////
// A small utility program. Given a list of numbers in a file, output some   //
// basic statistics of these numbers.                                        //
///////////////////////////////////////////////////////////////////////////////

#include <iostream>
#include <fstream>
#include <cassert>
#include <cstdlib>
#include <cmath>
#include <algorithm>
#include <vector>
#include <numeric>
#include "util.hpp"

int main(int argc, char** argv){
    using namespace std;
    if (argc<2){
        cerr<<"usage: listdis file [-bn binnum]"<<'\n';
        exit(0);
    }
    
    ifstream ifile;ifile.open(argv[1]);
    if (ifile.fail()){
        cerr<<"usage: listdis file [-bn binnum]"<<'\n';
        exit(0);
    }
    
    string ts;
    double tf;
    vector <double> flist;
    for (;ifile>>ts;)
        if (fromString(ts,tf))
            flist.push_back(tf);
    ifile.close();
    
    size_t size=flist.size();
    if (size<10){
        cerr<<"too few numbers in the file "<<argv[1]<<'\n';
        exit(0);
    }
    
    sort(flist.begin(),flist.end());
    
    if (flist[0]==flist[size-1]){
        cerr<<"All "<<flist[0]<<" in the file "<<argv[1]<<'\n';
        return 0;
    }
    double mean=accumulate(flist.begin(),flist.end(),0.0)/size;
    
    double sdv=0;
    for (size_t i=0;i<size;i++)
        sdv+=(flist[i]-mean)*(flist[i]-mean);
    sdv=sqrt(sdv/(size-1));
    
    cerr<<"min:       "<<flist[0]<<'\n'
        <<"max:       "<<flist[size-1]<<'\n'
        <<"mean:      "<<mean<<'\n'
        <<"sdv:       "<<sdv<<'\n'
        <<"quartiles: "<<flist[0]<<' '<<flist[size/4]<<' '<<flist[size/2]<<' '
                       <<flist[size*3/4]<<' '<<flist[size-1]<<'\n'
        <<"size:      "<<size<<'\n';
    
    size_t binnum;
    if (!getarg2(argc,argv,"-bn",binnum))
        binnum=int(sqrt(double(size))/2.);
    if (binnum<2) return 1;
    
    const double binsize=(flist[size-1]-flist[0])/(binnum-1);
    cerr<<"binnum:    "<<binnum<<'\n'
        <<"binsize:   "<<binsize<<'\n';
    
    vector <size_t> tsum(binnum) ; for (size_t i=0;i<binnum;i++)tsum[i]=0;
    
    for (size_t i=0;i<size;i++)
        tsum[size_t((flist[i]-flist[0]+binsize*0.5)/binsize)]++;
    
    for (size_t i=0;i<binnum;i++)
        cout<<flist[0]+binsize*i<<' '<<double(tsum[i])/size<<'\n';
    
    return 0;
}
