#include <iostream>
#include <fstream>
#include <cassert>
#include <numeric>
#include <algorithm>
#include <vector>
#include <cstdlib>
#include "util.hpp"
#include "seq.hpp"

int main(int argc, char** argv)
{
    using namespace std;
    using namespace cao;
    
    ifstream index;index.open("index/dom.ls");
    assert(index); 
    string ts;
    vector <string> ind;
    for (;index>>ts;)
        ind.push_back(ts);
    index.close();
    
    for (size_t i=0;i<ind.size();i++){
        string domName=ind[i];
        
        ifstream ifile;ifile.open(("maln/" + domName).c_str());
        assert(ifile); 
        ofstream ofile;ofile.open(("caln/"+domName).c_str());
        assert(ofile); 
        
        Seq tseq,domSeq;
        vector <Seq> seqs;
        if (ifile>>tseq){
            seqs.push_back(tseq);
            for (;ifile>>tseq;)
                seqs.push_back(tseq);
            
            vector <Seq> hits;
            for (size_t j=0;j<seqs.size();j++){
                if (seqs[j].name_==domName)
                    domSeq=seqs[j];
                else
                    hits.push_back(seqs[j]);
            }
            ////////////////////////////////////////////////////////////////////
            // output domain and consensus sequences
            Seq outDomSeq,outConSeq; 
            outDomSeq.name_ = domSeq.name_;
            outConSeq.name_ = "consensus";
            
            ////////////////////////////////////////////////////////////////////
            // no blast hits with sequence id > 50%
            if (hits.size()<1){ 
                for (size_t j=0;j<domSeq.length();j++){
                    if (abc().isGap(domSeq[j]))
                        continue;
                    outDomSeq.seq_.push_back(domSeq[j]);
                    outConSeq.seq_.push_back(domSeq[j]);
                }
                ofile<<outDomSeq<<outConSeq;
                continue;
            }
            
            ////////////////////////////////////////////////////////////////////
            // when it has blast hits with sequence id > 50%
            ////////////////////////////////////////////////////////////////////
            
            for (size_t j=0;j<domSeq.length();j++){
                
                if (abc().isResolved(domSeq[j])){
                    outDomSeq.seq_.push_back(domSeq[j]);
                    outConSeq.seq_.push_back(domSeq[j]);
                }else { // domSeq[j] is 'B' or 'X' or 'Z' or '-'
                    
                    size_t AAcounts[abc().size()]; 
                    // abc().size() is the number of resolved states of AA, 20
                    for (size_t k=0;k<abc().size();k++)
                        AAcounts[k]=0;
                    
                    for (size_t k=0;k<hits.size();k++){
                        if (abc().isResolved(hits[k][j])){
                            AAcounts[hits[k][j]]++;
                        }
                    }
                    
                    // 'B' or 'X' or 'Z' must be replaced
                    // '-' only to be replaced when there are > 50% hits
                    if (abc().isGap(domSeq[j])){
                        if (accumulate(AAcounts, AAcounts+abc().size(),size_t()) 
                              *2 < hits.size())
                            continue;
                    }
                    
                    outDomSeq.seq_.push_back(domSeq[j]);
                    if (*max_element(AAcounts, AAcounts+abc().size())==0)
                        outConSeq.seq_.push_back(domSeq[j]);
                    else {
                        size_t AAtype
                            = max_element(AAcounts, AAcounts+abc().size())
                              - AAcounts;
                        outConSeq.seq_.push_back(AAtype);
                    }
                }// finish replacing 'B' or 'X' or 'Z' or '-'
            }
            ofile<<outDomSeq<<outConSeq;
            
        }
        ///////////////////////////////////////////////////////////////////////
        // no alignment all
        ///////////////////////////////////////////////////////////////////////
        else { 
            ifstream sfile;sfile.open(("seq/"+domName).c_str());
            assert(sfile); 
            sfile>>tseq;
            sfile.close();
            ofile<<tseq;
            tseq.name_="consensus";
            ofile<<tseq;
        }
        
        ifile.close();
        ofile.close();
        // break;
    }
    return 0;
}
