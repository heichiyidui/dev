////////////////////////////////////////////////////////////////////////////////
// Copyright (C) 2009-10  Kuang Lin kuanggong@gmail.com                       //
////////////////////////////////////////////////////////////////////////////////

////////////////////////////////////////////////////////////////////////////////
// This file defines the class Seq. Protein Sequences (in FASTA format)       //
////////////////////////////////////////////////////////////////////////////////

#ifndef __LK_SEQ_HPP__
#define __LK_SEQ_HPP__

#include <iostream>
#include <fstream>
#include <vector>
#include <algorithm>
#include "util.hpp"
#include "alphabet.hpp"

////////////////////////////////////////////////////////////////////////////////
// begining of the Seq class                                                  //
////////////////////////////////////////////////////////////////////////////////

////////////////////////////////////////////////////////////////////////////////
// Currently we have only FASTA format reading and writing here. It is very   //
// limited. Other important formats, such as the Clustal format, need to be   //
// implemented.                                                               //
////////////////////////////////////////////////////////////////////////////////

namespace cao
{
    
    class Seq{
        public:
            Seq(){}
            Seq(const Seq&seq):name_(seq.name_),seq_(seq.seq_){}
            Seq& operator = (const Seq& seq)
                {seq_=seq.seq_; name_=seq.name_; return *this;}
            ~Seq(){}
            
            string const name() const {return name_;}
            size_t const size() const {return seq_.size();}
            size_t const length()const{return seq_.size();}
            
            size_t const operator [] (size_t ind) const {return seq_[ind];}
            size_t& operator [] (size_t ind) {return seq_[ind];}
            
            friend cao::ostream& operator << (ostream& os, const Seq& seq);
            
            friend cao::istream& operator >> (istream& is, Seq& seq);
            
            std::string name_;
            std::vector <size_t> seq_;
            
            bool operator == (const Seq& other) const {
                return this->seq_ == other.seq_;
            }
            
            bool operator != (const Seq& other) const {
                return !(*this==other);
            }
            
            bool operator <  (const Seq& other) const {
                return this->seq_ < other.seq_;
            }
    };

////////////////////////////////////////////////////////////////////////////////
// end of the Seq class                                                       //
////////////////////////////////////////////////////////////////////////////////

    std::ostream& operator<<(std::ostream&os, const Seq&seq){
        os<<'>'<<seq.name_;
        for (size_t i=0U;i<seq.size();i++){
            if (i%80==0)
                os<<'\n';
            os<<abc().letter(seq[i]);
        }
        os<<'\n';
        return os;
    }// friend output operator
    
    std::istream& operator>>(std::istream&is, Seq&seq){
        string str;
        bool isSeqFound=false;
        for (;getline(is,str);){
            if (str[0]!='>')
                continue;
            
            seq.name_=str.substr(1,str.length()-1); // without the first '>'
        
            string ts;
            for (;is.peek()!='>' && !(is.eof());){
                getline(is,str);
                ts+=str;
            }
            
            toupper(ts);
            
            seq.seq_.clear();
            for (size_t i=0;i<ts.size();i++){
                if (abc().isIncluded(ts[i]))
                    seq.seq_.push_back(abc().state(ts[i]));
            }
            isSeqFound=true;
            break;
        }
        
        if (!isSeqFound)
            is.setstate(ios::badbit);
        return is;
    }// friend input operator 

////////////////////////////////////////////////////////////////////////////////
// end of the IO functions                                                    //
////////////////////////////////////////////////////////////////////////////////
    
}// end of namespace cao

#endif  // __LK_SEQ_HPP__
