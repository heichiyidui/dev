////////////////////////////////////////////////////////////////////////////////
// Copyright (C) 2009-10  Kuang Lin kuanggong@gmail.com                       //
////////////////////////////////////////////////////////////////////////////////

////////////////////////////////////////////////////////////////////////////////
//                                                                            //
// Alphabet class, the dictionary of int and char for sequence presentations. //
//                                                                            //
// Protein residues (amino acids) or nucleotides are presented as characters  //
// in files of fasta and other human-readable formats. In computing, they     //
// are integers.                                                              //
//                                                                            //
// We (might) need different dictionaries for the translation tasks.          //
//                                                                            //
// Here I copied a lot from Bio++ (kimura.univ-montp2.fr/BioPP/index.html)    //
////////////////////////////////////////////////////////////////////////////////

#ifndef __ALPHABET_HPP__
#define __ALPHABET_HPP__

#include <cstdlib>
#include <iostream>
#include <cmath>
#include <vector>
#include <map>
#include <algorithm>

////////////////////////////////////////////////////////////////////////////////
// The Alphabet class                                                         //
////////////////////////////////////////////////////////////////////////////////

namespace cao{
    
    using namespace std;
    
    class Alphabet{
        
        public:
            
            Alphabet(string alphabetType="IUPAC Protein Ambiguous"){
            
                if (alphabetType!="IUPAC Protein Ambiguous"){
                    cerr<<"In initializing alphabet: unrecognized type: "
                        <<alphabetType<<'\n'
                        <<"Try \"IUPAC Protein Ambiguous\" instead...\n" ;
                    exit(0);
                }// well, only one type of alphabet recognized now...
                
                alphabetType_=alphabetType;
                
                nStates_=24;
                nResolvedStates_=20;
                
                letters_  = "ARNDCQEGHILKMFPSTWYVBZX-"; 
                    // note they are all capital letters here.
            
                gaps_     = "-";
                unknowns_ = "X0O?"; 
                    // gaps_[0] and unknowns_[0] should be always in letters_
                    // the first nResolvedStates_ characters in letters_
                    // should be the resolvable ones
            
                string tlist[]={"ALA","ARG","ASN","ASP","CYS", // ARNDC    
                                "GLN","GLU","GLY","HIS","ILE", // QEGHI    
                                "LEU","LYS","MET","PHE","PRO", // LKMFP    
                                "SER","THR","TRP","TYR","VAL", // STWYV 
                                "ASX","GLX","UNK","GAP"};      // BZX-
            
                for (size_t i=0;i<24;i++)
                    abbrs_.push_back(tlist[i]);
                
                // string stdSolventAbbr[]=
                // {"ADE", "CYT", "GUA", "INO", "THY", "URA",
                //  "WAT", "HOH", "TIP", "H2O", "DOD", "MOH"};
            
                abbrStateMap_["ALA"] = 0  ; // A
                abbrStateMap_["ARG"] = 1  ; // R
                abbrStateMap_["ASN"] = 2  ; // N
                abbrStateMap_["ASP"] = 3  ; // D
                abbrStateMap_["CYS"] = 4  ; // C
                abbrStateMap_["GLN"] = 5  ; // Q
                abbrStateMap_["GLU"] = 6  ; // E
                abbrStateMap_["GLY"] = 7  ; // G
                abbrStateMap_["HIS"] = 8  ; // H
                abbrStateMap_["ILE"] = 9  ; // I
                abbrStateMap_["LEU"] = 10 ; // L
                abbrStateMap_["LYS"] = 11 ; // K
                abbrStateMap_["MET"] = 12 ; // M
                abbrStateMap_["PHE"] = 13 ; // F
                abbrStateMap_["PRO"] = 14 ; // P
                abbrStateMap_["SER"] = 15 ; // S
                abbrStateMap_["THR"] = 16 ; // T
                abbrStateMap_["TRP"] = 17 ; // W
                abbrStateMap_["TYR"] = 18 ; // Y
                abbrStateMap_["VAL"] = 19 ; // V
            
                abbrStateMap_["ASX"] = 20 ; // B
                abbrStateMap_["GLX"] = 21 ; // Z
                abbrStateMap_["UNK"] = 22 ; // X
                abbrStateMap_["XXX"] = 22 ; // X
                abbrStateMap_["GAP"] = 23 ; // -
                
                abbrStateMap_["ALB"] = 0  ; // A
                abbrStateMap_["CSH"] = 4  ; // C
                abbrStateMap_["CYH"] = 4  ; // C
                abbrStateMap_["CSS"] = 4  ; // C
                abbrStateMap_["CYX"] = 4  ; // C
                abbrStateMap_["PCA"] = 6  ; // E
                abbrStateMap_["ILU"] = 9  ; // I
                abbrStateMap_["PR0"] = 14 ; // P
                abbrStateMap_["PRZ"] = 14 ; // P
                abbrStateMap_["HYP"] = 14 ; // P
                abbrStateMap_["TRY"] = 17 ; // W
                
            }// end of the constructor
    
            
            ////////////////////////////////////////////////////////////////////
            // about the alphabet                                             //
            ////////////////////////////////////////////////////////////////////
            
            string alphabetType() const {return alphabetType_;}
                //  a short description of the alphabet
            
            bool hasGap()     const{return gaps_.size()>0;}   
            bool hasUnknown() const{return unknowns_.size()>0;}
            
            size_t size()     const{return nResolvedStates_;}
                // say, 20 amino acids 
            size_t nStates()  const{return nStates_;}
                // say, 24 for amino acids + 'X' + 'B' + 'Z' + '-' 
            
            size_t gapState() const{return state(gaps_[0]);}
            
            ////////////////////////////////////////////////////////////////////
            // using the alphabet                                             //
            ////////////////////////////////////////////////////////////////////
            
            ////////////////////////////////////////////////////////////////////
            // is included
            // bool isIncluded (const size_t state) const;
            // bool isIncluded (const char letter)  const;
            // bool isIncluded (const string abbr)  const;
            
            bool isIncluded (const size_t state) const
                {return state<nStates_;}
    
            bool isIncluded (const char letter) const
                {char c=toupper(letter);
                return   letters_.find(c) != string::npos 
                      ||    gaps_.find(c) != string::npos  
                      ||unknowns_.find(c) != string::npos;}
    
            bool isIncluded (const string abbr) const
                {return abbrStateMap_.find(abbr) != abbrStateMap_.end();}
            
            ////////////////////////////////////////////////////////////////////
            // is resolved
            // bool isResolved (const size_t state) const;
            // bool isResolved (const char letter)  const;
            // bool isResolved (const string abbr)  const;
            
            bool isResolved (const size_t state) const
                {return state < nResolvedStates_;}
    
            bool isResolved (const char letter) const
                {return letters_.substr(0,nResolvedStates_).
                    find(toupper(letter)) != string::npos;}
    
            bool isResolved (const string abbr) const
                {return isIncluded(abbr) 
                    &&  isResolved(state(abbr));}
            
            ////////////////////////////////////////////////////////////////////
            // is gap
            // bool isGap (const size_t state)  const
            // bool isGap (const char  letter)  const
            
            bool isGap (const size_t state)  const
                {if (!isIncluded(state)) return false; 
                 return isGap(letter(state));}
                 
            bool isGap (const char letter) const
                {return gaps_.find(toupper(letter)) != string::npos;}
            
            ////////////////////////////////////////////////////////////////////
            // is unknown
            // bool isUnknown (const size_t state)  const
            // bool isUnknown (const char  letter)  const
            
            bool isUnknown (const size_t state)  const
                {if (!isIncluded(state)) return false; 
                 return isUnknown(letter(state));}
            
            bool isUnknown (const char   letter) const
                {return unknowns_.find(toupper(letter)) != string::npos;}
            
            ////////////////////////////////////////////////////////////////////
            // get letter, state and abbreviation
            // 
            // char   letter(const size_t state) const;
            // string abbr  (const size_t state) const;
            // size_t state (const char letter)  const;
            // size_t state (const string abbr)  const;
            
            char   letter(const size_t state) const{return letters_[state];}
            string abbr  (const size_t state) const{return abbrs_  [state];}
            
            size_t state(const char letter) const{
                char c=toupper(letter);
                if (unknowns_.find(c)!= string::npos) c=unknowns_[0];
                if (gaps_    .find(c)!= string::npos) c=gaps_    [0];
                return letters_.find_first_of(c);}
            
            size_t state(const string abbr) const
                {return abbrStateMap_.find(abbr)->second;}
            
        ////////////////////////////////////////////////////////////////////////
        // the private members                                                //
        ////////////////////////////////////////////////////////////////////////
            
        private:
            
            string alphabetType_;
            
            size_t nStates_;
            size_t nResolvedStates_;
            
            string letters_;
            string gaps_;
            string unknowns_;
            
            vector <string> abbrs_;  
            map <string, size_t> abbrStateMap_;
                // abbreviations of amino acid names
                // or residue name in PDB ATOM (coordinate) entries.
                // for CATH it is fine to use the standard set.
                // for PDB and SCOP, residues have many different variations
                // and different abbreviations for the same amino acid (in 
                // differently modified forms)
            
    };// end of class Alphabet
    
    inline Alphabet& abc(){
        static Alphabet * abc=new Alphabet("IUPAC Protein Ambiguous");
        return *abc;
    } // this is a reference to a global constant
      // the function leaks, because 'abc' was never destructed.
      // however, this amount of leaking is not serious.
        
} // end of namespace cao

////////////////////////////////////////////////////////////////////////////////
// end of the Alphabet class                                                  //
////////////////////////////////////////////////////////////////////////////////


#endif // __ALPHABET_HPP__
