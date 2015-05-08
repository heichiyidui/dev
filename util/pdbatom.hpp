////////////////////////////////////////////////////////////////////////////////
// Copyright (C) 2009-10 Kuang Lin kuanggong@gmail.com                        //
////////////////////////////////////////////////////////////////////////////////

////////////////////////////////////////////////////////////////////////////////
// This file defines ATOM entry from PDB files,(acs, bbs, cas etc.)           //
////////////////////////////////////////////////////////////////////////////////

#ifndef __LK_PDBATOM_HPP__
#define __LK_PDBATOM_HPP__

#include <sstream>
#include <iomanip>
#include "xyz.hpp"
#include "util.hpp"

////////////////////////////////////////////////////////////////////////////////
// the PdbAtom class                                                          //
////////////////////////////////////////////////////////////////////////////////
namespace cao
{
    using namespace std;
    class PdbAtom{
        
        public:
        
            PdbAtom():recordName_("ATOM  "), atomNum_(0), atomName_(" CA "),
                altLoc_(' '), resName_("ALA"), chainId_(' '), resNum_(0),
                insertCode_(' '), x(0.0F),y(0.0F),z(0.0F), 
                occu_(1.0F),tempFactor_(0.0F),element_("  "),charge_("  "){};
            
            bool operator == (const PdbAtom& other) const {
                if (   recordName_ != other.recordName_ 
                    || atomNum_    != other.atomNum_ 
                    || atomName_   != other.atomName_ 
                    || altLoc_     != other.altLoc_ 
                    || resName_    != other.resName_ 
                    || chainId_    != other.chainId_ 
                    || resNum_     != other.resNum_ 
                    || insertCode_ != other.insertCode_ 
                    || x           != other.x
                    || y           != other.y
                    || z           != other.z )
                    return false;
                return true;
            }
            
            bool operator != (const PdbAtom& other) const {
                return !(*this==other);
            }
            
            ////////////////////////////////////////////////////////////////////
            // input and output 
            friend cao::ostream& operator << (ostream& os, const PdbAtom& atom);
            friend cao::istream& operator >> (istream& is, PdbAtom& atom);
            
            PdbAtom& operator = (const XYZ& xyz)
                { x=xyz.x; y=xyz.y; z=xyz.z; return *this;}
            
            XYZ toXYZ() const {return XYZ(x,y,z);}
            
            ////////////////////////////////////////////////////////////////////
            // items of a PDB ATOM entry                                      //
            ////////////////////////////////////////////////////////////////////
            
            string  recordName_;// "ATOM  " or "HETATM" 
            
            int     atomNum_;   // atom serial number, can be negative!
            string  atomName_;  // atom name, like " CA ", " N  "
            
            char    altLoc_;    // alternative location indicator
            string  resName_;   // residue abbreviations, like "ALA", "VAL" etc.
            
            char    chainId_;   // chain identifier
            int     resNum_;    // residue sequence number, can be negative!
            
            char    insertCode_;// code for insertion of residues
            float   x,y,z;      // coordinates 
            
            float   occu_;      // occupancy
            float   tempFactor_;// temperature factor
            
            string  element_;   // element, right justified, " N", " C" etc.
            string  charge_;    // charge on the atom, "2+", "1-" etc.
            
    };
////////////////////////////////////////////////////////////////////////////////
// the end of PdbAtom class                                                   //
////////////////////////////////////////////////////////////////////////////////
    
////////////////////////////////////////////////////////////////////////////////
// friend output and input operators                                          //
////////////////////////////////////////////////////////////////////////////////
    
    ostream& operator << (ostream& os, const PdbAtom& atom){
        os<<atom.recordName_;
    
        os.setf(ios::right,ios::adjustfield);
        os<<setw(5)<<atom.atomNum_<<' '<<setw(4)<<atom.atomName_
          <<atom.altLoc_ <<atom.resName_<<' '<<atom.chainId_ 
          <<setw(4)<<atom.resNum_  <<atom.insertCode_<<"   ";
    
        os.setf(ios::fixed,ios::floatfield);
        os.precision(3);
        os<<setw(8)<<atom.x<<setw(8)<<atom.y<<setw(8)<<atom.z;
    
        os.precision(2);
        os<<setw(6)<<atom.occu_<<setw(6)<<atom.tempFactor_  // the last four are
          <<"          "<<atom.element_<<atom.charge_<<'\n';// often not used...
        
        return os;
    }// friend output operator
     // output are often less error-prone 
    
    istream& operator >> (istream& is, PdbAtom& atom){
        string str;
        for (;getline(is,str);){
            if (str.length()<60)
                continue;
            if (str.substr(0,6) != "ATOM  " && str.substr(0,6) != "HETATM") 
                continue; // we are only bothered about ATOM and HETATM entries.
            
            atom.recordName_=str.substr(0,6);
            atom.atomName_  =str.substr(12,4);
            atom.altLoc_    =str[16];
            atom.resName_   =str.substr(17,3);
            atom.chainId_   =str[21];
            atom.insertCode_=str[26];
        
            bool isFine = true; // if the reading of xyx is going ok
        
            isFine = isFine && fromString(str.substr(6 ,5),atom.atomNum_);
            isFine = isFine && fromString(str.substr(22,4),atom.resNum_);
        
            isFine = isFine && fromString(str.substr(30,8),atom.x );
            isFine = isFine && fromString(str.substr(38,8),atom.y);
            isFine = isFine && fromString(str.substr(46,8),atom.z);
            
            if (!isFine){
                is.setstate(ios::badbit);
                return is;
            }
        
            if (str.length()>64){
                ::fromString(str.substr(54,6),atom.occu_);
                ::fromString(str.substr(60,6),atom.tempFactor_);
            }
            if (str.length()>=80){
                atom.element_=str.substr(76,2);
                atom.charge_ =str.substr(78,2);
            }// the standard line should be 80 charactors, 
             // but some files don't have the last stuffs
            
            break;
        }
        return is;
    }
    // input might be more trouble here. 
    
////////////////////////////////////////////////////////////////////////////////
// some minor functions                                                       //
////////////////////////////////////////////////////////////////////////////////
    
    inline bool isFromSameResidue(const PdbAtom&a1,const PdbAtom&a2)
        {return   a1.resName_    == a2.resName_
               && a1.chainId_    == a2.chainId_
               && a1.resNum_     == a2.resNum_ 
               && a1.insertCode_ == a2.insertCode_;}
    
    inline bool operator < (const PdbAtom & a1, const PdbAtom & a2)
        {if (a1.resNum_ != a2.resNum_) 
            return a1.resNum_ < a2.resNum_;
        return a1.atomNum_ < a2.atomNum_;}
    
    inline const float CNOSMass(string const & atomName){
        switch (atomName[1]){
            case 'C': {return 12.0107;}
            case 'N': {return 14.0067;}
            case 'O': {return 15.9994;}
            case 'S': {return 32.0655;}
            default : {return 0.0F;}
        }
    }// the relative MASS of C, N, O, or S atoms
    
} // end of namespace cao
#endif  // __LK_PDBATOM_HPP__
