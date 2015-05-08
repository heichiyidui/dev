////////////////////////////////////////////////////////////////////////////////
// Copyright (C) 2001-11  Kuang Lin kuanggong@gmail.com                       //
////////////////////////////////////////////////////////////////////////////////

////////////////////////////////////////////////////////////////////////////////
// this file define amino acid names and some statics                         //
////////////////////////////////////////////////////////////////////////////////

#ifndef LK_AACID_HPP
#define LK_AACID_HPP

#include <string>
#include <cctype>
#include <cstring>
#include <set>
#include <algorithm>

namespace Aacid{

using namespace std;

////////////////////////////////////////////////////////////////////////////////
// some statistics of amino acids                                             //
// CATH 3.4.0 H representatives 2477 domains, 326881 Amino Acids              //
////////////////////////////////////////////////////////////////////////////////

float const RES_OCU[20]={  //the Residue Occurrences
    0.08128,0.05291,0.04349,0.05817,0.01346, // ARNDC
    0.03891,0.07167,0.07071,0.02208,0.05711, // QEGHI
    0.09290,0.06155,0.01809,0.03976,0.04490, // LKMFP
    0.05838,0.05552,0.01384,0.03511,0.07016};// STWYV
    
const float LOG_RES_OCU[20]={ 
    -2.5099,-2.9392,-3.1353,-2.8443,-4.3084, // ARNDC
    -3.2464,-2.6356,-2.6492,-3.8129,-2.8628, // QEGHI
    -2.3763,-2.7879,-4.0123,-3.2248,-3.1032, // LKMFP
    -2.8407,-2.8910,-4.2811,-3.3493,-2.6569};// STWYV 
    
const float SIDECHAIN_MASS[20]={
     24 ,102,66 ,68 ,56 ,  // ARNDC
     78 ,80 ,12 ,88 ,60 ,  // QEGHI
     60 ,74 ,80 ,96 ,48 ,  // LKMFP
     40 ,52 ,134,112,48 }; // STWYV
     
const size_t AA_ATOM_NUM[20]={ // how many heavy atoms in a residue
     5,11, 8, 8, 6, // ARNDC
     9, 9, 4,10, 8, // QEGHI
     8, 9, 8,11, 7, // LKMFP
     6, 7,14,12, 7};// STWYV
     
const std::string AA_ATOM_NAME[20][14]={// name of heavy atoms in a residue
/* A */{" C  "," CA "," CB "," N  "," O  ",},
/* R */{" C  "," CA "," CB "," CG "," CD "," CZ "," N  "," NE "," NH1"," NH2"," O  ",},
/* N */{" C  "," CA "," CB "," CG "," N  "," ND2"," O  "," OD1",},
/* D */{" C  "," CA "," CB "," CG "," N  "," O  "," OD1"," OD2",},
/* C */{" C  "," CA "," CB "," N  "," O  "," SG ",},
/* Q */{" C  "," CA "," CB "," CG "," CD "," N  "," NE2"," O  "," OE1",},
/* E */{" C  "," CA "," CB "," CG "," CD "," N  "," O  "," OE1"," OE2",},
/* G */{" C  "," CA "," N  "," O  ",},
/* H */{" C  "," CA "," CB "," CG "," CD2"," CE1"," N  "," ND1"," NE2"," O  ",},
/* I */{" C  "," CA "," CB "," CG1"," CG2"," CD1"," N  "," O  ",},
/* L */{" C  "," CA "," CB "," CG "," CD1"," CD2"," N  "," O  ",},
/* K */{" C  "," CA "," CB "," CG "," CD "," CE "," N  "," NZ "," O  ",},
/* M */{" C  "," CA "," CB "," CG "," CE "," N  "," O  "," SD ",},
/* F */{" C  "," CA "," CB "," CG "," CD1"," CD2"," CE1"," CE2"," CZ "," N  "," O  ",},
/* P */{" C  "," CA "," CB "," CG "," CD "," N  "," O  ",},
/* S */{" C  "," CA "," CB "," N  "," O  "," OG ",},
/* T */{" C  "," CA "," CB "," CG2"," N  "," O  "," OG1",},
/* W */{" C  "," CA "," CB "," CG "," CD1"," CD2"," CE2"," CE3"," CZ2"," CZ3"," CH2"," N  "," NE1"," O  ",},
/* Y */{" C  "," CA "," CB "," CG "," CD1"," CD2"," CE1"," CE2"," CZ "," N  "," O  "," OH ",},
/* V */{" C  "," CA "," CB "," CG1"," CG2"," N  "," O  ",}};

bool isAtomInAA(const string atomName, const size_t state){
    
    static set <string> names[20]; // the names of heavy atoms in amino acids
    static bool isInited=false;    // if names[] is initialized or not
    
    if (!isInited){
        for (size_t i=0;i<20;i++){
            names[i].insert(" C  ");
            names[i].insert(" CA ");
            names[i].insert(" CB ");
            names[i].insert(" N  ");
            names[i].insert(" O  ");
        }
        names[1].insert(" CG ");names[1].insert(" CD ");names[1].insert(" CZ ");
        names[1].insert(" NE ");names[1].insert(" NH1");names[1].insert(" NH2");
        
        names[2].insert(" CG ");names[2].insert(" ND2");names[2].insert(" OD1");
        
        names[3].insert(" CG ");names[3].insert(" OD1");names[3].insert(" OD2");
        
        names[4].insert(" SG ");
        
        names[5].insert(" CG ");names[5].insert(" CD ");
        names[5].insert(" NE2");names[5].insert(" OE1");
        
        names[6].insert(" CG ");names[6].insert(" CD ");
        names[6].insert(" OE2");names[6].insert(" OE1");
        
        names[7].erase(" CB "); // GLY
        
        names[8].insert(" CG ");names[8].insert(" CD2");names[8].insert(" CE1");
        names[8].insert(" ND1");names[8].insert(" NE2");
        
        names[9].insert(" CG1");names[9].insert(" CG2");names[9].insert(" CD1");
        
        names[10].insert(" CG ");names[10].insert(" CD2");
        names[10].insert(" CD1");// L
        
        names[11].insert(" CG ");names[11].insert(" CD ");
        names[11].insert(" CE ");names[11].insert(" NZ ");// K
        
        names[12].insert(" CG ");names[12].insert(" CE ");
        names[12].insert(" SD ");// M
        
        names[13].insert(" CG ");names[13].insert(" CD1");
        names[13].insert(" CD2");names[13].insert(" CE1");
        names[13].insert(" CE2");names[13].insert(" CZ ");// F
        
        names[14].insert(" CG ");names[14].insert(" CD ");// P
        
        names[15].insert(" OG ");// S
        
        names[16].insert(" CG2");names[16].insert(" OG1");// T
        
        names[17].insert(" CG ");names[17].insert(" CD1");
        names[17].insert(" CD2");names[17].insert(" CE2");
        names[17].insert(" CE3");names[17].insert(" CZ2");
        names[17].insert(" CZ3");names[17].insert(" CH2");
        names[17].insert(" NE1");// W
        
        names[18].insert(" CG ");names[18].insert(" CD1");
        names[18].insert(" CD2");names[18].insert(" CE1");
        names[18].insert(" CE2");names[18].insert(" CZ ");
        names[18].insert(" OH ");// Y
        
        names[19].insert(" CG1");names[19].insert(" CG2");//V
        
        isInited=true; // finished initialization here
    }
    return names[state].count(atomName);
}

////////////////////////////////////////////////////////////////////////////////
// get side-chain center and c from back bone                                 //
////////////////////////////////////////////////////////////////////////////////

//            XYZ v1=n-ca;
//            XYZ v2=c-ca;
//            XYZ v3=v1*v2;
//            v2=v1*v3;
//            
//            v1/=abs(v1); v2/=abs(v2); v3/=abs(v3);
//            
//            AAprojects[type].x+=cb.dot(v1);
//            AAprojects[type].y+=cb.dot(v2);
//            AAprojects[type].z+=cb.dot(v3);


    const float AA_CB[20][3]={
        /* A */{-0.264058,0.383927,0.603235},  /* A */
        /* R */{-1.28093, 2.23063, 2.00804 },  /* R */  
        /* N */{-0.634521,1.29839, 1.08347 },  /* N */ 
        /* D */{-0.683145,1.24647, 1.14152 },  /* D */ 
        /* C */{-0.549897,1.08546, 1.12241 },  /* C */ 
        /* Q */{-0.935605,1.73116, 1.40806 },  /* Q */ 
        /* E */{-0.937709,1.71671, 1.49664 },  /* E */ 
        /* G */{0,        0,       0       },  /* G */               
        /* H */{-0.822343,1.63728, 1.25214 },  /* H */ 
        /* I */{-0.63782, 1.19429, 1.19915 },  /* I */  
        /* L */{-0.804231,1.46024, 1.06092 },  /* L */ 
        /* K */{-1.0844,  1.91647, 1.58154 },  /* K */   
        /* M */{-1.00481, 1.76572, 1.37508 },  /* M */  
        /* F */{-0.929199,1.75043, 1.29308 },  /* F */ 
        /* P */{0.618194, 0.609921,1.08212 },  /* P */ 
        /* S */{-0.406099,0.584326,1.06051 },  /* S */
        /* T */{-0.445167,0.822609,1.12858 },  /* T */
        /* W */{-1.28204, 1.87678, 1.46564 },  /* W */  
        /* Y */{-1.0487,  1.97657, 1.38437 },  /* Y */   
        /* V */{-0.603196,0.850487,0.998482}}; /* V */

} // end of namespcae Aacid

////////////////////////////////////////////////////////////////////////////////
// the end of namespace Aacid                                                 //
////////////////////////////////////////////////////////////////////////////////

#endif  // LK_AACID_HPP
