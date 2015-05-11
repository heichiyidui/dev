////////////////////////////////////////////////////////////////////////////////
// Copyright (C) 2009  Kuang Lin kuanggong@gmail.com                          //
////////////////////////////////////////////////////////////////////////////////

////////////////////////////////////////////////////////////////////////////////
//                                                                            //
// I need some subroutines to compute the Delaunay triangulation, (Delaunay   //
// tetrahedralization) of a set of points in 3D space. In fact, I need the    //
// neighbour lists of atoms or residues of proteins.                          //
//                                                                            //
// If the TetGen package by Si Hang can get the job done, so be it.           //
//                                                                            //
// TetGen License                                                             //
// --------------                                                             //
//                                                                            //
// The software (TetGen) is licensed under the terms of the  MIT  license     //
// with the following exceptions:                                             //
//                                                                            //
// Distribution of  modified  versions  of this code is permissible UNDER     //
// THE CONDITION THAT  THIS CODE AND ANY MODIFICATIONS  MADE TO IT IN THE     //
// SAME SOURCE FILES  tetgen.h AND tetgen.cxx  REMAIN UNDER  COPYRIGHT OF     //
// THE  ORIGINAL AUTHOR,  BOTH  SOURCE AND OBJECT  CODE  ARE MADE  FREELY     //
// AVAILABLE  WITHOUT   CHARGE,   AND  CLEAR   NOTICE  IS  GIVEN  OF  THE     //
// MODIFICATIONS.                                                             //
//                                                                            //
// Distribution of this code for  any  commercial purpose  is permissible     //
// ONLY BY DIRECT ARRANGEMENT WITH THE COPYRIGHT OWNER.                       //
//                                                                            //
// The full  license text is reproduced below.                                //
//                                                                            //
// This means that TetGen is no free software, but for private, research,     //
// and  educational purposes it  can be  used at  absolutely no  cost and     //
// without further arrangements.                                              //
//                                                                            //
//                                                                            //
// For details, see http://tetgen.berlios.de                                  //
//                                                                            //
// And of course, the predicates.cxx is from Jonathan Richard Shewchuk.       //
// They are Routines for Arbitrary Precision Floating-point Arithmetic        //
// and Fast Robust Geometric Predicates.                                      //
//                                                                            //
////////////////////////////////////////////////////////////////////////////////

////////////////////////////////////////////////////////////////////////////////
//                                                                            //
// For a set of points S in space R^d (here three-dimentional space R^3),     //
// each point s has a Voronoi cell, also called Dirichlet cell, V(s)          //
// consisting of all points closer to s than any other points in S.           //
//                                                                            //
// Now the Delaunay triangulation (tetrahedralization in the case of three-   //
// dimentional space) of S, DT(S), is the tetrahedralization such that no     //
// point in S is inside the circumscribed spheres of the tetrahedron.         //
//                                                                            //
// But basically, you connect the points with common Voronoi cell walls...    //
//                                                                            //
// It is a well defined measure of closeness, used for protein structure in   //
// different applications.                                                    //
//                                                                            //
////////////////////////////////////////////////////////////////////////////////

#ifndef __COMPUTE_TETR_HPP__
#define __COMPUTE_TETR_HPP__

#include <vector>
#include <set>
#include <iostream>

#include "xyz.hpp"
#include "tetgen.h"

//g++ -c -O0 predicates.cxx
//g++ -c -O3 tetgen.cxx

////////////////////////////////////////////////////////////////////////////////
//                                                                            //
// What options we need to check here? ... No, we don't want any parameters.  //
// A vector of XYZ in, a map of tetrahedra connections out.                   //
//                                                                            //
////////////////////////////////////////////////////////////////////////////////

namespace cao
{
    using namespace std;
    
    inline vector <pair<size_t, size_t> >computeTetra(const vector<XYZ>& xyzs){
    
        tetgenio in,out;
        
        in.numberofpoints=xyzs.size();
        in.pointlist = new REAL[in.numberofpoints * 3];
    
        for (int i=0;i< in.numberofpoints;i++){
            in.pointlist[i*3  ]=xyzs[i].x;
            in.pointlist[i*3+1]=xyzs[i].y;
            in.pointlist[i*3+2]=xyzs[i].z;
        }
        
        char a[]="Q"; // "Q" for the quiet mode
        tetrahedralize(a, &in, &out); // get tetrahedra
        
        vector <pair<size_t,size_t> > output;
        
//        cout<<"out.numberoftetrahedra = "<<out.numberoftetrahedra<<'\n';
        
        for (int i=0;i<out.numberoftetrahedra;i++){
            size_t ti,tj;
            ti=out.tetrahedronlist[i*4  ]; tj=out.tetrahedronlist[i*4+1];
            output.push_back(pair<size_t,size_t> (ti,tj));
            output.push_back(pair<size_t,size_t> (tj,ti));
            ti=out.tetrahedronlist[i*4  ]; tj=out.tetrahedronlist[i*4+2];
            output.push_back(pair<size_t,size_t> (ti,tj));
            output.push_back(pair<size_t,size_t> (tj,ti));
            ti=out.tetrahedronlist[i*4  ]; tj=out.tetrahedronlist[i*4+3];
            output.push_back(pair<size_t,size_t> (ti,tj));
            output.push_back(pair<size_t,size_t> (tj,ti));
            
            ti=out.tetrahedronlist[i*4+1]; tj=out.tetrahedronlist[i*4+2];
            output.push_back(pair<size_t,size_t> (ti,tj));
            output.push_back(pair<size_t,size_t> (tj,ti));
            ti=out.tetrahedronlist[i*4+1]; tj=out.tetrahedronlist[i*4+3];
            output.push_back(pair<size_t,size_t> (ti,tj));
            output.push_back(pair<size_t,size_t> (tj,ti));
            
            ti=out.tetrahedronlist[i*4+2]; tj=out.tetrahedronlist[i*4+3];
            output.push_back(pair<size_t,size_t> (ti,tj));
            output.push_back(pair<size_t,size_t> (tj,ti));
        }
    
        set <pair<size_t,size_t> > tset(output.begin(),output.end());
        vector <pair<size_t,size_t> > output2(tset.begin(),tset.end());
        return output2;
    }
} // end of namespace cao

#endif // __COMPUTE_TETR_HPP__

