////////////////////////////////////////////////////////////////////////////////
// Copyright (C) 2002-10  Kuang Lin kuanggong@gmail.com                       //
////////////////////////////////////////////////////////////////////////////////

////////////////////////////////////////////////////////////////////////////////
// XYZ class, a class of 3-dimensional vectors for operations on coordinates  //
////////////////////////////////////////////////////////////////////////////////

#ifndef __LK_XYZ_HPP__
#define __LK_XYZ_HPP__

#include <cmath>
#include <vector>

namespace cao
{
    
    class XYZ{
        
        public:
            
            XYZ(float a=0.0F,float b=0.0F, float c=0.0F,size_t t=0)
                :x(a),y(b),z(c),type(t){}
            
            XYZ(const XYZ & v)
                :x(v.x),y(v.y),z(v.z),type(v.type){}
            
            XYZ & operator = (const XYZ & v){if(this == &v) return *this;
                x=v.x; y=v.y; z=v.z; type=v.type; return *this;}
                
            float x,y,z;  // atom coordinates
            
            size_t type;  // atom type
                          // OK, for most geometry computations, "type" is
                          // useless. But why should I bother with our small 
                          // problems and our current computting power? 
            
            static const float tolerance = 0.001; 
                // In PDB, coordinates are described with precission up
                // to 0.001 Angstrom. Any distance less than this will be 
                // considered zero. 
            
            ////////////////////////////////////////////////////////////////////
            // operators                                                      //
            ////////////////////////////////////////////////////////////////////
            
            bool operator == (const XYZ & v) const 
                {return (*this-v).abs() <  tolerance;}
                
            bool operator != (const XYZ & v) const 
                {return !(*this==v);}
            
            bool operator <  (const XYZ & v) const
                {return this->abs() < v.abs();}
            
            XYZ operator - () const {return XYZ(-x,-y,-z,type);}
            
            XYZ operator + (const XYZ&v)const{return XYZ(x+v.x, y+v.y, z+v.z);}
            XYZ operator - (const XYZ&v)const{return XYZ(x-v.x, y-v.y, z-v.z);}
            
            XYZ operator * (float s) const {return XYZ(x*s, y*s, z*s, type);}
            
            XYZ operator * (const XYZ &v)const 
                {return XYZ( y*v.z-z*v.y, z*v.x-x*v.z, x*v.y-y*v.x );}
                // the cross product 
            
            XYZ operator / (float s) const {s=1.0F/s; 
                return XYZ(x*s, y*s, z*s, type);} 
                // note: NO checking of s being zero!
            
            XYZ& operator += (const XYZ& v){x+=v.x; y+=v.y; z+=v.z;return*this;}
            
            XYZ& operator -= (const XYZ& v){x-=v.x; y-=v.y; z-=v.z;return*this;}
            
            XYZ& operator *= (float s){x*=s; y*=s; z*=s; return*this;}
            
            XYZ& operator /= (float s){s=1.0F/s; // again, 
                x*=s; y*=s; z*=s; return*this;}  // note: NO checking of zero 
            
            ////////////////////////////////////////////////////////////////////
            // member functions                                               //
            ////////////////////////////////////////////////////////////////////
            
            void  zero() { x=0.0F; y=0.0F; z=0.0F; type=0; }
                // zero the vector
            
            float absSq() const{ return x*x + y*y + z*z ; }
                // the squred absolute value, faster without the sqrt() call
                // enough for comparison of magnitude. 
            
            float abs()   const{ return sqrt(x*x + y*y + z*z); }
                // the absolute value, the magnitude of the vector
            
            void normalize(){
                float s=abs();
                if (s > 0.0F) {s=1.0F/s; x*=s; y*=s; z*=s;}}
                // normalize the vector, no error message on zero vectors.
            
            float dot(const XYZ &v )const {return x*v.x + y*v.y + z*v.z;}
                // the dot product
            
            XYZ cross(const XYZ &v) const 
                {return XYZ( y*v.z-z*v.y, z*v.x-x*v.z, x*v.y-y*v.x );}
                // the cross product 
            
    }; // end of class XYZ  
    // well, this class feels like a struct   
    
    ////////////////////////////////////////////////////////////////////////////
    // non-member functions                                                   //
    ////////////////////////////////////////////////////////////////////////////
    
    inline float abs  (const XYZ &v)
        {return sqrt(v.x*v.x + v.y*v.y + v.z*v.z);}
        // the magnitude of a vector
    
    inline float absSq(const XYZ &v)
        {return v.x*v.x + v.y*v.y + v.z*v.z;}
        // the squared value of the magnitude
    
} // end of namespace cao

#endif // the end of __LK_XYZ_HPP__

