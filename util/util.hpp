////////////////////////////////////////////////////////////////////////////////
// Copyright (C) 2005-12  Kuang Lin kuanggong@gmail.com                       //
////////////////////////////////////////////////////////////////////////////////

////////////////////////////////////////////////////////////////////////////////
// Some utility functions. I guess the name could mess things up... because   //
// everybody has a util.hpp or util.h                                         //
////////////////////////////////////////////////////////////////////////////////

#ifndef __util_hpp__
#define __util_hpp__

#include <string>
#include <sstream>
#include <iostream>
#include <vector>
#include <algorithm>

using namespace std;
    
////////////////////////////////////////////////////////////////////////////////
// converting between numbers and strings                                     //
////////////////////////////////////////////////////////////////////////////////
    
inline string& toupper(string &s) 
    {transform(s.begin(),s.end(),s.begin(),(int(*)(int))toupper);return s;}

inline string& tolower(string &s)
    {transform(s.begin(),s.end(),s.begin(),(int(*)(int))tolower);return s;}

template <class T> inline T fromString(const string& s, 
                                       ios_base& (*f)(ios_base&)=std::dec)
    {istringstream is(s);T t; is >> f >> t; return t;} // f can be hex or oct

template <class T> inline bool fromString(const string& s, T& t,
                                          ios_base& (*f)(ios_base&)=std::dec)
    {istringstream is(s);return !(is >> f >> t).fail();}

template <class T> inline string toString(T t)  // very simple version, need
    {ostringstream os; os<<t; return os.str();} // different precisions here

////////////////////////////////////////////////////////////////////////////////
// The tokenizers                                                             //
////////////////////////////////////////////////////////////////////////////////

inline void strSplit(const string& str, vector <string>& tokens, 
                     const string& delimiters=" "){
    tokens.clear();
    string::size_type lastPos = str.find_first_not_of(delimiters,0);
    string::size_type pos     = str.find_first_of    (delimiters,lastPos);
    while (string::npos!=pos || string::npos!=lastPos){
        tokens.push_back(str.substr(lastPos,pos-lastPos));
        lastPos = str.find_first_not_of(delimiters,pos);
        pos     = str.find_first_of    (delimiters,lastPos);
    }
}// end of string tokenizer

inline vector <string> tokenize(const string& str, const string& delims=" "){
    vector <string> tokens; 
    string::size_type lastPos = str.find_first_not_of(delims,0);
    string::size_type pos     = str.find_first_of    (delims,lastPos);
    while (string::npos!=pos || string::npos!=lastPos){
        tokens.push_back(str.substr(lastPos,pos-lastPos));
        lastPos = str.find_first_not_of(delims,pos);
        pos     = str.find_first_of    (delims,lastPos);
    }
    return tokens;
}// end of string tokenizer 

inline void strTabSplit(const string&str, vector <string>& tokens){
    tokens.clear();
    if (str.size()==0){
        tokens.push_back("");
        return;
    }
    string ts;
    std::istringstream iss(str);
    while ( getline(iss, ts, '\t') )
        tokens.push_back(ts);
    if (str[str.size()-1]=='\t')
        tokens.push_back("");
    return;
}// end of string tokenizer with '\t'

////////////////////////////////////////////////////////////////////////////////
// output...                                                                  //
////////////////////////////////////////////////////////////////////////////////

template <class T> class outln
    {public: void operator ()(T& val) const{cout<<val<<'\n';}};

template <class T> class outsp
    {public: void operator ()(T& val) const{cout<<val<<' ' ;}};
    
template <class T> void inline outl(T val){cout<<val<<'\n';}
template <class T> void inline outs(T val){cout<<val<<' ';}

// example: 
//    for_each(stringList.begin(),stringList.end(),outl <string>);
//    for_each(stringList.begin(),stringList.end(),outln<string>());

////////////////////////////////////////////////////////////////////////////////
// getting arguments                                                          //
////////////////////////////////////////////////////////////////////////////////

template <class T> inline bool getarg2(int argc,char*argv[],string str,T&val){
    for (int i=0;i<argc-1;i++)
        if (str==string(argv[i]))
            if (fromString(string(argv[i+1]),val))
                return true;
    return false;
}
        
inline bool getarg(int argc, char*argv[],string label){
    for (int i=1;i<argc;i++)
        if (label==string(argv[i]))
            return true;
    return false;
}
// exapmle:
//     getarg2(argc,argv,"threshold",thres);
//     getarg (argc,argv,"-quiet");

#endif  // __util_hpp__
