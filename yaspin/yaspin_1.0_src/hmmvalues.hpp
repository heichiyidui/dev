////////////////////////////////////////////////////////////////////////////////
// Copyright (C) 2002  Kuang Lin kxlin@nimr.mrc.ac.uk                         //
////////////////////////////////////////////////////////////////////////////////

////////////////////////////////////////////////////////////////////////////////
// a simple class to held parameters for a HMM model                          //
////////////////////////////////////////////////////////////////////////////////

#ifndef LK_HMMVALUES_HPP
#define LK_HMMVALUES_HPP 

#include <cmath>
#include <iostream>

class HmmValues{
public:
////////////////////////////////////////////////////////////////////////////////
// constructors and the destructor                                            //
////////////////////////////////////////////////////////////////////////////////
    HmmValues(int stateNumber=0,int charcNumber=0);
    HmmValues(const HmmValues& va);
    HmmValues&operator = (const HmmValues&va);

    ~HmmValues(){destroy();}
////////////////////////////////////////////////////////////////////////////////
// input and output                                                           //
////////////////////////////////////////////////////////////////////////////////
    class SizeError{}; // exception to throw when try to init with wrong sizes
    
    int StateNumber ()const {return mStateNumber;}
    int CharcNumber ()const {return mCharcNumber;}
    
    float Ptransition   (int i,int j)const {return mPtransition[i][j];}
    float LogPtransition(int i,int j)const {return mLogPtransition[i][j];}
    float Pemission     (int i,int j)const {return mPemission[i][j];}
    float LogPemission  (int i,int j)const {return mLogPemission[i][j];}
    
    static const float LOG0=-999999.;  // a big negtive number for log(0)
    
    void Ptransition(int&i,int&j,float p){
        mPtransition[i][j]=p;
        if(p>0.)mLogPtransition[i][j]=std::log(p);
        else mLogPtransition[i][j]=LOG0;
    }
        
    void Pemission(int&i,int&j,float p){
        mPemission[i][j]=p;
        if(p>0.)mLogPemission[i][j]=std::log(p);
        else mLogPemission[i][j]=LOG0;
    }
        
////////////////////////////////////////////////////////////////////////////////
// private members                                                            //
////////////////////////////////////////////////////////////////////////////////
private:
    void init();   // alocate memory
    void destroy();// release memory
    
    int mStateNumber;
    int mCharcNumber;
    
    float ** mPtransition;
    float ** mLogPtransition;
    float ** mPemission; 
    float ** mLogPemission;
};
////////////////////////////////////////////////////////////////////////////////
// end of the class HmmValues                                                 //
////////////////////////////////////////////////////////////////////////////////

#endif // the end of LK_HMMVALUES_HPP 
