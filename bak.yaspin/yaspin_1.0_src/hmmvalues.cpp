////////////////////////////////////////////////////////////////////////////////
// a simple class to held parameters for a HMM model                          //
////////////////////////////////////////////////////////////////////////////////

#include "hmmvalues.hpp"

////////////////////////////////////////////////////////////////////////////////
// constructors and the destructor                                            //
////////////////////////////////////////////////////////////////////////////////

HmmValues::HmmValues(int stateNumber,int charcNumber){
    if (stateNumber<0||charcNumber<0){
        std::cerr<<"Error in init HMM structure, too small \n";
        throw SizeError();
    }
    if (stateNumber>1000||charcNumber>1000){
        std::cerr<<"Error in init HMM structure, too big \n";
        throw SizeError();
    }
    
    mStateNumber=stateNumber; mCharcNumber = charcNumber;
    init();
    for (int i=0;i<mStateNumber;i++){
        for (int j=0;j<mStateNumber;j++){
            mPtransition[i][j]=1./mStateNumber;
            mLogPtransition[i][j]=std::log(mPtransition[i][j]);
        }
        for (int j=0;j<mCharcNumber;j++){
            mPemission[i][j]=1./mCharcNumber;
            mLogPemission[i][j]=std::log(mPemission[i][j]);
        }
    }
    return ;
}// end of HmmValues::HmmValues(int stateNumber,int charcNumber)

HmmValues::HmmValues(const HmmValues& va){
    mStateNumber=va.StateNumber();
    mCharcNumber=va.CharcNumber();
    init();
    for (int i=0;i<mStateNumber;i++){
        for (int j=0;j<mStateNumber;j++){
            mPtransition   [i][j]=va.Ptransition   (i,j);
            mLogPtransition[i][j]=va.LogPtransition(i,j);
        }
        for (int j=0;j<mCharcNumber;j++){
            mPemission   [i][j]=va.Pemission   (i,j);
            mLogPemission[i][j]=va.LogPemission(i,j);
        }
    }
    return;
}// end of HmmValues::HmmValues(const HmmValues& va)

HmmValues& HmmValues::operator = (const HmmValues&va){
    if (this==&va)
        return *this;
    destroy();
    mStateNumber=va.StateNumber();
    mCharcNumber=va.CharcNumber();
    init();
    for (int i=0;i<mStateNumber;i++){
        for (int j=0;j<mStateNumber;j++){
            mPtransition   [i][j]=va.Ptransition   (i,j);
            mLogPtransition[i][j]=va.LogPtransition(i,j);
        }
        for (int j=0;j<mCharcNumber;j++){
            mPemission   [i][j]=va.Pemission   (i,j);
            mLogPemission[i][j]=va.LogPemission(i,j);
        }
    }
    return *this;
}// end of HmmValues&HmmValues::operator = (const HmmValues&va)

////////////////////////////////////////////////////////////////////////////////
// private members                                                            //
////////////////////////////////////////////////////////////////////////////////

void HmmValues::init(){
    mPtransition    =new float*[mStateNumber+1];
    mLogPtransition =new float*[mStateNumber+1];
    mPemission      =new float*[mStateNumber+1];
    mLogPemission   =new float*[mStateNumber+1];
    
    for (int i=0;i<mStateNumber+1;i++){
        mPtransition   [i] =new float[mStateNumber+1];
        mLogPtransition[i] =new float[mStateNumber+1];
    }
    for (int i=0;i<mStateNumber+1;i++){
        mPemission     [i] =new float[mCharcNumber+1];
        mLogPemission  [i] =new float[mCharcNumber+1];
    }
}// end of HmmValues::init()

void HmmValues::destroy(){
    for (int i=0;i<mStateNumber+1;i++){
        delete[] mPtransition   [i];
        delete[] mLogPtransition[i];
        delete[] mPemission     [i];
        delete[] mLogPemission  [i];
    }
    delete[] mPtransition   ;
    delete[] mLogPtransition;
    delete[] mPemission     ;
    delete[] mLogPemission  ;
}// end of HmmValues::destroy()

////////////////////////////////////////////////////////////////////////////////
// end of the class HmmValues                                                 //
////////////////////////////////////////////////////////////////////////////////
