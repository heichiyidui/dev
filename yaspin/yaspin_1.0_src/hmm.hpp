////////////////////////////////////////////////////////////////////////////////
// Copyright (C) 2002  Kuang Lin kxlin@nimr.mrc.ac.uk                         //
////////////////////////////////////////////////////////////////////////////////

////////////////////////////////////////////////////////////////////////////////
// a simple implementation of some hidden Markov model algorithms             //
////////////////////////////////////////////////////////////////////////////////

#ifndef LK_HMM_HPP
#define LK_HMM_HPP

#include "hmmvalues.hpp"

////////////////////////////////////////////////////////////////////////////////
// begin of the class Hmm                                                     //
////////////////////////////////////////////////////////////////////////////////
class Hmm{
public:
    HmmValues mCurrent;
    HmmValues mReestimate;
////////////////////////////////////////////////////////////////////////////////
// constructors and the destructor                                            //
////////////////////////////////////////////////////////////////////////////////
    Hmm(HmmValues va=HmmValues(),int seqLength=0);
    Hmm(const Hmm& h);
    Hmm&operator = (const Hmm&h);
    ~Hmm(){destroy();}
////////////////////////////////////////////////////////////////////////////////
// output and input functions for class Hmm                                   //
////////////////////////////////////////////////////////////////////////////////
    class SizeError{}; // the exception report that some initializer is wrong
    
    int SeqLength  ()const{ return mSeqLength  ;}
    
    int   BestPath  (int i)const{ return mBestPath[i];}
    
    float HmmF(int i,int j)const{ return mHmmF[i][j];}
    float HmmB(int i,int j)const{ return mHmmB[i][j];}
    float HmmV(int i,int j)const{ return mHmmV[i][j];}
    
    float HmmSf(int i) const {return mHmmSf[i];}
    float HmmSb(int i) const {return mHmmSb[i];}
    
    void SeqLength (int seqLength);
////////////////////////////////////////////////////////////////////////////////
// algorithms for class Hmm                                                   //
////////////////////////////////////////////////////////////////////////////////
    float Viterbi (int *seq);// with log of Probabilities
    float forward (int *seq);// with Scaling of probabilities
    float backward(int *seq);// with Scaling of probabilities
    
    float Viterbi (float**pemission);// 
    float forward (float**pemission);
    float backward(float**pemission);
////////////////////////////////////////////////////////////////////////////////
// private members                                                            //
////////////////////////////////////////////////////////////////////////////////
private:
    
    void init();   // alocate memory
    void destroy();// release memory
    
    int mSeqLength;
    
    float*  mHmmSf;
    float*  mHmmSb;
    
    float** mHmmF;
    float** mHmmB;
    float** mHmmV;
    
    int*   mBestPath;
    
    int**  mTrack;
    
    float* mHmmFs; // these pointers are not directly accessable 
    float* mHmmBs; // they will be used via mHmmF etc.
    float* mHmmVs; // a benefit is we don't need mSeqLength for destroy()
    int  * mTracks;
};
////////////////////////////////////////////////////////////////////////////////
// end of class Hmm                                                           //
////////////////////////////////////////////////////////////////////////////////

#endif // LK_HMM_HPP
