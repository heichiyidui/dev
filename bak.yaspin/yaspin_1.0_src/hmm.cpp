////////////////////////////////////////////////////////////////////////////////
// Copyright (C) 2002  Kuang Lin kxlin@nimr.mrc.ac.uk                         //
////////////////////////////////////////////////////////////////////////////////

#include "hmm.hpp"

////////////////////////////////////////////////////////////////////////////////
// constructures                                                              //
////////////////////////////////////////////////////////////////////////////////
Hmm::Hmm(HmmValues va,int seqLength){
    mCurrent=va;
    if (seqLength<0)  throw SizeError();
    mSeqLength=seqLength; init();
}// enf of Hmm::Hmm(HmmValues,int)

Hmm::Hmm(const Hmm& h){
    mCurrent =h.mCurrent; mReestimate = h.mReestimate;
    mSeqLength=h.SeqLength();
    init();
    for (int i=0;i<mSeqLength+1;i++){
        mBestPath[i] = h.BestPath(i);
        mHmmSf   [i] = h.HmmSf(i);
        mHmmSb   [i] = h.HmmSb(i);
        for (int j=0;j<(mCurrent.StateNumber()+1);j++){
            mHmmF[i][j] = h.HmmF(i,j);
            mHmmB[i][j] = h.HmmB(i,j);
            mHmmV[i][j] = h.HmmV(i,j);
        }
    }
}// end of Hmm::Hmm(Hmm&)

Hmm&Hmm::operator = (const Hmm&h){
    if (this==&h)
        return *this;
    destroy();
    mCurrent = h.mCurrent; mReestimate=h.mReestimate;
    mSeqLength  =h.SeqLength();
    init(); 
    for (int i=0;i<mSeqLength+1;i++){
        mBestPath[i] = h.BestPath(i);
        mHmmSf   [i] = h.HmmSf(i);
        mHmmSb   [i] = h.HmmSb(i);
        for (int j=0;j<(mCurrent.StateNumber()+1);j++){
            mHmmF[i][j] = h.HmmF(i,j);
            mHmmB[i][j] = h.HmmB(i,j);
            mHmmV[i][j] = h.HmmV(i,j); // copy are incomplete, tracks etc. left
        }
    }
    return *this;
}// end of Hmm&Hmm::operator = (const Hmm&h)

////////////////////////////////////////////////////////////////////////////////
// input functions for class Hmm                                              //
////////////////////////////////////////////////////////////////////////////////
void Hmm::SeqLength (int seqLength){
    if (seqLength<0)throw SizeError();
    destroy();
    mSeqLength=seqLength;
    init();
    return;
}
////////////////////////////////////////////////////////////////////////////////
// algorithms for class Hmm                                                   //
////////////////////////////////////////////////////////////////////////////////
    
float Hmm::Viterbi(int *seq){ // the Viterbi algorithm
    // init
    mHmmV[0][0]=0.; // log(1)
    for (int i=1;i<mCurrent.StateNumber();i++)
        mHmmV[0][i]=HmmValues::LOG0;
    
    // recursion
    for (int i=1;i<=mSeqLength;i++){
        for (int j=0;j<mCurrent.StateNumber();j++) mHmmV[i][j]=0.;
        
        for (int j=0;j<mCurrent.StateNumber();j++){
            float max=HmmValues::LOG0;
            for (int k=0;k<mCurrent.StateNumber();k++){
                if (max < mHmmV[i-1][k] + mCurrent.LogPtransition(k,j)){
                    max = mHmmV[i-1][k] + mCurrent.LogPtransition(k,j);
                    mTrack[i][j]=k;
                }
            }
            mHmmV[i][j]=max + mCurrent.LogPemission(j,seq[i-1]);
        }
    }
    
    // termination
    float max=HmmValues::LOG0;
    int LastI=0;
    for (int i=1;i<mCurrent.StateNumber();i++){
        if (max < mHmmV[mSeqLength][i] + mCurrent.LogPtransition(i,0)){
            max = mHmmV[mSeqLength][i] + mCurrent.LogPtransition(i,0);
            LastI=i;
        }
    }
    
    // trackback
    if (LastI==0){
        std::cout<<"failed to terminate in Hmm::Viterbi()\n";
        throw SizeError();
    }
    for (int i=mSeqLength;i>=1;i--){
        mBestPath[i-1]=LastI;
        LastI=mTrack[i][LastI];
    }
    return max; // found in termination
}// end of float Hmm::Viterbi(int *seq)


float Hmm::Viterbi(float**pemission){
    using namespace std;
    // init
    mHmmV[0][0]=0; // log(1)
    for (int i=1;i<mCurrent.StateNumber();i++)
        mHmmV[0][i]=HmmValues::LOG0;
    
    // recursion
    for (int i=1;i<=mSeqLength;i++){
        
        for (int j=0;j<mCurrent.StateNumber();j++){
            float max=HmmValues::LOG0;
            for (int k=0;k<mCurrent.StateNumber();k++){
                if (max < mHmmV[i-1][k] + mCurrent.LogPtransition(k,j)){
                    max = mHmmV[i-1][k] + mCurrent.LogPtransition(k,j);
                    mTrack[i][j]=k;
                }
            }
            if (pemission[i-1][j]!=0)
                mHmmV[i][j]=max+log(pemission[i-1][j]);
            else
                mHmmV[i][j]=max+HmmValues::LOG0;
        }
    }
    
    // termination
    float max=HmmValues::LOG0;
    int LastI=0;
    for (int i=1;i<mCurrent.StateNumber();i++){
        if (max < mHmmV[mSeqLength][i] + mCurrent.LogPtransition(i,0)){
            max = mHmmV[mSeqLength][i] + mCurrent.LogPtransition(i,0);
            LastI=i;
        }
    }
    // trackback
    if (LastI==0){
        std::cerr<<"failed to terminate in Hmm::Viterbi()\n";
        throw SizeError();
    }
    for (int i=mSeqLength;i>=1;i--){
        mBestPath[i-1]=LastI;
        LastI=mTrack[i][LastI];
    }
    
    return max; // found in termination
}// end of float Hmm::Viterbi(float**pemission)

float Hmm::forward(int* seq){
    // init
    mHmmF[0][0]=1;
    for (int i=1;i<mCurrent.StateNumber();i++)mHmmF[0][i]=0;
    
    // recursion
    for (int i=1;i<=mSeqLength;i++){
        for (int j=0;j<mCurrent.StateNumber();j++) mHmmF[i][j]=0;
        
        for (int j=0;j<mCurrent.StateNumber();j++)
            for (int k=0;k<mCurrent.StateNumber();k++)
                mHmmF[i][j]+=mCurrent.Pemission(j,seq[i-1])
                            *mCurrent.Ptransition(k,j)
                            *mHmmF[i-1][k];
        mHmmSf[i]=0;
        for (int j=0;j<mCurrent.StateNumber();j++)
            mHmmSf[i]+=mHmmF[i][j];
        for (int j=0;j<mCurrent.StateNumber();j++)
            mHmmF[i][j]/=mHmmSf[i];
    }
    
    // termination
    float LogP=0;
    for (int i=0;i<mCurrent.StateNumber();i++){
        LogP+=mHmmF[mSeqLength][i]*mCurrent.Ptransition(i,0);
    }
    LogP=log(LogP);
    
    for (int i=1;i<=mSeqLength;i++)
        LogP+=log(mHmmSf[i]);
    return LogP;
}// end of float Hmm::forward(int* seq)

float Hmm::forward(float**pemission){
    // init
    mHmmF[0][0]=1;
    for (int i=1;i<mCurrent.StateNumber();i++)mHmmF[0][i]=0;
    
    // recursion
    for (int i=1;i<=mSeqLength;i++){
        for (int j=0;j<mCurrent.StateNumber();j++) mHmmF[i][j]=0;
        
        for (int j=0;j<mCurrent.StateNumber();j++){
            for (int k=0;k<mCurrent.StateNumber();k++)
                mHmmF[i][j]+=mCurrent.Ptransition(k,j)*mHmmF[i-1][k];
        
            mHmmF[i][j]*=pemission[i-1][j];
        }
        
        mHmmSf[i]=0;
        for (int j=0;j<mCurrent.StateNumber();j++)
            mHmmSf[i]+=mHmmF[i][j];
        for (int j=0;j<mCurrent.StateNumber();j++)
            mHmmF[i][j]/=mHmmSf[i];
    }
    
    
    // termination
    float LogP=0.;
    for (int i=0;i<mCurrent.StateNumber();i++){
        LogP+=mHmmF[mSeqLength][i]*mCurrent.Ptransition(i,0);
    }
    LogP=log(LogP);
    
    for (int i=1;i<=mSeqLength;i++)
        LogP+=log(mHmmSf[i]);
    return LogP;
}// end of float Hmm::forward(float**pemission)

float Hmm::backward(int* seq){ 
    // init
    mHmmSb[mSeqLength]=0.;
    for (int i=0;i<mCurrent.StateNumber();i++){
        
        mHmmB[mSeqLength][i]=mCurrent.Pemission(i,seq[mSeqLength-1])
                            *mCurrent.Ptransition(i,0);
        mHmmSb[mSeqLength] += mHmmB[mSeqLength][i];
    }
    for (int i=0;i<mCurrent.StateNumber();i++)
        mHmmB[mSeqLength][i]/=mHmmSb[mSeqLength];
    
    // recursion
    for (int i=mSeqLength-1;i>0;i--){
        for (int j=0;j<mCurrent.StateNumber();j++)
            mHmmB[i][j]=0.;
        
        for (int j=0;j<mCurrent.StateNumber();j++)
            for (int k=0;k<mCurrent.StateNumber();k++)
                mHmmB[i][j]+=mCurrent.Pemission(j,seq[i-1])
                            *mCurrent.Ptransition(j,k)
                            *mHmmB[i+1][k];
        mHmmSb[i]=0;
        for (int j=0;j<mCurrent.StateNumber();j++)
            mHmmSb[i]+=mHmmB[i][j];
        for (int j=0;j<mCurrent.StateNumber();j++)
            mHmmB[i][j]/=mHmmSb[i];
    }
    // termination
    float LogP=0;
    for (int i=0;i<mCurrent.StateNumber();i++){
        LogP+=mCurrent.Ptransition(0,i)*mHmmB[1][i];
    }
    LogP=log(LogP);
    
    for (int i=1;i<=mSeqLength;i++)
        LogP+=log(mHmmSb[i]);
    return LogP;
}// end of float Hmm::backward(int* seq)

float Hmm::backward(float** pemission){ 
    // init
    mHmmSb[mSeqLength]=0.;
    for (int i=0;i<mCurrent.StateNumber();i++){
        mHmmB[mSeqLength][i]=pemission[mSeqLength-1][i]
                            *mCurrent.Ptransition(i,0);
        
        mHmmSb[mSeqLength] += mHmmB[mSeqLength][i];
    }
    
    for (int i=0;i<mCurrent.StateNumber();i++)
        mHmmB[mSeqLength][i]/=mHmmSb[mSeqLength];
    
    // recursion
    for (int i=mSeqLength-1;i>0;i--){
        for (int j=0;j<mCurrent.StateNumber();j++)
            mHmmB[i][j]=0.;
        
        for (int j=0;j<mCurrent.StateNumber();j++)
            for (int k=0;k<mCurrent.StateNumber();k++)
                mHmmB[i][j]+=pemission[i-1][j]
                            *mCurrent.Ptransition(j,k)
                            *mHmmB[i+1][k];
        mHmmSb[i]=0;
        for (int j=0;j<mCurrent.StateNumber();j++)
            mHmmSb[i]+=mHmmB[i][j];
        for (int j=0;j<mCurrent.StateNumber();j++)
            mHmmB[i][j]/=mHmmSb[i];
    }
    // termination
    float LogP=0;
    for (int i=0;i<mCurrent.StateNumber();i++){
        LogP+=mCurrent.Ptransition(0,i)*mHmmB[1][i];
    }
    LogP=log(LogP);
    
    for (int i=1;i<=mSeqLength;i++)
        LogP+=log(mHmmSb[i]);
    return LogP;
}// end of float Hmm::backward(float** pemission)

////////////////////////////////////////////////////////////////////////////////
// private members                                                            //
////////////////////////////////////////////////////////////////////////////////

void Hmm::init(){
    mBestPath=new int  [mSeqLength+1];
    mTrack   =new int* [mSeqLength+1];
    
    mHmmV = new float * [mSeqLength+1];
    mHmmB = new float * [mSeqLength+1];
    mHmmF = new float * [mSeqLength+1];
    
    mTracks= new int [(mSeqLength+1)*mCurrent.StateNumber()];
    mHmmVs = new float [(mSeqLength+1)*(mCurrent.StateNumber()+1)];
    mHmmBs = new float [(mSeqLength+1)*(mCurrent.StateNumber()+1)];
    mHmmFs = new float [(mSeqLength+1)*(mCurrent.StateNumber()+1)];
    
    if (mHmmFs==NULL) throw SizeError();
    
    for (int i=0;i<mSeqLength+1;i++){
        mTrack[i] = & mTracks[i*mCurrent.StateNumber()];
        mHmmV [i] = & mHmmVs[i*(mCurrent.StateNumber()+1)];
        mHmmB [i] = & mHmmBs[i*(mCurrent.StateNumber()+1)];
        mHmmF [i] = & mHmmFs[i*(mCurrent.StateNumber()+1)];
    }
    mHmmSf = new float [mSeqLength+1];
    mHmmSb = new float [mSeqLength+1];
}// end of void Hmm::init()

void Hmm::destroy(){
    delete[] mBestPath;delete[] mTrack;   delete[] mTracks;
    delete[] mHmmV;    delete[] mHmmF;    delete[] mHmmB;
    delete[] mHmmVs;   delete[] mHmmFs;   delete[] mHmmBs;
    delete[] mHmmSf;   delete[] mHmmSb;
}// end of void Hmm::destroy()
////////////////////////////////////////////////////////////////////////////////
// end of the class Hmm                                                       //
////////////////////////////////////////////////////////////////////////////////

