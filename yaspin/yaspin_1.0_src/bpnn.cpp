////////////////////////////////////////////////////////////////////////////////
// Copyright (C) 2002 Kuang Lin  kxlin@nimr.mrc.ac.uk                         //
////////////////////////////////////////////////////////////////////////////////

////////////////////////////////////////////////////////////////////////////////
// This is the implement file of class BpNn, a neural network using the       //
// standard sigmoid function or the soft max transition function, supervised- //
// trainned using the backpropagation algorithm.                              //
////////////////////////////////////////////////////////////////////////////////

#include "bpnn.hpp"

////////////////////////////////////////////////////////////////////////////////
// constructors and the destructor                                            //
////////////////////////////////////////////////////////////////////////////////
BpNn::BpNn(const BpNn&nn){              // copy initialize
    init(0,0,0,true);
    *this=nn;
}// end of BpNn::BpNn(const BpNn&)

BpNn& BpNn::operator = (const BpNn&nn){ // copy assigner
    if (&nn==this)
        return *this;
    destroy();
    init(nn.NetI(),nn.NetM(),nn.NetO());
    NetAlpha(nn.NetAlpha());
    NetEta  (nn.NetEta());
    NetBias (nn.NetBias());
    NetGain (nn.NetGain());
    
    float *tfp=new float[nn.NetSize()];
    nn.OutputWeights(tfp);
    InputWeights(tfp);
    delete[]tfp;
    return *this;
}// enf of BpNn& BpNn::operator = (const BpNn&nn)

BpNn::~BpNn(){ 
    destroy();
    if (mTrainSize>0){
        delete[]mTrainBuffer;
        delete[]mInputBuffer;
        delete[]mTargetBuffer;
    }
    return;
}// end of ~BpNn(void)


////////////////////////////////////////////////////////////////////////////////
//  Input parameters and neural network                                       //
////////////////////////////////////////////////////////////////////////////////
BpNn& BpNn::NetAlpha(float netAlpha){
    if (netAlpha<0)  throw SizeError();
    mNetAlpha=netAlpha; 
    return *this;
}// end of NetAlpha(float)

BpNn& BpNn::NetEta(float netEta){
    if (netEta<=0)  throw SizeError();
    mNetEta=netEta;
    return *this;
}// end of NetEta(float)

BpNn& BpNn::NetBias(float netBias){
    mNetBias=netBias;
    return *this;
}// end of NetBias()

BpNn& BpNn::NetGain(float netGain){
    mNetGain=netGain;
    return *this;
}// end of NetGain()

BpNn& BpNn::CrossFold(int crossFold){
    if (crossFold<2) throw SizeError();
    mCrossFold = crossFold;
    return *this;
}// end of CrossFold(int)

void BpNn::RandomWeights(float weightRange){
    if (mNetI<0||mNetM<0||mNetO<=0){
        throw SizeError();
    }
    if (weightRange<=0){
        throw SizeError();
        weightRange=fabs(weightRange);
    }
    for (int i=1;i<mNetM+1;i++)
        for (int j=0;j<mNetI+1;j++)
            mWeight[0][i][j] = (float(random())/RAND_MAX)
                              *2*weightRange - weightRange;
    for (int i=1;i<mNetO+1;i++)
        for (int j=0;j<mNetM+1;j++)
            mWeight[1][i][j] = (float(random())/RAND_MAX)
                              *2*weightRange - weightRange;
}// end of RandomWeights(float)

void BpNn::InputNn(const char* nnFileName){
    using namespace std;
    ifstream NnFile(nnFileName);
    if (!NnFile.is_open()){
        cerr<<"Failed to open the nn file: "<<nnFileName<<'\n';
        throw SizeError();
    }

    char TempString[255];
    int NetI,NetM,NetO,TempInt;
    float TempFloat;
    for (;NnFile>>TempString;){
        if (0==strcmp(TempString,"NetI"))
            NnFile>>NetI;
        else if (0==strcmp(TempString,"NetM"))
            NnFile>>NetM;
        else if (0==strcmp(TempString,"NetO")){
            NnFile>>NetO; 
        }
        else if (0==strcmp(TempString,"NetAlpha")){
            NnFile>>TempFloat;
            NetAlpha(TempFloat);
        }
        else if (0==strcmp(TempString,"NetEta")){
            NnFile>>TempFloat;
            NetEta(TempFloat);
        }
        else if (0==strcmp(TempString,"NetBias")){
            NnFile>>TempFloat;
            NetBias(TempFloat);
        }
        else if (0==strcmp(TempString,"NetGain")){
            NnFile>>TempFloat;
            NetGain(TempFloat);
        }
        else if (0==strcmp(TempString,"NetCrossFold")){
            NnFile>>TempInt;
            CrossFold(TempInt);
        }
        else if (0==strcmp(TempString,"Weights"))
            break;
        else{
            throw SizeError(); // if we find some strange thing
        }
    }
    if (NnFile.eof())
        throw SizeError(); // we can't find the weigths
    
    init(NetI,NetM,NetO);// set the structure
    
    for (int i=0;i<(mNetI+1)*(mNetM+1)+(mNetM+1)*(mNetO+1);i++)
        NnFile>>mWeights[i];
    
    NnFile.close();
    return;
}// end of InputNn()

void BpNn::InputWeights(const float* weights){
    for (int i=1;i<mNetM+1;i++)
        for (int j=0;j<mNetI+1;j++)
            mWeight[0][i][j]=weights[(i-1)*(mNetI+1)+j];
    for (int i=1;i<mNetO+1;i++)
        for (int j=0;j<mNetM+1;j++)
            mWeight[1][i][j]
                =weights[(mNetI+1)*mNetM+(i-1)*(mNetM+1)+j];
}// end of InputWeights(float*)
 
void BpNn::InputTrainSet(const char*trainFileName,int trainSize){
    if (trainSize<mCrossFold)
        throw SizeError();
    
    if (mTrainSize>0){
        delete[]mTrainBuffer;
        delete[]mInputBuffer;
        delete[]mTargetBuffer;
        mTrainSize=0;
    }
    
    std::ifstream InTrainRes(trainFileName);
    if (!InTrainRes.is_open()){
        std::cerr<<"Failed to read a training set from "<<trainFileName<<"!\n";
        throw SizeError();
    }

    mTrainSize=trainSize;
    mInputBuffer =new float *[trainSize];
    mTargetBuffer=new float *[trainSize];
    mTrainBuffer =new float [trainSize*(mNetI+mNetO)];
    if (NULL==mTrainBuffer){
        std::cerr<<"Failed to allocate enough memory for the training set \n";
        throw SizeError();
    }
    
    for (int l=0;l<mTrainSize;l++){
        mInputBuffer [l] = & mTrainBuffer[l*(mNetI+mNetO)];
        mTargetBuffer[l] = & mTrainBuffer[l*(mNetI+mNetO)+mNetI];
    }

    std::cout<<"Reading train data from the file "<<trainFileName<<"...\n";
    char ts[555]; // a temp string
    for (int l=0;l<mTrainSize;l++){
        for (int i=0;i<mNetI;i++)
            InTrainRes>>mInputBuffer[l][i];
            InTrainRes.getline(ts,555);
        for (int i=0;i<mNetO;i++)
            InTrainRes>>mTargetBuffer[l][i];
            InTrainRes.getline(ts,555);
    }

    InTrainRes>>ts;
    std::cout<<"Read TrainFile last word: "<<ts<<'\n';
    InTrainRes.close();
}// end of InputTrainSet(char*,int)

////////////////////////////////////////////////////////////////////////////////
//  Output                                                                    //
////////////////////////////////////////////////////////////////////////////////
void BpNn::OutputNn(const char* nnFileName) const {
    std::ofstream NnFile(nnFileName);
    if (!NnFile.is_open()){
        std::cerr<<"Failed to open the output file "<<nnFileName<<'\n';
        throw SizeError();
    }
    
    NnFile<<"NetI     "<<mNetI<<"  \n"
          <<"NetM     "<<mNetM<<"  \n"
          <<"NetO     "<<mNetO<<" \n\n"
          <<"NetAlpha "<<mNetAlpha<<"  \n"
          <<"NetEta   "<<mNetEta  <<"  \n"
          <<"NetBias  "<<mNetBias <<"  \n"
          <<"NetGain  "<<mNetGain <<"  \n"
          <<"NetCrossFold "<<mCrossFold<<" \n\n";

    NnFile<<"Weights  \n";
    int i,j;
    for (i=0;i<mNetM+1;i++)for (j=0;j<mNetI+1;j++)NnFile<<mWeight[0][i][j]<<' ';
    NnFile<<'\n';
    for (i=0;i<mNetO+1;i++)for (j=0;j<mNetM+1;j++)NnFile<<mWeight[1][i][j]<<' ';
    NnFile<<'\n';

    NnFile.close();
    return;
}// end of OutputNn(char*)

void BpNn::OutputNn(void)const{
    std::cout<<"\nNetI     "<<mNetI<<"  \n"
             <<"NetM     "<<mNetM<<"  \n"
             <<"NetO     "<<mNetO<<"  \n\n"
             <<"NetAlpha "<<mNetAlpha<<"  \n"
             <<"NetEta   "<<mNetEta  <<"  \n"
             <<"NetBias  "<<mNetBias <<"  \n"
             <<"NetGain  "<<mNetGain <<"  \n"
             <<"NetCrossFold "<<mCrossFold <<"  \n\n";
}// end of OutputNn() 

void BpNn::OutputWeights(float* weights)const{
    int i,j;
    for (i=1;i<mNetM+1;i++)
        for (j=0;j<mNetI+1;j++)
    	    weights[(i-1)*(mNetI+1)+j]=mWeight[0][i][j];
            
    for (i=1;i<mNetO+1;i++)
        for (j=0;j<mNetM+1;j++)
    	    weights[(mNetI+1)*mNetM+(i-1)*(mNetM+1)+j]
                    =mWeight[1][i][j];
}// end of OutputWeights(float*)

int BpNn::NetSize(void) const{
    return (mNetI+1)*mNetM+(mNetM+1)*mNetO;
}// end of NetSize()

////////////////////////////////////////////////////////////////////////////////
//  Propagate, train and test of Net                                          //
////////////////////////////////////////////////////////////////////////////////
//1. with the sigmoid transition function

float const* BpNn::Propagate(const float* inputData){
    register int i,j;
    register float sum;
    for (i=1;i<mNetM+1;i++){
        sum = mNetBias*mWeight[0][i][0];
        for (j=1;j<mNetI+1;j++)
            sum += inputData[j-1]*mWeight[0][i][j];
        mOutput[0][i]=1.0/(1.0+exp(-mNetGain*sum));//the sigmoid function
    }
    for (i=1;i<mNetO+1;i++){
        sum = mNetBias*mWeight[1][i][0];
        for (j=1;j<mNetM+1;j++)
            sum += mOutput[0][j]*mWeight[1][i][j];
        mOutput[1][i]=1.0/(1.0+exp(-mNetGain*sum));
    }
    return &mOutput[1][1];
}// end of Propagate(float*)

float BpNn::Test(void){
    if (mTrainSize<=0)
        throw SizeError(); // forgot to input the training set etc.
    float TrainError=0;
    for (int i=0;i<mTrainSize;i++){
        Propagate(mInputBuffer[i]);
        TrainError+=ComputeMseError(mTargetBuffer[i]);
    }
    TrainError/=mTrainSize;
    std::cout<<"Train Err : "<<TrainError<<'\n';
    return TrainError;
}// end of Test(void)

float BpNn::Test(const char* testFileName,const int   testSize){
    std::ifstream TestFile(testFileName);
    if (!TestFile.is_open()){
        std::cerr<<"Failed to get the test set from "<<testFileName<<'\n';
        throw SizeError();
    }

    float TestError=0;
    float Input[mNetI],Target[mNetO];
    
    for (int i=0;i<testSize;i++){
        if (TestFile.eof())
            throw SizeError();
        
        int j;
        for (j=0;j<mNetI;j++)
            TestFile>>Input[j];
        char ts[255]; 
        TestFile.getline(ts,255);
        for (j=0;j<mNetO;j++)
            TestFile>>Target[j];
        TestFile.getline(ts,255);
        Propagate(Input);
        TestError+=ComputeMseError(Target);
    }
    TestFile.close();
    TestError/=testSize;
    std::cout<<"Test Err : "<<TestError<<'\n';
    return TestError;
}// end of Test(char*,int)

float BpNn::Train(void){// train with k fold cross_v
    if (mCrossFold>mTrainSize || mCrossFold<1)
        throw SizeError(); // forget to input the training set etc.
    
    float allCrossError=0;
    for (int CrossBeg=0;CrossBeg<mCrossFold;CrossBeg++){
        
        int beg=0,end =CrossBeg*mTrainSize/mCrossFold;
        for (int i=beg;i<end;i++){
            Propagate(mInputBuffer[i]);
            ComputeMseError(mTargetBuffer[i]);
            BackPropagate();
            AdjustWeights(mInputBuffer[i]);
        }
 
        beg=(CrossBeg+1)*mTrainSize/mCrossFold; end=mTrainSize;
        for (int i=beg;i<end;i++){
            Propagate(mInputBuffer[i]);
            ComputeMseError(mTargetBuffer[i]);
            BackPropagate();
            AdjustWeights(mInputBuffer[i]);
        }
        
        beg=CrossBeg*mTrainSize/mCrossFold;
        end=(CrossBeg+1)*mTrainSize/mCrossFold;
        for (int i=beg;i<end;i++){
            Propagate(mInputBuffer[i]);
            allCrossError+=ComputeMseError(mTargetBuffer[i]);
        }
    }
    allCrossError/=(mTrainSize/mCrossFold)*mCrossFold;
    
    return allCrossError;
}// end of Train(void)

////////////////////////////////////////////////////////////////////////////////
// 2. with the soft max transition function

const float* BpNn::SoftMaxPro(const float* inputData){
    register int i,j; register float sum;
    
    for (i=1;i<mNetM+1;i++){
        sum = mNetBias*mWeight[0][i][0];
        for (j=1;j<mNetI+1;j++)
            sum += inputData[j-1]*mWeight[0][i][j];
        mOutput[0][i]=1.0/(1.0+exp(-mNetGain*sum));
        // the sigmoid function for the first layer
    }
    
    sum=0.;
    for (i=1;i<mNetO+1;i++){
        mOutput[1][i]=mNetBias*mWeight[1][i][0];
        for (j=1;j<mNetM+1;j++)
            mOutput[1][i] += mOutput[0][j]*mWeight[1][i][j];
        mOutput[1][i]=exp(mOutput[1][i]);
        sum+=mOutput[1][i];
    }
    sum=1./sum;
    for (i=1;i<mNetO+1;i++) 
        mOutput[1][i]=mOutput[1][i]*sum;
    return &mOutput[1][1];  
}// end of SoftMaxPro(float*)

float BpNn::SoftMaxTrain(void){ // k fold cross validation training
    if (mCrossFold>mTrainSize || mCrossFold<1)
        throw SizeError(); // forget to input the training set etc.
  
    float allCrossError=0;
    for (int CrossBeg=0;CrossBeg<mCrossFold;CrossBeg++){
        
        int beg=0,end =CrossBeg*mTrainSize/mCrossFold;
        for (int i=beg;i<end;i++){
            SoftMaxPro(mInputBuffer[i]);
            ComputeCeError(mTargetBuffer[i]);
            BackPropagate();
            AdjustWeights(mInputBuffer[i]);
        }
 
        beg=(CrossBeg+1)*mTrainSize/mCrossFold; 
        end=(mTrainSize/mCrossFold)*mCrossFold; // some last samples will not 
        for (int i=beg;i<end;i++){              // be used
            SoftMaxPro(mInputBuffer[i]);
            ComputeCeError(mTargetBuffer[i]);
            BackPropagate();
            AdjustWeights(mInputBuffer[i]);
        }
        
        beg=CrossBeg*mTrainSize/mCrossFold;
        end=(CrossBeg+1)*mTrainSize/mCrossFold;
        for (int i=beg;i<end;i++){
            SoftMaxPro(mInputBuffer[i]);
            allCrossError+=ComputeCeError(mTargetBuffer[i]);
        }
    }
    allCrossError/=(mTrainSize/mCrossFold)*mCrossFold;
    return allCrossError;
}// end of SoftMaxTrain(void)

float BpNn::SoftMaxTest(void){
    float CeErr=0.;
    float QuadErr=0.;
    for (int i=0;i<mTrainSize;i++){
        SoftMaxPro(mInputBuffer[i]);
        CeErr  +=ComputeCeError (mTargetBuffer[i]);
	QuadErr+=ComputeMseError(mTargetBuffer[i]);
    }
    CeErr  /=mTrainSize;
    QuadErr/=mTrainSize;
    std::cout<<"Cross Entropy Train Err: "<<CeErr<<' ';
    std::cout<<"Quad Train Err: "<<QuadErr<<'\n';
    return CeErr;
}// end of SoftMaxTest

float BpNn::SoftMaxTest(const char*testFileName,const int  testSize){
    std::ifstream TestFile(testFileName);
    if (!TestFile.is_open()){
        std::cout<<"Open test file "<<testFileName<<" Error\n";
        throw SizeError();
    }
    
    float TestError=0,QuadError=0;
    float Input[mNetI], Target[mNetO];
    
    for (int i=0;i<testSize;i++){
        if (TestFile.eof())
            throw SizeError();
        for (int j=0;j<mNetI;j++){
            TestFile>>Input[j];
        }
        char ts[255]; 
        TestFile.getline(ts,255);
        
        for (int j=0;j<mNetO;j++){
            TestFile>>Target[j];
        }
        TestFile.getline(ts,255);
        
        SoftMaxPro(Input);
        TestError+=ComputeCeError(Target);
        QuadError+=ComputeMseError(Target);
    }
    TestFile.close();
    
    TestError/=testSize;
    QuadError/=testSize;
    std::cout<<"Cross Entropy Test Err: "<<TestError
             <<"  Squired difference Test Err: "<<QuadError<<'\n';
    return TestError;
}// end of SoftMaxTest(char*,int)

// //////////////////////////////////////////////////////////////////////
// // private members                                                  //
// //////////////////////////////////////////////////////////////////////
// void BpNn::init(int netI,int netM,int netO);// allocat memory
// void BpNn::destroy(void);  //release memory, to be used in the destructor
// 

void BpNn::init(int netI,int netM,int netO,bool defaultParameters){
    if (netI<0||netM<0||netO<0)  throw SizeError();
    
    if (defaultParameters){
        mNetAlpha =0.5;   mNetEta   =0.001;
        mNetBias  =1;     mNetGain  =1;
        mTrainSize=0;     mCrossFold=10;
    }

    mNetI=netI;mNetM=netM;mNetO=netO; 
    // begin to alocate memory
    
    mOutputs = new float [mNetM+1 + mNetO+1];
    mOutput[0]=&mOutputs[0];
    mOutput[1]=&mOutputs[mNetM+1];
    mOutput[0][0] = mNetBias; // the first unit was used for the threshold
    mOutput[1][0] = mNetBias; // for the sigmoid function
    
    mErrors = new float [mNetM+1 + mNetO+1];
    mError[0]=&mErrors[0];
    mError[1]=&mErrors[mNetM+1];
    mError [0][0] = 0;
    mError [1][0] = 0;

    mWeights =new float [(mNetM+1)*(mNetI+1) + (mNetM+1)*(mNetO+1)];
    mdWeights=new float [(mNetM+1)*(mNetI+1) + (mNetM+1)*(mNetO+1)];
    
    mWeight [0]   = new float * [mNetM+1];
    mdWeight[0]   = new float * [mNetM+1];
    for (int i=0;i<(mNetM+1);i++){
        mWeight [0][i] =&mWeights [(mNetI+1)*i];
        mdWeight[0][i] =&mdWeights[(mNetI+1)*i];
    }

    mWeight [1]   = new float * [mNetO+1];
    mdWeight[1]   = new float * [mNetO+1];
    for (int i=0;i<(mNetO+1);i++){
        mWeight [1][i] = &mWeights [(mNetM+1)*(mNetI+1)+(mNetM+1)*i];
        mdWeight[1][i] = &mdWeights[(mNetM+1)*(mNetI+1)+(mNetM+1)*i];
    }
    
    for (int j=0;j<mNetI+1;j++)
        mWeight[0][0][j]=0;
    for (int j=0;j<mNetM+1;j++)
        mWeight[1][0][j]=0;
    
    for (int i=0;i<mNetM+1;i++)
        for (int j=0;j<mNetI+1;j++)
            mdWeight[0][i][j]=0;
    for (int i=0;i<mNetO+1;i++)
        for (int j=0;j<mNetM+1;j++)
            mdWeight[1][i][j]=0;
}// end of init(int,int,int)

void BpNn::destroy(void){
    delete[] mOutputs;
    delete[] mErrors;
    
    delete[]mWeight [0];delete[]mWeight [1];
    delete[]mdWeight[0];delete[]mdWeight[1];
    
    delete[]mWeights ;  delete[]mdWeights;
    return;    
}// end of destroy()

float BpNn::ComputeMseError(const float* target){
    // only after PropagateNn(float*)
    register float  Out,Err,TotalError=0;
    for (int i=1;i<mNetO+1;i++){
        Out = mOutput[1][i];
        Err = target[i-1]-Out;
        TotalError +=Err*Err;
        mError[1][i]=mNetGain*Out*(1.0-Out)*Err;
    }
    return TotalError;
}// end of ComputeMseError(float*)

float BpNn::ComputeCeError(float const*target){
    register float TotalError=0.;
    for (int i=1;i<mNetO+1;i++){
	mError[1][i]=target[i-1]-mOutput[1][i];
        TotalError -=target[i-1]*log(mOutput[1][i]);
    }
    return TotalError;
}// end of ComputeCeError(const float*target)

void BpNn::BackPropagate(void){
    // only after Compute*Error(float*)
    register float Err;
    for (int i=1;i<mNetM+1;i++){
        Err = 0.;
        for (int j=1;j<mNetO+1;j++)
            Err += mWeight[1][j][i]*mError[1][j];
        mError[0][i]=mNetGain*mOutput[0][i]*(1.0-mOutput[0][i])*Err;
    }     // of course Error[0][1]===0;
}// end of BackPropagate()
 
void BpNn::AdjustWeights(const float* inputData){
    // only after BackpropagateNn()
    register int i,j;
    register float DWeight;
    
    for (i=1;i<mNetM+1;i++){
        DWeight=mdWeight[0][i][0]; // when j==0;
        mWeight [0][i][0]+= mNetEta*mError[0][i]*mNetBias
                           +mNetAlpha*DWeight;
        mdWeight[0][i][0] = mNetEta*mError[0][i]*mNetBias;
                
        for (j=1;j<mNetI+1;j++){
            DWeight=mdWeight[0][i][j];
            mWeight [0][i][j]+= mNetEta*mError[0][i]*inputData[j-1]
                               +mNetAlpha*DWeight;
            mdWeight[0][i][j] = mNetEta*mError[0][i]*inputData[j-1];
        }
    }
    for (i=1;i<mNetO+1;i++)
        for (j=0;j<mNetM+1;j++){
            DWeight=mdWeight[1][i][j];
            mWeight [1][i][j]+= mNetEta*mError[1][i]*mOutput[0][j]
                               +mNetAlpha*DWeight;
            mdWeight[1][i][j] = mNetEta*mError[1][i]*mOutput[0][j];
        }
    return;
}// end of AdjustWeights()

////////////////////////////////////////////////////////////////////////////////
//  for batch training                                                        //
////////////////////////////////////////////////////////////////////////////////

void BpNn::ZeroDWeights(){
    register int i,j;
    for (i=1;i<mNetM+1;i++){
        mdWeight[0][i][0] = 0;
        for (j=1;j<mNetI+1;j++) mdWeight[0][i][j] = 0;
    }
    for (i=1;i<mNetO+1;i++)
        for (j=0;j<mNetM+1;j++){
            mdWeight[1][i][j] = 0;
        }
}// end of ZeroDWeights()

void BpNn::AddDWeights(const float* inputData){// only after BackpropagateNn()
    register int i,j;
    for (i=1;i<mNetM+1;i++){
        mdWeight[0][i][0] += mNetEta*mError[0][i]*mNetBias;
        for (j=1;j<mNetI+1;j++)
            mdWeight[0][i][j] += mNetEta*mError[0][i]*inputData[j-1];
    }
    for (i=1;i<mNetO+1;i++)
        for (j=0;j<mNetM+1;j++)
            mdWeight[1][i][j] += mNetEta*mError[1][i]*mOutput[0][j];
    return;
}// end of AddDWeights()

void BpNn::IncludeDWeights(){
    register int i,j;
    for (i=1;i<mNetM+1;i++){
        mWeight [0][i][0]+= mNetEta*mdWeight[0][i][0];
        for (j=1;j<mNetI+1;j++)
            mWeight [0][i][j]+= mNetEta*mdWeight[0][i][j];
    }
    for (i=1;i<mNetO+1;i++)
        for (j=0;j<mNetM+1;j++)
            mWeight [1][i][j]+= mNetEta*mdWeight[1][i][j];
}// end of IncludeDWeights()
        
float BpNn::BatchSoftMaxTrain(){ // k-fold cross-validation training
    if (mCrossFold>mTrainSize || mCrossFold<1)
        throw SizeError(); // forget to input the training set etc.
    
    float allCrossError=0;
    for (int CrossBeg=0;CrossBeg<mCrossFold;CrossBeg++){
        
        ZeroDWeights();
        
        int beg=0,end =CrossBeg*mTrainSize/mCrossFold;
        for (int i=beg;i<end;i++){
            SoftMaxPro(mInputBuffer[i]);
            ComputeCeError(mTargetBuffer[i]);
            BackPropagate();
            AddDWeights(mInputBuffer[i]);
        }
 
        beg=(CrossBeg+1)*mTrainSize/mCrossFold; 
        end=(mTrainSize/mCrossFold)*mCrossFold; // some last samples will not 
        for (int i=beg;i<end;i++){              // be used
            SoftMaxPro(mInputBuffer[i]);
            ComputeCeError(mTargetBuffer[i]);
            BackPropagate();
            AddDWeights(mInputBuffer[i]);
        }
        
        IncludeDWeights();
        
        beg=CrossBeg*mTrainSize/mCrossFold;
        end=(CrossBeg+1)*mTrainSize/mCrossFold;
        for (int i=beg;i<end;i++){
            SoftMaxPro(mInputBuffer[i]);
            allCrossError+=ComputeCeError(mTargetBuffer[i]);
        }
    }
    allCrossError/=(mTrainSize/mCrossFold)*mCrossFold;
    return allCrossError;
}// end of BatchSoftMaxTrain()

////////////////////////////////////////////////////////////////////////////////
//   end of BpNn                                                              //
////////////////////////////////////////////////////////////////////////////////
