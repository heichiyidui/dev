////////////////////////////////////////////////////////////////////////////////
// Copyright (C) 2002 Kuang Lin  kxlin@nimr.mrc.ac.uk                         //
////////////////////////////////////////////////////////////////////////////////

////////////////////////////////////////////////////////////////////////////////
// This is the implement file of class BpNn, a neural network using the       //
// standard sigmoid function or the soft max transition function, supervised- //
// trainned using the backpropagation algorithm.                              //
////////////////////////////////////////////////////////////////////////////////

#ifndef LK_BPNN_HPP
#define LK_BPNN_HPP

#include <cmath>
#include <cstdlib>
#include <cstring>
#include <fstream>
#include <iostream>

////////////////////////////////////////////////////////////////////////////////
//  the begin of class BpNn                                                   //
////////////////////////////////////////////////////////////////////////////////
class BpNn{
public:
////////////////////////////////////////////////////////////////////////////////
// constructors and the destructor                                            //
////////////////////////////////////////////////////////////////////////////////
    BpNn(){init(0,0,0,true);}
    BpNn(int netI,int netM,int netO){init(netI,netM,netO,true);}
    BpNn(char const*nnFileName){init(0,0,0,true);InputNn(nnFileName);}
    
    BpNn(const BpNn&); // copy initialize, only info that will be stored in file
        // will be copied
    BpNn& operator = (const BpNn&); // copy assigner, like the initialzer, only 
        // info that will be stored in file got copied.
    
    ~BpNn(); 
////////////////////////////////////////////////////////////////////////////////
//  Input parameters and neural network                                       //
////////////////////////////////////////////////////////////////////////////////
    
    BpNn& NetAlpha (const float netAlpha);// set the momentum factor
    BpNn& NetEta   (const float netEta  );// set the learninng rate
    BpNn& NetBias  (const float netBias );// set BIAS (threshold)
    BpNn& NetGain  (const float netGain );// set GAIN of Sigmoid function
    BpNn& CrossFold(const int crossFold );// set the number of fold 
        // for cross validation training
    
    void RandomWeights(const float weightRange);
        // randomly init weights to the scale [-weightRange, +weightRange]
    
    void InputNn (const char* nnFileName); // read Nn from a file
    void InputWeights(const float * inputWeights);// input weight matrices
    void InputTrainSet(char const* trainFileName,int trainSize);
        // read the training set into memory
    
////////////////////////////////////////////////////////////////////////////////
//  Output                                                                    //
////////////////////////////////////////////////////////////////////////////////
    class SizeError{};// exception class, report error in neural network size
    class HugeNumber{};// report error in propagation, numbers is too big/small
    
    void OutputNn(const char* outputFileName)const;// output Nn to a file
    void OutputNn(void)const;// put Nn to the standard output (no weights)
    
    int NetI() const {return mNetI;}// size of the input  layer
    int NetM() const {return mNetM;}// size of the hidden layer
    int NetO() const {return mNetO;}// size of the output layer
    
    float NetAlpha()const{return mNetAlpha;} // the momentum factor
    float NetEta  ()const{return mNetEta;}   // the learning rate
    
    float NetBias ()const{return mNetBias;}  // the threshold and 
    float NetGain ()const{return mNetGain;}  // the gain of the sigmoid function
    
    void OutputWeights(float* outputWeights)const;// copy weigths to array
    int  NetSize(void)const; // how many connection weights in the net
    
////////////////////////////////////////////////////////////////////////////////
//  Propagate, train and test of Net                                          //
////////////////////////////////////////////////////////////////////////////////
    // with the sigmoid transition function
    float const* Propagate(const float* inputData);// propagate the net
        // return a pointer to the output values
    float Train(void);// k fold cross_v training, sigmoid function
        // with the training set in memory
    float Test(void); // return the training error 
        // (after cross validatioin training)
    float Test(const char*testFileName,const int  testSize);
        // test on a set from a file
    
    // with the soft max trainsition function
    const float* SoftMaxPro(const float* input);// Propagate the net
        // return a pointer to the output values
    float SoftMaxTrain(void); // k fold cross validation training, set in mem
    float SoftMaxTest (void); // training error
    float SoftMaxTest (const char*testFileName,const int  testSize);  
        // test on a set from a file, cross entropy error 
    
    // batch training
    float BatchSoftMaxTrain(); // k-fold cross-validation, batch training
////////////////////////////////////////////////////////////////////////////////
//  private members                                                           //
////////////////////////////////////////////////////////////////////////////////
private:
    void init(int netI,int netM,int netO,bool defaultParameters=false);
        // allocat memory, set defaultParameters
    void destroy(void);  //release memory
    
    float ComputeMseError(const float* target); // Mean-Squired-Error
    float ComputeCeError(const float*target);   // cross entropy error
    
    void BackPropagate(void); // after Compute*Error
    void AdjustWeights(const float* inputData);//after BackPropagate()
    
    // for batch training
    void ZeroDWeights(); // no change of weights
    void AddDWeights(const float* inputData);  //add the value for this case
    void IncludeDWeights(); // include changes in this batch
    
    int mNetI,mNetM,mNetO;   // size of each layer

    float    mNetAlpha ;     // momentum factor of training
    float    mNetEta   ;     // learning rate   of training
    float    mNetBias  ;     // Bias (threshold) is often 1
    float    mNetGain  ;     // Gain of S function, often 1
    
    float*   mOutputs;       // all output values
    float*   mOutput[2];     // Output values of each layer
    
    float*   mErrors;        // all error terms
    float*   mError[2];      // error terms of each layer
    
    float*   mWeights;       // all conection weights
    float**  mWeight[2];     // pointers to conection weights of each layer
    float*   mdWeights;      // all last weight change
    float**  mdWeight[2];    // pointers to last weight change of each layer
    

    int      mTrainSize ;    // the size of training set
    float*   mTrainBuffer;   // values of training input and target
    float**  mInputBuffer;   // pointers to the input  value arrays
    float**  mTargetBuffer;  // pointers to the target value arrays
    
    int      mCrossFold ;    // often 10 fold cross-validation
};
////////////////////////////////////////////////////////////////////////////////
//   end of BpNn                                                              //
////////////////////////////////////////////////////////////////////////////////

#endif  // LK_BPNN_HPP

