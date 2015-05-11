#include <cmath>
#include <vector>
#include "cstdlib"
#include "snn.hpp"

// train and test the simple neural network

// usage: trainNn 6 1e-4 30 1 6_4_30_1.Nn
//        trainNn WingLength LearningRate NumberHiddenUnits RandomSeed

////////////////////////////////////////////////////////////////////////////////
// functions for reading domain tables                                        //
////////////////////////////////////////////////////////////////////////////////

size_t readTableLength(std::string domName){
    std::ifstream ifile;
    ifile.open(("/home/klinbrc/scratch/dev/yaspin/nntable/"+domName).c_str());
    assert(ifile); 
    size_t ti; ifile>>ti;
    ifile.close();
    return ti;
}// get domain table length

int readTable(std::string domName,size_t wingLength,
              float*input,float*output,bool*isAAInDom,size_t*state){
    using namespace std;
    ifstream ifile;
    ifile.open(("/home/klinbrc/scratch/dev/yaspin/nntable/"+domName).c_str());
    assert(ifile); 
    size_t length;
    ifile>>length;
    
    for (size_t i=0;i<wingLength;i++){
        for (size_t j=0;j<21;j++)
            input[i*21+j]=-1.0;
        input[i*21+20]=1.0;
    }
    for (size_t i=length+wingLength;i<length+wingLength*2;i++){
        for (size_t j=0;j<21;j++)
            input[i*21+j]=-1.0;
        input[i*21+20]=1.0;
    }
    
    bool tBool; char SS; size_t contNum; float tf;
    
    for (size_t i=0;i<length;i++){
        ifile>>tBool>>SS>>contNum;
        isAAInDom[i]=tBool;
        
        // CGHIDEF 7 states to 21 states
        size_t ssState=0;
        switch (SS){
            case 'C':{ ssState=0  ;break;}
            case 'G':{ ssState=1  ;break;}
            case 'H':{ ssState=2  ;break;}
            case 'I':{ ssState=3  ;break;}
            case 'D':{ ssState=4  ;break;}
            case 'E':{ ssState=5  ;break;}
            case 'F':{ ssState=6  ;break;}
            case 'U':{            ;break;}// unknown for AA not in domain
            default :{
                cout<<"WHAT SS? "<<SS<<'\n';
                break;
            }
        }
        if (contNum>6) ssState+=7;
        if (contNum>8) ssState+=7;
        state[i]=ssState;
        
        for (size_t j=0;j<21;j++)
            output[i*21+j]=0.0;
        output[i*21+ssState]=1.0;
        
        for (size_t j=0;j<20;j++){
            ifile>>tf;
            input[(i+wingLength)*21+j]=tf/5.0;
        }
        input[(i+wingLength)*21+20]=-1.0;
    }
    ifile.close();
    return 0;
}

////////////////////////////////////////////////////////////////////////////////
// the main program                                                           //
////////////////////////////////////////////////////////////////////////////////
int main(int argc, char** argv)
{
    using namespace std;
    
    const size_t WING_LENGTH=atoi(argv[1]);
    // the number of neighbour residues on one side
    const float LEARNING_RATE=atof(argv[2]);
    // the learning rate
    const size_t N_MID=atoi(argv[3]);
    // the number of units in the hidden (middle) layer of the network    
    const size_t RSEED=atoi(argv[4]);
    
    srandom(time(0)+RSEED);
    random();random();random();random();random();
    
    ////////////////////////////////////////////////////////////////////////////
    // start reading the training set                                         //
    ////////////////////////////////////////////////////////////////////////////
    cout<<"reading the training set"<<'\n';
    ifstream index;
    
    index.open("/home/klinbrc/scratch/dev/yaspin/index/train.ls");
    assert(index); 
    string ts;
    vector <string> ind;
    for (;index>>ts;)
        ind.push_back(ts);
    index.close();
    
    random_shuffle(ind.begin(),ind.end());
    
    vector <float *> input;  // input to the network
    vector <float *> output; // output
    vector <bool  *> isInDom;// is the residue in the domain
    vector <size_t*> state;  // residue SS and burial state
    vector <size_t > length; // domain lengths
    
    for (size_t l=0;l<ind.size();l++){
        string domName=ind[l];
        
        size_t domLength=readTableLength(domName);
        
        float* domInput =new float[(domLength+WING_LENGTH*2)*21];
        float* domOutput=new float[domLength*21];
        bool*  isAAInDom=new bool[domLength];
        size_t*domState =new size_t[domLength];
        
        readTable(domName, WING_LENGTH, domInput,domOutput, isAAInDom,domState);
        
        input.push_back(domInput);
        output.push_back(domOutput);
        isInDom.push_back(isAAInDom);
        state.push_back(domState);
        length.push_back(domLength);
    }
    
    ////////////////////////////////////////////////////////////////////////////
    // start reading the test set                                             //
    ////////////////////////////////////////////////////////////////////////////
    cout<<"reading the test set"<<'\n';
    
    index.open("/home/klinbrc/scratch/dev/yaspin/index/testtr.ls");
    assert(index); 
    vector <string> testInd;
    for (;index>>ts;)
        testInd.push_back(ts);
    index.close();
        
    vector <float *> testInput;  // input to the network
    vector <float *> testOutput; // output
    vector <bool  *> testIsInDom;// is the residue in the domain
    vector <size_t*> testState;  // residue SS and burial state
    vector <size_t > testLength; // domain lengths
    
    for (size_t l=0;l<testInd.size();l++){
        string domName=testInd[l];
        size_t domLength=readTableLength(domName);
        
        float* domInput =new float[(domLength+WING_LENGTH*2)*21];
        float* domOutput=new float[domLength*21];
        bool*  isAAInDom=new bool[domLength];
        size_t*domState =new size_t[domLength];
        
        readTable(domName, WING_LENGTH, domInput,domOutput, isAAInDom,domState);
        
        testInput.push_back(domInput);
        testOutput.push_back(domOutput);
        testIsInDom.push_back(isAAInDom);
        testState.push_back(domState);
        testLength.push_back(domLength);
    }
    
    ////////////////////////////////////////////////////////////////////////////
    // finished reading the training and test set, start training             //
    ////////////////////////////////////////////////////////////////////////////
    
    size_t nin=21*(WING_LENGTH*2+1),nout=21;
    
    Snn net(nin,N_MID,nout);
    
    const size_t PRIOR_STATE_DUP[21]=
        {1,14,2,12,13,4,15,4,32,4,40,33,8,28,12,55,5,68,76,12,61};
    // the snow-ball burn in iteratioins
    for (size_t iter=20;iter>0;iter--){
        size_t priorStateDup[21];
        for (size_t i=0;i<21;i++)
            priorStateDup[i]=size_t(pow(PRIOR_STATE_DUP[i],iter/20.0));
        for (size_t l=0;l<ind.size();l++){// for each domain
            for (size_t i=0;i<length[l];i++){// for each amino acid
                if (!isInDom[l][i])
                    continue;
                for (size_t j=0;j<priorStateDup[state[l][i]];j++){
                    net.propagate(input[l]+i*21);
                    net.backpropagate(input[l]+i*21,output[l]+i*21);
                    // update the weights in training
                    net.updateWeights(0.5,LEARNING_RATE); 
                    net.nullDWeights();
                }
            }
        }
    }
    
    /// the real training iterations  
    double lastTestErr=1e+8;
    for (size_t iter=0;iter<500;iter++){ // for each iteration
        double trainErr=0;
        for (size_t l=0;l<ind.size();l++){// for each domain
            for (size_t i=0;i<length[l];i++){// for each amino acid
                if (!isInDom[l][i])
                    continue;
                net.propagate(input[l]+i*21);
                net.backpropagate(input[l]+i*21,output[l]+i*21);
                trainErr+=net.computeCrossEntropy(output[l]+i*21);
                // update the weights in training
                net.updateWeights(0.5,LEARNING_RATE); 
                net.nullDWeights();
            }
        }
        
        cout<<iter<<'\t';
        cout<<trainErr<<'\t';
        
        double testErr=0;
        for (size_t l=0;l<testInd.size();l++){
            for (size_t i=0;i<testLength[l];i++){
                if (!testIsInDom[l][i])
                    continue;
                net.propagate(testInput[l]+i*21);
                testErr+=net.computeCrossEntropy(testOutput[l]+i*21);
            }
        }
        cout<<testErr<<'\n';
        if (testErr > lastTestErr)
            break;
        lastTestErr=testErr;
    }
    ////////////////////////////////////////////////////////////////////////////
    // finished training, start release memory and save network               //
    ////////////////////////////////////////////////////////////////////////////
    
    for (size_t i=0;i<ind.size();i++){
        delete[] input[i];
        delete[] output[i];
        delete[] isInDom[i];
        delete[] state[i];
    }
    
    ind.clear();
    input.clear();  
    output.clear(); 
    isInDom.clear();
    state.clear();  
    length.clear(); 
    
    for (size_t i=0;i<testInd.size();i++){
        delete[] testInput[i];
        delete[] testOutput[i];
        delete[] testIsInDom[i];
        delete[] testState[i];
    }
    
    testInd.clear();
    testInput.clear();
    testOutput.clear();
    testIsInDom.clear();
    testState.clear();
    testLength.clear();
    
    net.writeNet(argv[5]); 
    return 0;
}
