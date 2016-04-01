#include <cmath>
#include <vector>
#include "cstdlib"
#include "snn.hpp"

// test the simple neural network

////////////////////////////////////////////////////////////////////////////////
// functions for reading domain tables                                        //
////////////////////////////////////////////////////////////////////////////////

size_t readTableLength(std::string domName){
    std::ifstream ifile;
    ifile.open(("/home/kuang/dev/yaspin/nntable/"+domName).c_str());
    assert(ifile); 
    size_t ti; ifile>>ti;
    ifile.close();
    return ti;
}// get domain table length

int readTable(std::string domName,size_t wingLength,
              float*input,float*output,bool*isAAInDom,size_t*state){
    using namespace std;
    ifstream ifile;
    ifile.open(("/home/kuang/dev/yaspin/nntable/"+domName).c_str());
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
    
    const size_t WING_LENGTH=9;
    // the number of neighbour residues on one side
    
    ////////////////////////////////////////////////////////////////////////////
    // start reading the test set                                             //
    ////////////////////////////////////////////////////////////////////////////
//    cout<<"reading the test set"<<'\n';
    ifstream index;
    
    index.open("/home/kuang/dev/yaspin/index/testtr.ls");
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
    // finished reading the test set, start testing                           //
    ////////////////////////////////////////////////////////////////////////////
    
    const float LOG_STATE_DIS[21]={
        -1.16053, -3.78586, -1.85348, -3.66232, -3.74353, -2.58915, -3.85147, 
        -2.68007, -4.64713, -2.53524, -4.87285, -4.65618, -3.1996 , -4.51839, 
        -3.64623, -5.17406, -2.74341, -5.383  , -5.5012 , -3.66944, -5.27308};

    Snn net; net.readNet("9_4_30_4.Nn");
    
    for (size_t l=0;l<ind.size();l++){// for each domain
        for (size_t i=0;i<length[l];i++){// for each amino acid
            if (!isInDom[l][i])
                continue;
            const float*tp=net.propagate(input[l]+i*21);
            for (size_t j=0;j<21;j++){
                cout<<log(tp[j])-LOG_STATE_DIS[j]<<' ';
            }
            cout<<'\n';
        }
    }
    
    ////////////////////////////////////////////////////////////////////////////
    // finished testing, start release memory                                 //
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
    
    return 0;
}
