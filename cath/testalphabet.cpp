////////////////////////////////////////////////////////////////////////////////
// the test program for the Alphabet class and the protein alphabet routine.  //
////////////////////////////////////////////////////////////////////////////////

#include "alphabet.hpp"

int main(int argc, char** argv)
{
    // 1. get an alphabet
    // 2. get a protein alphabet
    // 3. read a sequence file. 
    // 4. tranlate all chars and ints correctly
    
    using namespace std;
    
    cao::Alphabet abc; // IUPAC Protein Ambiguous 
    string name1="ARNDCQEGHILKMFPSTWYVBZX-a0783z_?x0C";
    
    cout<<"     string:           "<<name1<<'\n';
    cout<<"     isIncludedLetter: ";
    for (size_t i=0;i<name1.length();i++)
        cout<<abc.isIncluded(name1[i]);
    cout<<'\n';
    cout<<"     isResolvedLetter: ";
    for (size_t i=0;i<name1.length();i++)
        cout<<abc.isResolved(name1[i]);
    cout<<'\n';
    cout<<"     isGapLetter:      ";
    for (size_t i=0;i<name1.length();i++)
        cout<<abc.isGap(name1[i]);
    cout<<'\n';
    cout<<"     isUnknownLetter:  ";
    for (size_t i=0;i<name1.length();i++)
        cout<<abc.isUnknown(name1[i]);
    cout<<'\n';
    
    cout<<"     letter(state):    ";
    for (size_t i=0;i<name1.length();i++)
        if (abc.isIncluded(name1[i]))
            cout<<abc.letter(abc.state(name1[i]));
    cout<<'\n';
    
    cout<<"     abbr(state):      ";
    for (size_t i=0;i<name1.length();i++)
        if (abc.isIncluded(name1[i]))
            cout<<abc.abbr(abc.state(name1[i]))<<' ';
    cout<<'\n';
    
    cout<<"     letter(state(abbr(state))):\n"
        <<"                       ";
    for (size_t i=0;i<name1.length();i++)
        if (abc.isIncluded(name1[i]))
            cout<<abc.letter(abc.state(abc.abbr(abc.state(name1[i]))));
    cout<<'\n';
    
    
    
    return 0;
}
