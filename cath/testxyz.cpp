////////////////////////////////////////////////////////////////////////////////
// test the XYZ vector class                                                  //
////////////////////////////////////////////////////////////////////////////////
#include <iostream>
#include "xyz.hpp"

int main(int argc, char** argv)
{
    
    using namespace std;
    cao::XYZ v1(1,0,0);
    cao::XYZ v2(0,1,0.001);
    cao::XYZ v3=v1.cross(v2);
    cout<<v3.x<<' '<<v3.y<<' '<<v3.z<<'\n';
    cout<<v1.dot(v1)<<'\n';
    cout<<v1.dot(v2)<<'\n';
    
    cao::XYZ const vc(0,0,1);
    cout<<vc.z<<'\n';
//    v1 =v1*2.3F ;
//    cout<<-v1.x<<'\n';
//    cout<<abs(v1)<<'\n';
//    cout<<"OK"<<'\n';
    return 0;
}
