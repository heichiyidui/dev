#include <string>
#include <vector>
#include <cmath>
#include <algorithm>
#include <iterator>
#include <iostream>
#include <fstream>
#include <numeric>

int main(int argc,char**argv){
    using namespace std;
    ifstream ifile(argv[1]);
    float score;
    int type;
    vector <float> cont_scores;
    vector <float> non_cont_scores;
    
    for (;ifile>>score>>type;){
        if (type==1)
            cont_scores.push_back(score);
        else
            non_cont_scores.push_back(score);
    }
    ifile.close();
    
    double mean0 = accumulate(non_cont_scores.begin(),   \
                              non_cont_scores.end(),0.0) \
                            / non_cont_scores.size();
    double mean1 = accumulate(cont_scores.begin(),   \
                              cont_scores.end(),0.0) \
                            / cont_scores.size();
    
    // calculate standard deviations
    vector<double> diff0(non_cont_scores.size());
    
    transform(non_cont_scores.begin(), non_cont_scores.end(), \
              diff0.begin(), bind2nd(minus<double>(), mean0));
                
    double sq_sum0 = inner_product(diff0.begin(), diff0.end(), \
                                   diff0.begin(), 0.0);
                
    double stdev0 = sqrt(sq_sum0 / non_cont_scores.size());

    vector<double> diff1(cont_scores.size());
    
    transform(cont_scores.begin(), cont_scores.end(), \
              diff1.begin(), bind2nd(minus<double>(), mean1));
                
    double sq_sum1 = inner_product(diff1.begin(), diff1.end(), \
                                   diff1.begin(), 0.0);
                
    double stdev1 = sqrt(sq_sum1 / cont_scores.size());
    
    cout<<mean0<<' '<<stdev0<<' '<<mean1<<' '<<stdev1<<'\n';
    return (0);
}
