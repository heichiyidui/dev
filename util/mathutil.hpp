///////////////////////////////////////////////////////////////////////////////
// Copyright (C) 2007-10  Kuang Lin kuanggong@gmail.com                      //
///////////////////////////////////////////////////////////////////////////////

#ifndef __mathutil_hpp__
#define __mathutil_hpp__

#include <cstdlib>
#include <cmath>
#include <vector>
#include <algorithm>
#include <numeric>


///////////////////////////////////////////////////////////////////////////////
// factors                                                                   //
///////////////////////////////////////////////////////////////////////////////

inline double gammaln(double xx){ // log of gamma function
    // Stirling's approximation: ln(n!) = n*ln(n)-n
    // or n! = sqrt(2*pi*n)*pow(n/e,n);
    // Lanczos approximation from Numerical Recipes
    // xx > 0 !
    double x,y, tmp, ser; 
    const static double cof[6]={76.18009172947146,-86.50532032941677,
        24.01409824083091,-1.231739572450155,
        0.1208650973866179e-2,-0.5395239384953e-5};
    y=x=xx;     tmp=x+5.5;
    tmp-= (x+0.5)*log(tmp);
    ser=1.000000000190015;
    for (int j=0;j<6;j++) ser+=cof[j]/++y;
    return -tmp+log(2.5066282746310005*ser/x);
}// from Numerical Recipes


inline double factor(size_t n){ // from Numerical Recipes
    static size_t ntop=4;
    static double a[33]={1.0, 1.0, 2.0, 6.0, 24.0};
    if (n>32) return exp(gammaln(n+1.0));
    while (ntop<n){
        int j=ntop++;
        a[ntop]=a[j]*ntop;
    }
    return a[n];
}// again, from Numerical Recipes

inline double lnbico(int N,int K) // log binomial coefficient
    {return log(factor(N))-log(factor(K))-log(factor(N-K));}

////////////////////////////////////////////////////////////////////////////////
// Kullback-Leibler distances                                                 //
////////////////////////////////////////////////////////////////////////////////

inline float computeReDis(float const* p, float const* q, const size_t n){
    float re=0;
    for (size_t i=0;i<n;i++)
        if (q[i]>0 && p[i]>0)
            re+=p[i]*(log(p[i]/q[i]));
    return re;  
}// calculate relative entropy distance between two distributions

inline float computeKlDis(float const* p, float const* q, const size_t n,
                          const float lambda=0.5){
    float m,dis=0;
    for (size_t i=0;i<n;i++){
        m=p[i]*lambda+q[i]*(1-lambda);
        if ( m > 0 && p[i] > 0 && q[i] > 0 ){
            dis+=  lambda    * p[i]*log(p[i]/m) ;
            dis+= (1-lambda) * q[i]*log(q[i]/m);
        }
    }
    return dis;
}// calculate (symmetric) Kullback-Leibler divergence between two distributions

////////////////////////////////////////////////////////////////////////////////
// mean, standard deviation etc...                                            //
////////////////////////////////////////////////////////////////////////////////

inline float computeMean(const float* tArray, const size_t nT)
    {return std::accumulate(tArray,tArray+nT,0.0F)/nT;}

inline float computeVar(const float* tArray, const size_t nT){
    float mean = computeMean(tArray,nT); 
    float var=0.0F; 
    for (size_t i=0;i<nT;i++)
        var += (tArray[i] - mean) * (tArray[i] - mean);
    return var/= nT-1.0F ;
}   // estimate the variance

inline float computeStdDev(const float* tArray, const size_t nT)
    {return sqrt(computeVar(tArray,nT));} 
    // estimate the standard deviation

inline float computeCorrCoef(const float* xArray, 
                             const float* yArray, const size_t n){
    float meanX=computeMean(xArray,n), meanY=computeMean(yArray,n);
    float SSxx=0.0F, SSyy=0.0F, SSxy=0.0F;
    for (size_t i=0;i<n;i++){
        SSxx += xArray[i] * xArray[i];
        SSyy += yArray[i] * yArray[i];
        SSxy += xArray[i] * yArray[i];
    }
    return  (n*SSxy - n*meanX*n*meanY) /
       sqrt((n*SSxx - n*meanX*n*meanX) * (n*SSyy - n*meanY*n*meanY));
}// end of estimating Pearson's correlation coefficient
 // could output 'nan' because there is no checking of 
 // (n*SSxx - n*meanX*n*meanX) or (n*SSyy - n*meanY*n*meanY) being zero

inline float computeMse(const float* xArray, 
                        const float* yArray, const size_t n){
    float tsum(0.0F); 
    for (size_t i=0;i<n;i++) 
        tsum+=(xArray[i]-yArray[i])*(xArray[i]-yArray[i]);
    return tsum/n;
}// calculate MSE, mean squared error

inline float computeRmse(const float* xArray,
                         const float* yArray, const size_t n){
    return sqrt(computeMse(xArray,yArray,n));
}// calculate RMSE, root mean squared error

inline float rand01(){return float(rand())/RAND_MAX;} 
// get a random sample of the uniform distribution [0,1]

inline float randNorm(){
    float val=0;for(size_t i=0;i<100;i++)val+=rand01();return (val-50)/2.88;
}// get a random sample of the normal distribution N(0,1)

///////////////////////////////////////////////////////////////////////////////

inline float norm_isf(const float p){
    float const static c[] = {2.515517, 0.802853, 0.010328};
    float const static d[] = {1.432788, 0.189269, 0.001308};
    if (p<0.5){
        float t=sqrt(-2.0*log(p));
        return t - ((c[2]*t + c[1])*t + c[0]) / 
               (((d[2]*t + d[1])*t + d[0])*t + 1.0);
    } else{
        float t=sqrt(-2.0*log(1-p));
    
        return -t + ((c[2]*t + c[1])*t + c[0]) / 
               (((d[2]*t + d[1])*t + d[0])*t + 1.0);
    }
}// rank to normal scores, from Abramowitz and Stegun formula 26.2.23.
// http://www.johndcook.com/normal_cdf_inverse.html

#endif // __mathutil_hpp__ 
