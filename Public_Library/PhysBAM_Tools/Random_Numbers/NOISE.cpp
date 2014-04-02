//#####################################################################
// Copyright 2004-2007, Ron Fedkiw, Eran Guendelman, Andrew Selle.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#include <PhysBAM_Tools/Random_Numbers/NOISE.h>
#include <cmath>
#include <cstdlib>
using namespace PhysBAM;

static const int number_points=512;
template<class T> const int NOISE<T>::noise_seed=666;
template<class T> const int NOISE<T>::prime1=173;
template<class T> const int NOISE<T>::prime2=263;
template<class T> const int NOISE<T>::prime3=337;
template<class T> const T NOISE<T>::phi=(T).6180339;
template<class T> T NOISE<T>::points[number_points];
template<class T> bool NOISE<T>::initialized=false;

//#####################################################################
// Function Hash
//#####################################################################
template<class T> int NOISE<T>::
Hash(const int a, const int b,const int c)
{
    return (a+b+c) & (number_points-1);
}
//#####################################################################
// Function Initialize_Noise
//#####################################################################
template<class T> void NOISE<T>::
Initialize_Noise()
{  
    srand(noise_seed);
    for(int i=0;i<number_points;i++) points[i]=(T)rand()/RAND_MAX;
}
//#####################################################################
// Function Noise_Function1
//#####################################################################
template<class T> T NOISE<T>::
Noise_Function1(const VECTOR<T,3>& p)
{
    int xi,yi,zi,xa,xb,ya,yb,za,zb;
    T xf,yf,zf,p000,p001,p010,p011,p100,p101,p110,p111;
    xf=p.x;yf=p.y;zf=p.z;
    xi=(int)floor(xf);xa=(int)floor(prime1*(xi*phi-floor(xi*phi)));xb=(int)floor(prime1*((xi+1)*phi-floor((xi+1)*phi)));
    yi=(int)floor(yf);ya=(int)floor(prime2*(yi*phi-floor(yi*phi)));yb=(int)floor(prime2*((yi+1)*phi-floor((yi+1)*phi)));
    zi=(int)floor(zf);za=(int)floor(prime3*(zi*phi-floor(zi*phi)));zb=(int)floor(prime3*((zi+1)*phi-floor((zi+1)*phi)));
    p000=points[Hash(xa,ya,za)];p001=points[Hash(xa,ya,zb)];p010=points[Hash(xa,yb,za)];p011=points[Hash(xa,yb,zb)];
    p100=points[Hash(xb,ya,za)];p101=points[Hash(xb,ya,zb)];p110=points[Hash(xb,yb,za)];p111=points[Hash(xb,yb,zb)];
    xf-=xi;yf-=yi;zf-=zi;
    // bend the line to fake a spline interpolation
    if(xf>.5){xf=1-xf;xf=2*xf*xf;xf=1-xf;}else{xf=2*xf*xf;}
    if(yf>.5){yf=1-yf;yf=2*yf*yf;yf=1-yf;}else{yf=2*yf*yf;}
    if(zf>.5){zf=1-zf;zf=2*zf*zf;zf=1-zf;}else{zf=2*zf*zf;}
    return p000*(1-xf)*(1-yf)*(1-zf)+p001*(1-xf)*(1-yf)*zf+p010*(1-xf)*yf*(1-zf)+p011*(1-xf)*yf*zf+p100*xf*(1-yf)*(1-zf)+p101*xf*(1-yf)*zf+p110*xf*yf*(1-zf)+p111*xf*yf*zf;
}
//#####################################################################
// Function Noise_Function3
//#####################################################################
template<class T> VECTOR<T,3> NOISE<T>::
Noise_Function3(const VECTOR<T,3>& p)
{
    VECTOR<T,3> tmp,v;
    v.x=Noise_Function1(p);tmp.x=p.y+prime1;tmp.y=p.x+prime2;tmp.z=p.z+prime3;
    v.y=Noise_Function1(tmp);tmp.x+=prime2;tmp.y+=prime3;tmp.z+=prime1;
    v.z=Noise_Function1(tmp);
    return v;
}
//#####################################################################
// Function Noise1
//#####################################################################
template<class T> T NOISE<T>::
Noise1(const VECTOR<T,3>& p_input,const int octaves,const T persistance)
{
    if(!initialized) Initialize_Noise();T result=0,factor=1;VECTOR<T,3> p=p_input;
    for(int i=0;i<octaves;i++,factor*=persistance){result+=Noise_Function1(p)*factor;p*=2;}
    return result/Normalization_Factor(octaves,persistance);
}
//#####################################################################
// Function Noise3
//#####################################################################
template<class T> void NOISE<T>::
Noise3(const VECTOR<T,3>& p_input,VECTOR<T,3>& v,const int octaves,const T persistance)
{
    VECTOR<T,3> p=p_input;
    if (!initialized) Initialize_Noise();VECTOR<T,3> c;T factor=1;v.x=v.y=v.z=0;
    for(int i=0;i<octaves;i++,factor*=persistance){v=Noise_Function3(p);c*=factor;v+=c;p*=2;}
}
//#####################################################################
// Function Normalization_Factor
//#####################################################################
template<class T> T NOISE<T>::
Normalization_Factor(const int octaves,const T persistance)
{  
    return (1-persistance)/(1-pow(persistance,octaves));
}
//#####################################################################
// Function Fractional_Brownian
//#####################################################################
template<class T> T NOISE<T>::
Fractional_Brownian(const VECTOR<T,3>& p,const int octaves,const T lacunarity,const T gain)
{
    if(!initialized) Initialize_Noise();
    T sum=0,lambda=1,o=1;
    for(int i=0;i<octaves;i++){
        sum+=o*(2*Noise_Function1(lambda*p)-1);
        lambda*=(T)lacunarity;o*=gain;}
    return (T).5*(sum+1);
}
//#####################################################################
// Function Turbulence
//#####################################################################
template<class T> T NOISE<T>::
Turbulence(const VECTOR<T,3>& p,const int octaves,const T lacunarity,const T gain)
{
    if(!initialized) Initialize_Noise();
    T sum=0,lambda=1,o=1;
    for(int i=0;i<octaves;i++){
        sum+=o*abs(2*Noise_Function1(lambda*p)-1);
        lambda*=(T)lacunarity;o*=gain;}
    return (T).5*(sum+1);
}
//#####################################################################
// Function Fractional_Brownian
//#####################################################################
template<class T> VECTOR<T,3> NOISE<T>::
Fractional_Brownian_3D(const VECTOR<T,3>& p,const int octaves,const T lacunarity,const T gain)
{
    if(!initialized) Initialize_Noise();
    VECTOR<T,3> sum;T lambda=1,o=1;
    for(int i=0;i<octaves;i++){
        sum+=o*((T)2*Noise_Function3(lambda*p)-VECTOR<T,3>::All_Ones_Vector());
        lambda*=(T)lacunarity;o*=gain;}
    return (T).5*(sum+1);
}
//#####################################################################
// Function Turbulence
//#####################################################################
template<class T> VECTOR<T,3> NOISE<T>::
Turbulence_3D(const VECTOR<T,3>& p,const int octaves,const T lacunarity,const T gain)
{
    if(!initialized) Initialize_Noise();
    VECTOR<T,3> sum;T lambda=1,o=1;
    for(int i=0;i<octaves;i++){
        sum+=o*abs((T)2*Noise_Function3(lambda*p)-VECTOR<T,3>::All_Ones_Vector());
        lambda*=(T)lacunarity;o*=gain;}
    return (T).5*(sum+1);
}
//#####################################################################
template class NOISE<float>;
template class NOISE<double>;
