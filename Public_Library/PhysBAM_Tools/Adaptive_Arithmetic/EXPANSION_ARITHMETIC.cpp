//#####################################################################
// Copyright 2008, Jeffrey Hellrung, Eftychios Sifakis.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class EXPANSION_ARITHMETIC
//#####################################################################
#include <PhysBAM_Tools/Adaptive_Arithmetic/EXPANSION_ARITHMETIC.h>
#include <PhysBAM_Tools/Arrays/ARRAY.h>
#include <PhysBAM_Tools/Arrays_Computations/HEAPIFY.h>
#include <PhysBAM_Tools/Log/DEBUG_UTILITIES.h>
using namespace PhysBAM;
template<> const float EXPANSION_ARITHMETIC<float>::splitter=4097.f;
template<> const double EXPANSION_ARITHMETIC<double>::splitter=134217729.l;
#define INEXACT volatile
//#####################################################################
// Primitive operation helpers
//#####################################################################
namespace{
template<class T> inline void Fast_Two_Sum(const T a,const T b,INEXACT T& x,INEXACT T& y)
{x=(T)(a+b);INEXACT T bvirt=x-a;y=b-bvirt;}

template<class T> inline void Two_Sum(const T a,const T b,INEXACT T& x,INEXACT T& y)
{x=(T)(a+b);INEXACT T bvirt=(T)(x-a);T avirt=x-bvirt,bround=b-bvirt,around=a-avirt;y=around+bround;}
    
template<class T> inline void Split(const T a,T& ahi,T& alo)
{INEXACT T c=(T)(EXPANSION_ARITHMETIC<T>::splitter*a),abig=(T)(c-a);ahi=c-abig;alo=a-ahi;}

template<class T> inline void Two_Product_Presplit(const T a,const T b,const T bhi,const T blo,INEXACT T& x,T& y)
{x=(T)(a*b);INEXACT T c=(T)(EXPANSION_ARITHMETIC<T>::splitter*a),abig=(T)(c-a);
T ahi=c-abig,alo=a-ahi,err1=x-ahi*bhi,err2=err1-alo*bhi,err3=err2-ahi*blo;y=(alo*blo)-err3;}
}
//#####################################################################
// Function Initialize
//##################################################################### 
template<class T> void EXPANSION_ARITHMETIC<T>::
Check_IEEE_Compliance()
{
  INEXACT T epsilon=(T)1;
  bool every_other=true;
  T half=(T).5,splitter_debug=(T)1,check=(T)1,lastcheck;

  do{lastcheck=check;epsilon*=half;if(every_other) splitter_debug*=(T)2;every_other=!every_other;check=(T)(1+epsilon);}
  while((check!=(T)1)&&(check!=lastcheck));
  splitter_debug+=(T)1;if(splitter_debug!=splitter) PHYSBAM_FATAL_ERROR();
}
//#####################################################################
// Function Add_Expansions
//##################################################################### 
template<class T> void EXPANSION_ARITHMETIC<T>::
Add_Expansions(const ARRAY<T>& e,const ARRAY<T>& f,ARRAY<T>& h)
{
  INEXACT T Qnew,hh;
  T enow=e(1),fnow=f(1),Q;
  int eindex=1,findex=1;

  h.Remove_All();h.Preallocate(e.m+f.m);
  if((fnow>enow)==(fnow>-enow)){Q=enow;if(++eindex<=e.m)enow=e(eindex);}
  else{Q=fnow;if(++findex<=f.m)fnow=f(findex);}
  if((eindex<=e.m)&&(findex<=f.m)){
      if((fnow>enow)==(fnow>-enow)){Fast_Two_Sum(enow,Q,Qnew,hh);if(++eindex<=e.m)enow=e(eindex);}
      else{Fast_Two_Sum(fnow,Q,Qnew,hh);if(++findex<=f.m)fnow=f(findex);}
      Q=Qnew;
      if(hh!=0) h.Append(const_cast<T&>(hh)); // NOTE: if(hh!=0) is not equivalent to if(hh) (bool)--(float)0=true
      while((eindex<=e.m)&&(findex<=f.m)){
          if((fnow>enow)==(fnow>-enow)){Two_Sum(Q,enow,Qnew,hh);if(++eindex<=e.m)enow=e(eindex);}
          else{Two_Sum(Q,fnow,Qnew,hh);if(++findex<=f.m)fnow=f(findex);}
          Q=Qnew;
          if(hh!=0) h.Append(const_cast<T&>(hh));}}
  while(eindex<=e.m){
      Two_Sum(Q,enow,Qnew,hh);
      if(++eindex<=e.m)enow=e(eindex);
      Q=Qnew;
      if(hh!=0) h.Append(const_cast<T&>(hh));}
  while(findex<=f.m){
      Two_Sum(Q,fnow,Qnew,hh);
      if(++findex<=f.m)fnow=f(findex);
      Q=Qnew;
      if(hh!=0) h.Append(const_cast<T&>(hh));}
  if((Q!=0)||!h.m) h.Append(Q);
}
//#####################################################################
// Function Scale_Expansion
//##################################################################### 
template<class T> void EXPANSION_ARITHMETIC<T>::
Scale_Expansion(const ARRAY<T>& e,const T b,ARRAY<T>& h)
{
  INEXACT T Q,sum,product1;
  T hh,product0,enow,bhi,blo;
  int eindex;

  h.Remove_All();h.Preallocate(2*e.m);
  Split(b,bhi,blo);
  Two_Product_Presplit(e(1),b,bhi,blo,Q,hh);
  if(hh!=0) h.Append(hh);
  for(eindex=2;eindex<=e.m;eindex++){
      enow=e(eindex);
      Two_Product_Presplit(enow,b,bhi,blo,product1,product0);
      Two_Sum(Q,product0,sum,hh);
      if(hh!=0) h.Append(hh);
      Fast_Two_Sum(product1,sum,Q,hh);
      if(hh!=0) h.Append(hh);}
  if((Q!=0)||!h.m) h.Append(const_cast<T&>(Q));
}
//#####################################################################
// Function Compress_Expansion
//##################################################################### 
template<class T> void EXPANSION_ARITHMETIC<T>::
Compress_Expansion(ARRAY<T>& e)
{
  INEXACT T Qnew;
  T Q,q,enow, hnow;
  int eindex, hindex,top=1;

  ARRAY<T> h;h.Preallocate(e.m);
  Q=e.Last();
  for(eindex=e.m-1;eindex>=1;eindex--){
      enow=e(eindex);
      Fast_Two_Sum(Q,enow,Qnew,q);
      if(q!=0){h.Append(const_cast<T&>(Qnew));Q=q;} else Q=Qnew;}
  ARRAYS_COMPUTATIONS::Reverse_In_Place(h);
  for(hindex=1;hindex<=h.m;hindex++){
      hnow=h(hindex);
      Fast_Two_Sum(hnow,Q,Qnew,q);
      if(q!=0) h(top++)=q;
      Q=Qnew;}
  h.Resize(top);
  h(top)=Q;
  e.Exchange(h);
}
//#####################################################################
template class EXPANSION_ARITHMETIC<float>;
template class EXPANSION_ARITHMETIC<double>;
