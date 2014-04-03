//#####################################################################
// Copyright 2003-2009, Ronald Fedkiw, Geoffrey Irving, Nipun Kwatra, Duc Nguyen, Andrew Selle.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class CUBIC_MN_INTERPOLATION 
//#####################################################################
#ifndef __CUBIC_MN_INTERPOLATION__
#define __CUBIC_MN_INTERPOLATION__

#include <PhysBAM_Tools/Grids_Uniform_Interpolation/INTERPOLATION_UNIFORM.h>
namespace PhysBAM{

template<class T,class T2>
class CUBIC_MN_INTERPOLATION
{
private:
    T m0,m1,m3,n0,n1,n2,n3; // used for acceleration
public:

    CUBIC_MN_INTERPOLATION();
    ~CUBIC_MN_INTERPOLATION();
    void Set_Parameters(const T b=1./3,const T c=1./3);
    T2 Cubic_MN(const T2& u_0,const T2& u_1,const T2& u_2,const T2& u_3,const T alpha) const;

    /*T2 Cubic_MN(const T2& u_0,const T2& u_1,const T2& u_2,const T2& u_3,const T alpha) const
    {T alpha2=alpha*alpha;T alpha3=alpha2*alpha;
    ARRAY<T> weights;weights.Append(-alpha3/2.+alpha2-.5*alpha); //f_k-1
    weights.Append(alpha3/2.-5*alpha2/2.+1); //f_k
    weights.Append(-alpha3/2.+2*alpha2+alpha/2.); //f_k+1
    weights.Append(alpha3/2.-alpha2/2.);  //f_k+2
    return weights(1)*u_0+weights(2)*u_1+weights(3)*u_2+weights(4)*u_3;}*/

    /*ARRAY<T> Cubic_MN_Weights(const T alpha) const
    {T temp=abs(1+alpha);ARRAY<T> weights;weights.Append(((n0*temp+n1)*temp+n2)*temp+n3);
    temp=abs(alpha);weights.Append((m0*temp+m1)*temp*temp+m3);
    temp=abs(1-alpha);weights.Append((m0*temp+m1)*temp*temp+m3);
    temp=abs(2-alpha);weights.Append(((n0*temp+n1)*temp+n2)*temp+n3);
    for(int i=1;i<=weights.m;i++) weights(i)=max((T)0.,min((T)1.,weights(i)));
    return weights;}*/

    ARRAY<T> Cubic_MN_Weights(const T alpha) const;

//#####################################################################
};
}
#endif
