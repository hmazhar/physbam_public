//#####################################################################
// Copyright 2005, Eran Guendelman, Geoffrey Irving.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class RLE_LINEAR_PROFILE
//#####################################################################
#ifndef COMPILE_WITHOUT_RLE_SUPPORT
#ifndef __RLE_LINEAR_PROFILE__
#define __RLE_LINEAR_PROFILE__

#include <PhysBAM_Tools/Matrices/SYMMETRIC_MATRIX_2X2.h>
#include <PhysBAM_Tools/Vectors/VECTOR_2D.h>
namespace PhysBAM{

template<class T>
class RLE_LINEAR_PROFILE
{
public:
    T zeroth_moment,first_moment;

    RLE_LINEAR_PROFILE()
        :zeroth_moment(0),first_moment(0)
    {}

    void Update_Constant_Profile(const int len,const T scale,T& u)
    {u+=scale*zeroth_moment/(len+1);}

    void Update_Linear_Profile(const int dj,const int len,const T scale,T& u_lo,T& u_hi) // u_hi unused if len==0
    {if(len) Update_Long_Linear_Profile(dj,len,scale,u_lo,u_hi);
    else u_lo+=scale*zeroth_moment;}

    void Update_Long_Linear_Profile(const int dj,const int len,const T scale,T& u_lo,T& u_hi)
    {assert(len>0);SYMMETRIC_MATRIX<T,2> normal_matrix_inverse((T)(len*(2*len+1)),(T)(-3*len),(T)6);T matrix_scale=(T)2/(len*(len+1)*(len+2));
    VECTOR<T,2> change=scale*matrix_scale*(normal_matrix_inverse*VECTOR<T,2>(zeroth_moment,first_moment-dj*zeroth_moment));
    u_lo+=change.x;u_hi+=change.x+len*change.y;}

    void Add_Single_Value(const int offset,const T u)
    {zeroth_moment+=u;first_moment+=offset*u;}

    void Add_Constant_Profile(const int offset,const int len,const T u) // changes applied to offset, offset+1, ..., offset+len
    {T constant=(len+1)*u;zeroth_moment+=constant;first_moment+=(offset+(T).5*len)*constant;}

    void Add_Linear_Profile(const int offset,const int len,const T ustart,const T uslope)
    {T constant=(len+1)*(ustart+(T).5*len*uslope);zeroth_moment+=constant;
    first_moment+=offset*constant+(T)one_sixth*len*(len+1)*(3*ustart+(2*len+1)*uslope);}

    void Add_Linear_Product_Profile(const int offset,const int len,const T u1start,const T u1slope,const T u2start,const T u2slope)
    {T c1=u1start*u2start,c2=u1start*u2slope+u2start*u1slope,c3=u1slope*u2slope;
    T constant=(T)one_sixth*(len+1)*(6*c1+3*len*c2+len*(2*len+1)*c3);zeroth_moment+=constant;
    first_moment+=offset*constant+(T)one_twelfth*len*(len+1)*(6*c1+2*(2*len+1)*c2+3*len*(len+1)*c3);}

//#####################################################################
};
}
#endif
#endif
