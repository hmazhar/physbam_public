//#####################################################################
// Copyright 2002-2005, Ronald Fedkiw, Geoffrey Irving, Joseph Teran.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class HEXAHEDRON
//##################################################################### 
#ifndef __HEXAHEDRON__
#define __HEXAHEDRON__

#include <PhysBAM_Tools/Vectors/VECTOR_3D.h>
namespace PhysBAM{

template<class T>
class HEXAHEDRON
{
public:
    VECTOR<T,3> x1,x2,x3,x4,x5,x6,x7,x8; 
    
    HEXAHEDRON()
        :x1(-1,-1,-1),x2(-1,-1,1),x3(-1,1,-1),x4(-1,1,1),x5(1,-1,-1),x6(1,-1,1),x7(1,1,-1),x8(1,1,1)
    {}

    HEXAHEDRON(const VECTOR<T,3>& x1_input,const VECTOR<T,3>& x2_input,const VECTOR<T,3>& x3_input,const VECTOR<T,3>& x4_input,
        const VECTOR<T,3>& x5_input,const VECTOR<T,3>& x6_input,const VECTOR<T,3>& x7_input,const VECTOR<T,3>& x8_input)
        :x1(x1_input),x2(x2_input),x3(x3_input),x4(x4_input),x5(x5_input),x6(x6_input),x7(x7_input),x8(x8_input)
    {}

    static T Volume(const VECTOR<T,3>& x1,const VECTOR<T,3>& x2,const VECTOR<T,3>& x3,const VECTOR<T,3>& x4,
                    const VECTOR<T,3>& x5,const VECTOR<T,3>& x6,const VECTOR<T,3>& x7,const VECTOR<T,3>& x8)
    {return abs(Signed_Volume(x1,x2,x3,x4,x5,x6,x7,x8));}

    static T Signed_Volume(const VECTOR<T,3>& x1,const VECTOR<T,3>& x2,const VECTOR<T,3>& x3,const VECTOR<T,3>& x4,
                           const VECTOR<T,3>& x5,const VECTOR<T,3>& x6,const VECTOR<T,3>& x7,const VECTOR<T,3>& x8)
    {VECTOR<T,3> x28=x8-x2,x23=x3-x2,x24=x4-x2,x35=x5-x3,x58=x8-x3,x57=x7-x5,x25=x5-x2,x26=x6-x2,x15=x5-x1;
    return (T)one_sixth*(VECTOR<T,3>::Triple_Product(x28,x23,x24)-VECTOR<T,3>::Triple_Product(x35,x58,x57)+VECTOR<T,3>::Triple_Product(x25,x28,x26)
        +VECTOR<T,3>::Triple_Product(x23,x28,x25)-VECTOR<T,3>::Triple_Product(x25,x35,x15));}

//#####################################################################
    T Volume() const;
    T Signed_Volume() const;
//#####################################################################
};   
}
#endif

