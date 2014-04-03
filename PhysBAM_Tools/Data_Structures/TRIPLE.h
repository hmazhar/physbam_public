//#####################################################################
// Copyright 2003-2007, Geoffrey Irving, Frank Losasso, Neil Molino.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class TRIPLE
//##################################################################### 
#ifndef __TRIPLE__
#define __TRIPLE__

namespace PhysBAM{

template<class T1,class T2,class T3>
class TRIPLE
{
public:
    T1 x;T2 y;T3 z;

    TRIPLE(int input=0) 
        :x(T1()),y(T2()),z(T3())
    {}

    TRIPLE(const T1& x_input,const T2& y_input,const T3& z_input) 
        :x(x_input),y(y_input),z(z_input)
    {}

    bool operator==(const TRIPLE& t) const
    {return x==t.x && y==t.y && z==t.z;}

    bool operator!=(const TRIPLE& t) const
    {return !(*this==t);}

//#####################################################################
};  
template<class T1,class T2,class T3>
inline TRIPLE<T1,T2,T3> Tuple(const T1& x,const T2& y,const T3& z)
{return TRIPLE<T1,T2,T3>(x,y,z);}
}
#endif
