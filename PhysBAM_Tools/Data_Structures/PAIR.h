//#####################################################################
// Copyright 2003-2008, Geoffrey Irving, Michael Lentine, Frank Losasso, Joseph Teran.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class PAIR
//##################################################################### 
#ifndef __PAIR__
#define __PAIR__

#include <PhysBAM_Tools/Math_Tools/choice.h>
namespace PhysBAM{

template<class T1,class T2>
class PAIR
{
public:
    template<int i> struct GET:public IF<i==1,T1,T2>{};

    T1 x;T2 y;

    PAIR() 
        :x(T1()),y(T2()) 
    {}

    PAIR(const T1& x_input,const T2& y_input) 
        :x(x_input),y(y_input) 
    {}

    ~PAIR()
    {}

    bool operator==(const PAIR& p) const
    {return x==p.x && y==p.y;}

    bool operator!=(const PAIR& p) const
    {return !(*this==p);}

    bool operator<(const PAIR& p) const
    {return x<p.x || (x==p.x && y<p.y);}
    
    bool operator>(const PAIR& p) const
    {return x>p.x || (x==p.x && y>p.y);}

    template<int i> typename GET<i>::TYPE& Get()
    {return choice<i>(x,y);}

    template<int i> const typename GET<i>::TYPE& Get() const
    {return choice<i>(x,y);}

    void Get(T1& a,T2& b) const
    {a=x;b=y;}

//#####################################################################
};  
template<class T1,class T2>
inline PAIR<T1,T2> Tuple(const T1& x,const T2& y)
{return PAIR<T1,T2>(x,y);}
}
#endif
