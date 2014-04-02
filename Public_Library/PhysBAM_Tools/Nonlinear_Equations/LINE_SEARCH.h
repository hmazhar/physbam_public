//#####################################################################
// Copyright 2009.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class LINE_SEARCH
//##################################################################### 
#ifndef __LINE_SEARCH__
#define __LINE_SEARCH__
namespace PhysBAM{
template<class F> class NONLINEAR_FUNCTION;

template<class T>
struct LINE_SEARCH
{
    struct BRACKET
    {
        T a,m,b;
        T Fa,Fm,Fb;
    };

    void Update_Interval(BRACKET& s,T d,T Fd) const; // a < b < c, update with new point d, where a < d < c.
    int Compute_Quadratic_Minimum(const BRACKET& s,T& x,T tolerance) const; // 0=terminate, 1=fail, 2=good
    bool Line_Search_Quadratic_Golden_Section(NONLINEAR_FUNCTION<T(T)>& F,T a,T b,T& x,int max_iterations,T interval_tolerance,T quadratic_tolerance) const;
    bool Line_Search_Golden_Section(NONLINEAR_FUNCTION<T(T)>& F,T a,T b,T& x,int max_iterations,T interval_tolerance) const;
    T Best_Value(const BRACKET& s) const;
};
}
#endif
