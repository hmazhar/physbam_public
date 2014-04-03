//#####################################################################
// Copyright 2002-2006, Robert Bridson, Ronald Fedkiw, Frederic Gibou, Geoffrey Irving.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class QUADRATIC
//#####################################################################
#ifndef __QUADRATIC__
#define __QUADRATIC__

#include <PhysBAM_Tools/Math_Tools/exchange.h>
#include <PhysBAM_Tools/Math_Tools/sqr.h>
#include <PhysBAM_Tools/Nonlinear_Equations/NONLINEAR_FUNCTION.h>
namespace PhysBAM{

template<class T>
class QUADRATIC:public NONLINEAR_FUNCTION<T(T)>
{
public:
    T a,b,c; // coefficients
    int roots; // number of roots, -1 indicates a=b=c=0 - always a root!
    T root1,root2; // root1 < root2

    QUADRATIC();
    QUADRATIC(const T a_input,const T b_input,const T c_input);

    T operator()(const T x) const PHYSBAM_OVERRIDE
    {return (a*x+b)*x+c;}

    T Discriminant() const
    {return sqr(b)-4*a*c;}

    void Coefficients_From_Interpolation(T x1,T y1,T x2,T y2,T x3,T y3);
    void Compute_Roots();
    void Compute_Roots_In_Interval(const T xmin,const T xmax);
//#####################################################################
};
}
#endif
