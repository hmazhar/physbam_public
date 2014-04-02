//#####################################################################
// Copyright 2007, Geoffrey Irving.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class SOBOL
//#####################################################################
// See Bratley and Fox. 1988. Algorithm 659: Implementing Sobol's quasirandom sequence generator. ACM Trans. Math. Softw. 14, 88-100.
//#####################################################################
#include <PhysBAM_Tools/Math_Tools/integer_log.h>
#include <PhysBAM_Tools/Math_Tools/RANGE.h>
#include <PhysBAM_Tools/Random_Numbers/SOBOL.h>
using namespace PhysBAM;
//#####################################################################
namespace{
const int max_dimension=4;
const VECTOR<int,max_dimension> polynomial_degree(1,2,3,3);
const VECTOR<int,max_dimension> polynomial_value(0,1,1,2); // represents x+1, x^2+x+1, x^3+x+1, x^3+x^2+1
int m_initial[max_dimension][max_dimension]={{1},{1,1},{1,3,7},{1,3,3}}; // TODO: should be const, but gcc 4.0.2 is broken.
}
//#####################################################################
// Constructor
//#####################################################################
template<class TV> SOBOL<TV>::
SOBOL(const RANGE<TV>& box)
    :offset(box.min_corner),scales(box.Edge_Lengths()/(T)((TI)1<<max_bits))
{
    STATIC_ASSERT(d<=max_dimension);

    // compute direction numbers
    v.Resize(max_bits);
    for(int i=1;i<=d;i++){
        const int degree=polynomial_degree[i];
        const int polynomial=polynomial_value[i];
        // fill in initial values for m (taken from Numerical Recipes, since according to Bratley and Fox optimal values satisfy complicated conditions
        ARRAY<TI> m(v.m);
        for(int j=1;j<=degree;j++) m(j)=m_initial[i-1][j-1];
        // fill in rest of m using recurrence
        for(int j=degree+1;j<=v.m;j++){
            m(j)=(m(j-degree)<<degree)^m(j-degree);
            for(int k=1;k<degree;k++) 
                if(polynomial&(1<<(k-1))) m(j)^=m(j-k)<<k;}
        // compute direction vectors (stored as Vi * 2^v.m)
        for(int j=1;j<=v.m;j++)
            v(j)[i]=m(j)<<(v.m-j);}

    // start counting
    n=0;
}
//#####################################################################
// Destructor
//#####################################################################
template<class TV> SOBOL<TV>::
~SOBOL()
{}
//#####################################################################
// Function Get_Vector
//#####################################################################
template<class TV> TV SOBOL<TV>::
Get_Vector()
{
    int rightmost_zero_position=1+integer_log_exact(rightmost_bit((int)~n));
    PHYSBAM_ASSERT(rightmost_zero_position<=v.m,"Ran out of bits (this means floating point precision has already been exhausted)");
    const VECTOR<TI,d>& vc=v(rightmost_zero_position);
    for(int i=1;i<=d;i++) x[i]^=vc[i];
    n++;
    return offset+scales*TV(x);
}
//#####################################################################
template class SOBOL<VECTOR<float,1> >;
template class SOBOL<VECTOR<float,2> >;
template class SOBOL<VECTOR<float,3> >;
#ifndef COMPILE_WITHOUT_DOUBLE_SUPPORT
template class SOBOL<VECTOR<double,1> >;
template class SOBOL<VECTOR<double,2> >;
template class SOBOL<VECTOR<double,3> >;
#endif
