//#####################################################################
// Copyright 2004, Andrew Selle.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class PIECEWISE_CONSTANT_PDF
//#####################################################################
#ifndef __PIECEWISE_CONSTANT_PDF__
#define __PIECEWISE_CONSTANT_PDF__

#include <PhysBAM_Tools/Arrays/ARRAY.h>
#include <PhysBAM_Tools/Data_Structures/PAIR.h>
namespace PhysBAM{

template<class T>
class PIECEWISE_CONSTANT_PDF
{
public:
    int m;
    ARRAY<T> pdf; // probability distribution function
    ARRAY<T> cdf; // cumulative distribution function
    T normalization_constant;

    PIECEWISE_CONSTANT_PDF(const int m_input=0)
    {
        Initialize(m_input);
    }

    void Initialize(const int m_input)
    {m=m_input;pdf.Resize(m);cdf.Resize(m+1);}

    void Compute_Cumulative_Distribution_Function()
    {cdf(1)=0;
    for(int i=2;i<=cdf.m;i++)cdf(i)=cdf(i-1)+pdf(i-1); // compute unnormalized cdf
    normalization_constant=cdf(cdf.m); // store pdf normalization constant
    for(int i=1;i<=cdf.m;i++)cdf(i)/=normalization_constant;} // normalize to 1

    PAIR<int,T> Sample(T xi) const // returns (item# , xi')
    {int lower_bound_cdf_index=cdf.Binary_Search(xi)-1;
    if(!lower_bound_cdf_index){return PAIR<int,T>(0,0);}
    T xi_prime=(xi-cdf(lower_bound_cdf_index))/(cdf(lower_bound_cdf_index+1)-cdf(lower_bound_cdf_index));
    return PAIR<int,T>(lower_bound_cdf_index,xi_prime);}

//#####################################################################
};
}
#endif
