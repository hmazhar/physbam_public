//#####################################################################
// Copyright 2009.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#include <PhysBAM_Tools/Arrays_Computations/DOT_PRODUCT.h>
#include <PhysBAM_Tools/Ordinary_Differential_Equations/COMBINED_COLLISIONS_PARAMETER_SPACE.h>
using namespace PhysBAM;
//#####################################################################
// Constructor
//#####################################################################
template<class TV> COMBINED_COLLISIONS_PARAMETER_SPACE<TV>::
COMBINED_COLLISIONS_PARAMETER_SPACE(COMBINED_COLLISIONS<TV>& combined_collisions_input)
    :combined_collisions(combined_collisions_input),parameters(combined_collisions.collider_info.m)
{
    for(int i=1;i<=combined_collisions.collider_info.m;i++)
        parameters(i).Resize(combined_collisions.collider_info(i).scale.m);
}
//#####################################################################
// Destructor
//#####################################################################
template<class TV> COMBINED_COLLISIONS_PARAMETER_SPACE<TV>::
~COMBINED_COLLISIONS_PARAMETER_SPACE()
{
}
//#####################################################################
// Function Zero_Clone
//#####################################################################
template<class TV> PARAMETER_SPACE<typename TV::SCALAR>& COMBINED_COLLISIONS_PARAMETER_SPACE<TV>::
Zero_Clone() const
{
    return *new COMBINED_COLLISIONS_PARAMETER_SPACE(combined_collisions);
}
//#####################################################################
// Function Op
//#####################################################################
template<class TV> void COMBINED_COLLISIONS_PARAMETER_SPACE<TV>::
Op(T a,const P& x,T b,const P& y) // this=a*x+b*y
{
    const COMBINED_COLLISIONS_PARAMETER_SPACE<TV>& cx=dynamic_cast<const COMBINED_COLLISIONS_PARAMETER_SPACE<TV>&>(x);
    const COMBINED_COLLISIONS_PARAMETER_SPACE<TV>& cy=dynamic_cast<const COMBINED_COLLISIONS_PARAMETER_SPACE<TV>&>(y);
    for(int i=1;i<=parameters.m;i++)
        parameters(i)=a*cx.parameters(i)+b*cy.parameters(i);
}
//#####################################################################
// Function Copy
//#####################################################################
template<class TV> void COMBINED_COLLISIONS_PARAMETER_SPACE<TV>::
Copy(const P& x) // this=x
{
    const COMBINED_COLLISIONS_PARAMETER_SPACE<TV>& cx=dynamic_cast<const COMBINED_COLLISIONS_PARAMETER_SPACE<TV>&>(x);
    parameters=cx.parameters;
}
//#####################################################################
// Function Zero
//#####################################################################
template<class TV> void COMBINED_COLLISIONS_PARAMETER_SPACE<TV>::
Zero()
{
    for(int i=1;i<=parameters.m;i++)
        ARRAYS_COMPUTATIONS::Fill(parameters(i),(T)0);
}
//#####################################################################
// Function Dot
//#####################################################################
template<class TV> typename TV::SCALAR COMBINED_COLLISIONS_PARAMETER_SPACE<TV>::
Dot(const P& x) const
{
    const COMBINED_COLLISIONS_PARAMETER_SPACE<TV>& cx=dynamic_cast<const COMBINED_COLLISIONS_PARAMETER_SPACE<TV>&>(x);
    T dot=0;
    for(int i=1;i<=parameters.m;i++)
        dot+=ARRAYS_COMPUTATIONS::Dot_Product(parameters(i),cx.parameters(i));
    return dot;
}
template class COMBINED_COLLISIONS_PARAMETER_SPACE<VECTOR<float,1> >;
template class COMBINED_COLLISIONS_PARAMETER_SPACE<VECTOR<float,2> >;
template class COMBINED_COLLISIONS_PARAMETER_SPACE<VECTOR<float,3> >;
#ifndef COMPILE_WITHOUT_DOUBLE_SUPPORT
template class COMBINED_COLLISIONS_PARAMETER_SPACE<VECTOR<double,1> >;
template class COMBINED_COLLISIONS_PARAMETER_SPACE<VECTOR<double,2> >;
template class COMBINED_COLLISIONS_PARAMETER_SPACE<VECTOR<double,3> >;
#endif
