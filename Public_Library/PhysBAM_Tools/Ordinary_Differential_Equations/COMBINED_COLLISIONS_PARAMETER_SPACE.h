//#####################################################################
// Copyright 2009.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class COMBINED_COLLISIONS_PARAMETER_SPACE
//#####################################################################
#ifndef __COMBINED_COLLISIONS_PARAMETER_SPACE__
#define __COMBINED_COLLISIONS_PARAMETER_SPACE__
#include <PhysBAM_Tools/Arrays/ARRAY.h>
#include <PhysBAM_Tools/Nonlinear_Equations/NONLINEAR_FUNCTION.h>
#include <PhysBAM_Tools/Ordinary_Differential_Equations/COMBINED_COLLISIONS.h>
namespace PhysBAM{
template<class TV>
struct COMBINED_COLLISIONS_PARAMETER_SPACE:PARAMETER_SPACE<typename TV::SCALAR>
{
    typedef typename TV::SCALAR T;
    typedef PARAMETER_SPACE<T> P;
    COMBINED_COLLISIONS<TV>& combined_collisions;
    ARRAY<ARRAY<T> > parameters;

    COMBINED_COLLISIONS_PARAMETER_SPACE(COMBINED_COLLISIONS<TV>& combined_collisions_input);
    virtual ~COMBINED_COLLISIONS_PARAMETER_SPACE();

    virtual P& Zero_Clone() const;
    virtual void Op(T a,const P& x,T b,const P& y); // this=a*x+b*y
    virtual void Copy(const P& x); // this=x
    virtual void Zero();
    virtual T Dot(const P& x) const;
};
}
#endif
