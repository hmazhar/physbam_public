//#####################################################################
// Copyright 2009.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class COMBINED_COLLISIONS_NONLINEAR_FUNCTION
//#####################################################################
#ifndef __COMBINED_COLLISIONS_NONLINEAR_FUNCTION__
#define __COMBINED_COLLISIONS_NONLINEAR_FUNCTION__
#include <PhysBAM_Tools/Arrays/ARRAY.h>
#include <PhysBAM_Tools/Nonlinear_Equations/NONLINEAR_FUNCTION.h>
#include <PhysBAM_Tools/Ordinary_Differential_Equations/COMBINED_COLLISIONS.h>
#include <PhysBAM_Tools/Ordinary_Differential_Equations/COMBINED_COLLISIONS_PARAMETER_SPACE.h>
namespace PhysBAM{
template<class TV>
struct COMBINED_COLLISIONS_NONLINEAR_FUNCTION:public NONLINEAR_FUNCTION<typename TV::SCALAR(PARAMETER_SPACE<typename TV::SCALAR>)>
{
    typedef typename TV::SCALAR T;
    typedef PARAMETER_SPACE<T> P;typedef COMBINED_COLLISIONS_PARAMETER_SPACE<TV> CP;
    COMBINED_COLLISIONS<TV>& combined_collisions;
    T dt;
    T time;

    COMBINED_COLLISIONS_NONLINEAR_FUNCTION(COMBINED_COLLISIONS<TV>& combined_collisions_input,T dt_input,T time_input);
    virtual ~COMBINED_COLLISIONS_NONLINEAR_FUNCTION();

//#####################################################################
    virtual T operator()(const P& x) const;
    virtual void Gradient(const P& x,P& g) const;
    void Prepare_State(const CP& x) const;
    void Test_System(const CP* y=0) const;
//#####################################################################
};
}
#endif
