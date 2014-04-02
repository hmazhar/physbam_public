//#####################################################################
// Copyright 2009.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#include <PhysBAM_Tools/Log/LOG.h>
#include <PhysBAM_Tools/Ordinary_Differential_Equations/COMBINED_COLLISIONS_NONLINEAR_FUNCTION.h>
#include <PhysBAM_Tools/Random_Numbers/RANDOM_NUMBERS.h>
using namespace PhysBAM;
//#####################################################################
// Constructor
//#####################################################################
template<class TV> COMBINED_COLLISIONS_NONLINEAR_FUNCTION<TV>::
COMBINED_COLLISIONS_NONLINEAR_FUNCTION(COMBINED_COLLISIONS<TV>& combined_collisions_input,T dt_input,T time_input)
    :combined_collisions(combined_collisions_input),dt(dt_input),time(time_input)
{
}
//#####################################################################
// Destructor
//#####################################################################
template<class TV> COMBINED_COLLISIONS_NONLINEAR_FUNCTION<TV>::
~COMBINED_COLLISIONS_NONLINEAR_FUNCTION()
{
}
//#####################################################################
// Function operator()
//#####################################################################
template<class TV> typename TV::SCALAR COMBINED_COLLISIONS_NONLINEAR_FUNCTION<TV>::
operator()(const P& x) const
{
    Prepare_State(dynamic_cast<const CP&>(x));
    return combined_collisions.accumulator->Kinetic_Energy(combined_collisions.referenced_list);
}
//#####################################################################
// Function Gradient
//#####################################################################
template<class TV> void COMBINED_COLLISIONS_NONLINEAR_FUNCTION<TV>::
Gradient(const P& x,P& g) const // g = grad(x)
{
    const CP& cx=dynamic_cast<const CP&>(x);
    CP& cg=dynamic_cast<CP&>(g);
    Prepare_State(cx);
    for(int i=1;i<=cg.parameters.m;i++)
        for(int j=1;j<=cg.parameters(i).m;j++)
            cg.parameters(i)(j)=combined_collisions.colliders(i)->Kinetic_Energy_Gradient(j,*combined_collisions.accumulator,cx.parameters(i)(j),dt,time);
}
//#####################################################################
// Function Prepare_State
//#####################################################################
template<class TV> void COMBINED_COLLISIONS_NONLINEAR_FUNCTION<TV>::
Prepare_State(const CP& x) const
{
    combined_collisions.accumulator->Clear(combined_collisions.referenced_list);
    for(int i=1;i<=x.parameters.m;i++)
        for(int j=1;j<=x.parameters(i).m;j++)
            combined_collisions.colliders(i)->Accumulate_Impulse(j,*combined_collisions.accumulator,x.parameters(i)(j),dt,time);
}
//#####################################################################
// Function Test_System
//#####################################################################
template<class TV> void COMBINED_COLLISIONS_NONLINEAR_FUNCTION<TV>::
Test_System(const CP* y) const
{
    RANDOM_NUMBERS<T> random;
    random.Set_Seed(1234);
    CP x(combined_collisions);
    CP g(combined_collisions);
    CP ngl(combined_collisions);
    CP ngr(combined_collisions);
    CP ngc(combined_collisions);
    if(y) x.parameters=y->parameters;
    else
        for(int i=1;i<=x.parameters.m;i++)
            random.Fill_Uniform(x.parameters,(T).1,(T).2);
    T ke=(*this)(x);
    Gradient(x,g);
    T delta=(T)1e-6;
    for(int i=1;i<=x.parameters.m;i++)
        for(int j=1;j<=x.parameters(i).m;j++){
            x.parameters(i)(j)+=delta;
            T r=(*this)(x);
            x.parameters(i)(j)-=delta*2;
            T l=(*this)(x);
            x.parameters(i)(j)+=delta;
            ngl.parameters(i)(j)=(ke-l)/delta;
            ngr.parameters(i)(j)=(r-ke)/delta;
            ngc.parameters(i)(j)=(r-l)/(2*delta);}

    LOG::cout<<"KE: "<<ke<<std::endl;
    LOG::cout<<"Computed: "<<g.parameters<<std::endl;
    LOG::cout<<"Approximated (l): "<<ngl.parameters<<std::endl;
    LOG::cout<<"Approximated (r): "<<ngr.parameters<<std::endl;
    LOG::cout<<"Approximated (c): "<<ngc.parameters<<std::endl;
}
template class COMBINED_COLLISIONS_NONLINEAR_FUNCTION<VECTOR<float,1> >;
template class COMBINED_COLLISIONS_NONLINEAR_FUNCTION<VECTOR<float,2> >;
template class COMBINED_COLLISIONS_NONLINEAR_FUNCTION<VECTOR<float,3> >;
#ifndef COMPILE_WITHOUT_DOUBLE_SUPPORT
template class COMBINED_COLLISIONS_NONLINEAR_FUNCTION<VECTOR<double,1> >;
template class COMBINED_COLLISIONS_NONLINEAR_FUNCTION<VECTOR<double,2> >;
template class COMBINED_COLLISIONS_NONLINEAR_FUNCTION<VECTOR<double,3> >;
#endif
