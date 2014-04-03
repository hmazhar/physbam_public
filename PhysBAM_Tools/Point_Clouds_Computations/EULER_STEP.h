//#####################################################################
// Copyright 2009, Michael Lentine.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class EULER_STEP
//#####################################################################
#ifndef __EULER_STEP__
#define __EULER_STEP__

#include <PhysBAM_Tools/Arrays/ARRAY.h>
namespace PhysBAM{

namespace POINT_CLOUDS_COMPUTATIONS
{
    template<class TV,class T_INDICES>
    static void Euler_Step(const T_INDICES& indices,ARRAY_VIEW<TV> X,const ARRAY_VIEW<TV> V,const typename TV::SCALAR dt)
    {X.Subset(indices)+=dt*V.Subset(indices);}

    template<class TV,class T_INDICES>
    static void Euler_Step(const T_INDICES& indices,ARRAY_VIEW<TV> V,const ARRAY_VIEW<TV> F,const ARRAY_VIEW<typename TV::SCALAR> mass,const typename TV::SCALAR dt)
    {V.Subset(indices)+=dt/mass.Subset(indices)*F.Subset(indices);}

    template<class TV>
    static void Euler_Step(ARRAY_VIEW<TV> X,const ARRAY_VIEW<TV> V,const typename TV::SCALAR dt)
    {X+=dt*V;}

    template<class TV>
    static void Euler_Step(ARRAY_VIEW<TV> V,const ARRAY_VIEW<TV> F,const ARRAY_VIEW<typename TV::SCALAR> mass,const typename TV::SCALAR dt)
    {V+=dt/mass*F;}
}
}
#endif
