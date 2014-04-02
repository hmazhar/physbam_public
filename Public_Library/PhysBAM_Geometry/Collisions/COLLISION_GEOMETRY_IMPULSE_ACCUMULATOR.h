//#####################################################################
// Copyright 2006, Geoffrey Irving.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class COLLISION_GEOMETRY_IMPULSE_ACCUMULATOR
//#####################################################################
#ifndef __COLLISION_GEOMETRY_IMPULSE_ACCUMULATOR__
#define __COLLISION_GEOMETRY_IMPULSE_ACCUMULATOR__

#include <PhysBAM_Tools/Matrices/MATRIX_POLICY.h>
#include <PhysBAM_Tools/Utilities/NONCOPYABLE.h>
#include <PhysBAM_Tools/Utilities/PHYSBAM_OVERRIDE.h>
namespace PhysBAM{

template<class TV>
class COLLISION_GEOMETRY_IMPULSE_ACCUMULATOR:public NONCOPYABLE
{
private:
    typedef typename TV::SPIN T_SPIN;
public:

    virtual ~COLLISION_GEOMETRY_IMPULSE_ACCUMULATOR()
    {}

//#####################################################################
    virtual void Reset()=0;
    virtual void Add_Impulse(const TV& location,const TWIST<TV>& impulse)=0;
//#####################################################################
};
}
#endif
