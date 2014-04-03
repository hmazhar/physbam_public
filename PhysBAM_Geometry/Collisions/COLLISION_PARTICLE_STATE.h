//#####################################################################
// Copyright 2006, Geoffrey Irving.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#ifndef __COLLISION_PARTICLE_STATE__
#define __COLLISION_PARTICLE_STATE__

namespace PhysBAM{

template<class TV>
class COLLISION_PARTICLE_STATE
{
    typedef typename TV::SCALAR T;
public:
    bool enforce;
    TV normal;
    T VN;
    TV VT_body;
    T delta_VN;
    T friction;

    COLLISION_PARTICLE_STATE()
        :enforce(false),VN(0),delta_VN(0),friction(0)
    {}

//#####################################################################
};
}
#endif
