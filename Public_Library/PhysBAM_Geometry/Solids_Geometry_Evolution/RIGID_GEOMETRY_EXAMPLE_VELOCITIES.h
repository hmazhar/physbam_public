//#####################################################################
// Copyright 2002-2008, Robert Bridson, Ronald Fedkiw, Geoffrey Irving, Michael Lentine, Neil Molino, Craig Schroeder, Andrew Selle, Tamar Shinar, Eftychios Sifakis, Joseph Teran.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class RIGID_GEOMETRY_EXAMPLE_VELOCITIES
//#####################################################################
#ifndef __RIGID_GEOMETRY_EXAMPLE_VELOCITIES__
#define __RIGID_GEOMETRY_EXAMPLE_VELOCITIES__

#include <PhysBAM_Tools/Arrays/ARRAYS_FORWARD.h>
namespace PhysBAM{

template<class TV> class TWIST;
template<class TV> class ROTATION;
template<class TV> class FRAME;
template<class TV>
class RIGID_GEOMETRY_EXAMPLE_VELOCITIES
{
    typedef typename TV::SCALAR T;
    typedef typename TV::SPIN T_SPIN;
public:

    RIGID_GEOMETRY_EXAMPLE_VELOCITIES()
    {}

    virtual ~RIGID_GEOMETRY_EXAMPLE_VELOCITIES();

//#####################################################################
    virtual void Set_External_Positions(ARRAY_VIEW<TV> X,ARRAY_VIEW<ROTATION<TV> > rotation,const T time); // set external positions
    virtual void Set_External_Velocities(ARRAY_VIEW<TWIST<TV> > twist,const T velocity_time,const T current_position_time); // velocity and angular velocity
    virtual void Set_Kinematic_Positions(FRAME<TV>& frame,const T time,const int id);
    virtual bool Set_Kinematic_Velocities(TWIST<TV>& twist,const T time,const int id); // Return true if twist was set
//#####################################################################
};
}
#endif
