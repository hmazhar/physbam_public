//#####################################################################
// Copyright 2006-2008, Zhaosheng Bao, Ron Fedkiw, Geoffrey Irving, Sergey Koltakov, Michael Lentine, Craig Schroeder, Andrew Selle, Tamar Shinar, Eftychios Sifakis, Joseph Teran.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class KINEMATIC_EVOLUTION
//#####################################################################
#ifndef __KINEMATIC_EVOLUTION__
#define __KINEMATIC_EVOLUTION__

#include <PhysBAM_Tools/Arrays/ARRAYS_FORWARD.h>
#include <PhysBAM_Tools/Matrices/MATRIX_POLICY.h>
namespace PhysBAM{

template<class TV> class RIGID_GEOMETRY_COLLECTION;
template<class TV> class RIGIDS_EXAMPLE_FORCES_AND_VELOCITIES;
template<class TV> class RIGID_GEOMETRY_STATE;

template<class TV>
class KINEMATIC_EVOLUTION
{
    typedef typename TV::SCALAR T;
    typedef typename TV::SPIN T_SPIN;
    enum WORKAROUND {d=TV::m};
public:
    RIGID_GEOMETRY_COLLECTION<TV>& rigid_geometry_collection;
    ARRAY<RIGID_GEOMETRY_STATE<TV> > kinematic_current_state; // TODO: These are sparse; compact them
    ARRAY<RIGID_GEOMETRY_STATE<TV> > kinematic_next_state;
    bool use_kinematic_keyframes;

    KINEMATIC_EVOLUTION(RIGID_GEOMETRY_COLLECTION<TV>& rigid_body_collection_input,bool use_kinematic_keyframes_input);
    virtual ~KINEMATIC_EVOLUTION();

//#####################################################################
    void Set_External_Velocities(ARRAY_VIEW<TWIST<TV> > twist,const T velocity_time,const T current_position_time);
    virtual void Set_External_Velocities(TWIST<TV>& twist,const T time,const int id);
    virtual void Set_Kinematic_Velocities(TWIST<TV>& twist,const T frame_dt,const T time,const int id);
    void Set_External_Positions(ARRAY_VIEW<TV> X,ARRAY_VIEW<ROTATION<TV> > rotation,const T time);
    virtual void Set_External_Positions(TV& X,ROTATION<TV>& rotation,const T time,const int id);
    void Get_Current_Kinematic_Keyframes(const T dt,const T time);
    void Reset_Kinematic_Rigid_Bodies(const T time);
//#####################################################################
};
}
#endif
