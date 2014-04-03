//#####################################################################
// Copyright 2006-2009, Zhaosheng Bao, Ron Fedkiw, Geoffrey Irving, Sergey Koltakov, Nipun Kwatra, Michael Lentine, Craig Schroeder, Andrew Selle, Tamar Shinar, Eftychios Sifakis, Joseph Teran.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#include <PhysBAM_Geometry/Solids_Geometry/RIGID_GEOMETRY.h>
#include <PhysBAM_Geometry/Solids_Geometry_Evolution/KINEMATIC_EVOLUTION.h>
#include <PhysBAM_Geometry/Solids_Geometry_Evolution/RIGID_GEOMETRY_EXAMPLE_VELOCITIES.h>
using namespace PhysBAM;
template<class TV> KINEMATIC_EVOLUTION<TV>::
KINEMATIC_EVOLUTION(RIGID_GEOMETRY_COLLECTION<TV>& rigid_geometry_collection_input,bool use_kinematic_keyframes_input)
    :rigid_geometry_collection(rigid_geometry_collection_input),use_kinematic_keyframes(use_kinematic_keyframes_input)
{
}
template<class TV> KINEMATIC_EVOLUTION<TV>::
~KINEMATIC_EVOLUTION()
{
}
//#####################################################################
// Function Set_External_Velocities
//#####################################################################
template<class TV> void KINEMATIC_EVOLUTION<TV>::
Set_External_Velocities(ARRAY_VIEW<TWIST<TV> > twist,const T velocity_time,const T current_position_time)
{
    for(int i=1;i<=rigid_geometry_collection.kinematic_rigid_geometry.m;i++){int p=rigid_geometry_collection.kinematic_rigid_geometry(i);Set_External_Velocities(twist(p),velocity_time,p);}
    rigid_geometry_collection.rigid_geometry_example_velocities->Set_External_Velocities(twist,velocity_time,current_position_time);
}
//#####################################################################
// Function Set_External_Velocities
//#####################################################################
template<class TV> void KINEMATIC_EVOLUTION<TV>::
Set_External_Velocities(TWIST<TV>& twist,const T time,const int id)
{
    RIGID_GEOMETRY<TV>& rigid_geometry=rigid_geometry_collection.Rigid_Geometry(id);
    if(rigid_geometry.is_static){twist=TWIST<TV>();return;}
    if(!use_kinematic_keyframes){Set_Kinematic_Velocities(twist,(T)1e-3,time,id);return;} // Use 1e-3 for backward differencing if the example did not implement velocity calculations.
    RIGID_GEOMETRY_STATE<TV> interpolated_state;
    rigid_geometry.Interpolate_Between_States(kinematic_current_state(rigid_geometry.particle_index),kinematic_next_state(rigid_geometry.particle_index),time,interpolated_state);
    twist=interpolated_state.twist;
}
//#####################################################################
// Function Set_Kinematic_Velocities
//#####################################################################
template<class TV> void KINEMATIC_EVOLUTION<TV>::
Set_Kinematic_Velocities(TWIST<TV>& twist,const T frame_dt,const T time,const int id)
{
    RIGID_GEOMETRY<TV>* rigid_geometry=&rigid_geometry_collection.Rigid_Geometry(id);
    if(rigid_geometry_collection.rigid_geometry_example_velocities->Set_Kinematic_Velocities(twist,time,id)) return;
    RIGID_GEOMETRY_STATE<TV> previous_state,current_state;previous_state.time=time-frame_dt;current_state.time=time;
    int new_id=id;
    rigid_geometry_collection.rigid_geometry_example_velocities->Set_Kinematic_Positions(previous_state.frame,previous_state.time,new_id);
    rigid_geometry_collection.rigid_geometry_example_velocities->Set_Kinematic_Positions(current_state.frame,current_state.time,new_id);
    rigid_geometry->Compute_Velocity_Between_States(previous_state,current_state,current_state);
    twist=current_state.twist;
}
//#####################################################################
// Function Get_Current_Kinematic_Keyframes
//#####################################################################
template<class TV> void KINEMATIC_EVOLUTION<TV>::
Get_Current_Kinematic_Keyframes(const T dt,const T time)
{
    if(!use_kinematic_keyframes) return;
    kinematic_current_state.Remove_All();kinematic_next_state.Remove_All();
    kinematic_current_state.Resize(rigid_geometry_collection.particles.array_collection->Size());
    kinematic_next_state.Resize(rigid_geometry_collection.particles.array_collection->Size());
    for(int i=1;i<=rigid_geometry_collection.kinematic_rigid_geometry.m;i++){int p=rigid_geometry_collection.kinematic_rigid_geometry(i);
        kinematic_current_state(p).time=time;kinematic_next_state(p).time=time+dt;
        rigid_geometry_collection.rigid_geometry_example_velocities->Set_Kinematic_Positions(kinematic_current_state(p).frame,time,p);
        rigid_geometry_collection.rigid_geometry_example_velocities->Set_Kinematic_Positions(kinematic_next_state(p).frame,time+dt,p);
        Set_Kinematic_Velocities(kinematic_current_state(p).twist,dt,time,p);
        Set_Kinematic_Velocities(kinematic_next_state(p).twist,dt,time+dt,p);}
}
//#####################################################################
// Function Set_External_Positions
//#####################################################################
template<class TV> void KINEMATIC_EVOLUTION<TV>::
Set_External_Positions(TV& X,ROTATION<TV>& rotation,const T time,const int id)
{
    RIGID_GEOMETRY<TV>* rigid_geometry=&rigid_geometry_collection.Rigid_Geometry(id);
    int new_id=id;
    int index=rigid_geometry->particle_index;
    if(!use_kinematic_keyframes){
        FRAME<TV> frame;rigid_geometry_collection.rigid_geometry_example_velocities->Set_Kinematic_Positions(frame,time,new_id);
        X=frame.t;rotation=frame.r;return;}
    RIGID_GEOMETRY_STATE<TV> interpolated_state;
    rigid_geometry->Interpolate_Between_States(kinematic_current_state(index),kinematic_next_state(index),time,interpolated_state);
    X=interpolated_state.frame.t;rotation=interpolated_state.frame.r;
}
//#####################################################################
// Function Set_External_Positions
//#####################################################################
template<class TV> void KINEMATIC_EVOLUTION<TV>::
Set_External_Positions(ARRAY_VIEW<TV> X,ARRAY_VIEW<ROTATION<TV> > rotation,const T time)
{
    for(int i=1;i<=rigid_geometry_collection.kinematic_rigid_geometry.m;i++){
        int p=rigid_geometry_collection.kinematic_rigid_geometry(i);Set_External_Positions(X(p),rotation(p),time,p);}
    rigid_geometry_collection.rigid_geometry_example_velocities->Set_External_Positions(X,rotation,time);
}
//#####################################################################
// Function Reset_Kinematic_Rigid_Bodies
//#####################################################################
template<class TV> void KINEMATIC_EVOLUTION<TV>::
Reset_Kinematic_Rigid_Bodies(const T time)
{
    // Move kinematic bodies to their position at given time
    for(int i=1;i<=rigid_geometry_collection.kinematic_rigid_geometry.m;i++){
        RIGID_GEOMETRY<TV>& rigid_geometry=rigid_geometry_collection.Rigid_Geometry(rigid_geometry_collection.kinematic_rigid_geometry(i));
        Set_External_Positions(rigid_geometry.X(),rigid_geometry.Rotation(),time,rigid_geometry.particle_index);
        Set_External_Velocities(rigid_geometry.Twist(),time,rigid_geometry.particle_index);}
}
//#####################################################################
template class KINEMATIC_EVOLUTION<VECTOR<float,1> >;
template class KINEMATIC_EVOLUTION<VECTOR<float,2> >;
template class KINEMATIC_EVOLUTION<VECTOR<float,3> >;
#ifndef COMPILE_WITHOUT_DOUBLE_SUPPORT
template class KINEMATIC_EVOLUTION<VECTOR<double,1> >;
template class KINEMATIC_EVOLUTION<VECTOR<double,2> >;
template class KINEMATIC_EVOLUTION<VECTOR<double,3> >;
#endif
