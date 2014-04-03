//#####################################################################
// Copyright 2002-2008, Robert Bridson, Ronald Fedkiw, Geoffrey Irving, Michael Lentine, Neil Molino, Craig Schroeder, Andrew Selle, Tamar Shinar, Eftychios Sifakis, Joseph Teran.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#include <PhysBAM_Tools/Arrays/ARRAY_VIEW.h>
#include <PhysBAM_Tools/Log/DEBUG_UTILITIES.h>
#include <PhysBAM_Tools/Vectors/VECTOR_3D.h>
#include <PhysBAM_Geometry/Solids_Geometry_Evolution/RIGID_GEOMETRY_EXAMPLE_VELOCITIES.h>
using namespace PhysBAM;
//#####################################################################
// Destructor
//#####################################################################
template<class TV> RIGID_GEOMETRY_EXAMPLE_VELOCITIES<TV>::
~RIGID_GEOMETRY_EXAMPLE_VELOCITIES()
{}
//#####################################################################
// Function Set_External_Positions
//#####################################################################
template<class TV> void RIGID_GEOMETRY_EXAMPLE_VELOCITIES<TV>::
Set_External_Positions(ARRAY_VIEW<TV> X,ARRAY_VIEW<ROTATION<TV> > rotation,const T time)
{
    PHYSBAM_WARN_IF_NOT_OVERRIDDEN();
}
//#####################################################################
// Function Set_External_Velocities
//#####################################################################
template<class TV> void RIGID_GEOMETRY_EXAMPLE_VELOCITIES<TV>::
Set_External_Velocities(ARRAY_VIEW<TWIST<TV> > twist,const T velocity_time,const T current_position_time)
{
    PHYSBAM_WARN_IF_NOT_OVERRIDDEN();
}
//#####################################################################
// Function Set_Kinematic_Positions
//#####################################################################
template<class TV> void RIGID_GEOMETRY_EXAMPLE_VELOCITIES<TV>::
Set_Kinematic_Positions(FRAME<TV>& frame,const T time,const int id)
{
    PHYSBAM_WARN_IF_NOT_OVERRIDDEN();
}
//#####################################################################
// Function Set_Kinematic_Velocities
//#####################################################################
template<class TV> bool RIGID_GEOMETRY_EXAMPLE_VELOCITIES<TV>::
Set_Kinematic_Velocities(TWIST<TV>& twist,const T time,const int id)
{
    PHYSBAM_WARN_IF_NOT_OVERRIDDEN();
    return false;
}
//#####################################################################
template class RIGID_GEOMETRY_EXAMPLE_VELOCITIES<VECTOR<float,1> >;
template class RIGID_GEOMETRY_EXAMPLE_VELOCITIES<VECTOR<float,2> >;
template class RIGID_GEOMETRY_EXAMPLE_VELOCITIES<VECTOR<float,3> >;
#ifndef COMPILE_WITHOUT_DOUBLE_SUPPORT
template class RIGID_GEOMETRY_EXAMPLE_VELOCITIES<VECTOR<double,1> >;
template class RIGID_GEOMETRY_EXAMPLE_VELOCITIES<VECTOR<double,2> >;
template class RIGID_GEOMETRY_EXAMPLE_VELOCITIES<VECTOR<double,3> >;
#endif
