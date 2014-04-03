//#####################################################################
// Copyright 2004-2008, Ron Fedkiw, Geoffrey Irving, Michael Lentine, Craig Schroeder, Andrew Selle, Tamar Shinar, Eftychios Sifakis.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class RIGID_GEOMETRY_EVOLUTION_PARAMETERS
//#####################################################################
#include <PhysBAM_Tools/Read_Write/Utilities/FILE_UTILITIES.h>
#include <PhysBAM_Geometry/Solids_Geometry/RIGID_GEOMETRY_EVOLUTION_PARAMETERS.h>
namespace PhysBAM{
//#####################################################################
// Constructor
//#####################################################################
RIGID_GEOMETRY_EVOLUTION_PARAMETERS::
RIGID_GEOMETRY_EVOLUTION_PARAMETERS()
    :use_kinematic_keyframes(true)
{
}
//#####################################################################
// Destructor
//#####################################################################
RIGID_GEOMETRY_EVOLUTION_PARAMETERS::
~RIGID_GEOMETRY_EVOLUTION_PARAMETERS()
{
}
//#####################################################################
}
