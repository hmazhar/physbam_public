//#####################################################################
// Copyright 2004-2008, Ron Fedkiw, Geoffrey Irving, Michael Lentine, Craig Schroeder, Andrew Selle, Tamar Shinar, Eftychios Sifakis, Jonathan Su, Joseph Teran.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class RIGID_GEOMETRY_EVOLUTION_PARAMETERS
//#####################################################################
#ifndef __RIGID_GEOMETRY_EVOLUTION_PARAMETERS__
#define __RIGID_GEOMETRY_EVOLUTION_PARAMETERS__

#include <PhysBAM_Tools/Utilities/NONCOPYABLE.h>
namespace PhysBAM{

class RIGID_GEOMETRY_EVOLUTION_PARAMETERS:public NONCOPYABLE
{
public:
    bool use_kinematic_keyframes;

    RIGID_GEOMETRY_EVOLUTION_PARAMETERS();
    virtual ~RIGID_GEOMETRY_EVOLUTION_PARAMETERS();
//#####################################################################
};
}
#endif
