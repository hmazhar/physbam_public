//#####################################################################
// Copyright 2009, Jon Gretarsson.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class IMPLICIT_OBJECT_TRANSFORMED
//#####################################################################
#include <PhysBAM_Tools/Vectors/VECTOR_3D.h>
#include <PhysBAM_Geometry/Implicit_Objects_Uniform/MULTIBODY_LEVELSET_IMPLICIT_OBJECT.h>
#include <PhysBAM_Geometry/Registry/STRUCTURE_REGISTRY.h>
//#####################################################################
// Register this class as read-write
//#####################################################################
namespace PhysBAM{
bool Register_Multibody_Levelset_Implicit_Object(){
    static bool done=false;if(done) return true;done=true;
    STRUCTURE_REGISTRY<VECTOR<float,2> >::Register<MULTIBODY_LEVELSET_IMPLICIT_OBJECT<VECTOR<float,2> > >();
    STRUCTURE_REGISTRY<VECTOR<float,3> >::Register<MULTIBODY_LEVELSET_IMPLICIT_OBJECT<VECTOR<float,3> > >();
#ifndef COMPILE_WITHOUT_DOUBLE_SUPPORT
    STRUCTURE_REGISTRY<VECTOR<double,2> >::Register<MULTIBODY_LEVELSET_IMPLICIT_OBJECT<VECTOR<double,2> > >();
    STRUCTURE_REGISTRY<VECTOR<double,3> >::Register<MULTIBODY_LEVELSET_IMPLICIT_OBJECT<VECTOR<double,3> > >();
#endif
    return true;
}
bool registered_multibody_levelset_implicit_object=Register_Multibody_Levelset_Implicit_Object();
}
