//#####################################################################
// Copyright 2004, Frank Losasso, Andrew Selle.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class RENDERING_ETHER
//#####################################################################
#ifndef __RENDERING_ETHER__
#define __RENDERING_ETHER__

#include <PhysBAM_Geometry/Basic_Geometry/PLANE.h>
#include <PhysBAM_Rendering/PhysBAM_Ray_Tracing/Rendering_Objects/RENDERING_OBJECT.h>
namespace PhysBAM{

template<class T>
class RENDERING_ETHER:public RENDERING_OBJECT<T>
{
    using RENDERING_OBJECT<T>::Inside;using RENDERING_OBJECT<T>::Intersection; // silence -Woverloaded-virtual
public:
    using RENDERING_OBJECT<T>::priority;using RENDERING_OBJECT<T>::name;using RENDERING_OBJECT<T>::support_transparent_overlapping_objects;
    
    RENDERING_ETHER()
    {
        priority=-1;
        name="Rendering_Ether";
        support_transparent_overlapping_objects=true;
    }

    bool Inside(const VECTOR<T,3>& location) const PHYSBAM_OVERRIDE
    {return true;}

    bool Outside(const VECTOR<T,3>& location) const PHYSBAM_OVERRIDE
    {return false;}

    bool Intersection(RAY<VECTOR<T,3> >& ray) const PHYSBAM_OVERRIDE
    {return false;}

    TRIANGULATED_SURFACE<T>* Generate_Triangles() const PHYSBAM_OVERRIDE
    {return 0;}

//#####################################################################
};   
}
#endif

