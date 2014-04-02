//#####################################################################
// Copyright 2004, Geoffrey Irving, Frank Losasso.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class DUALCONTOUR_OCTREE  
//##################################################################### 
#ifndef COMPILE_WITHOUT_DYADIC_SUPPORT
#ifndef __DUALCONTOUR_OCTREE__
#define __DUALCONTOUR_OCTREE__

#include <PhysBAM_Tools/Arrays/ARRAY.h>
#include <PhysBAM_Tools/Utilities/NONCOPYABLE.h>
#include <PhysBAM_Tools/Vectors/VECTOR.h>
namespace PhysBAM{

template<class T,int d> class VECTOR;
template<class T> class LEVELSET_OCTREE;
template<class T> class TRIANGULATED_SURFACE;

template<class T>
class DUALCONTOUR_OCTREE:public NONCOPYABLE
{
    typedef VECTOR<T,3> TV;
public:
    ARRAY<VECTOR<int,4> > topology;
    int num_topology;
    ARRAY<TV> geometry;
    ARRAY<TV> normals;
    ARRAY<T> phi_nodes;
    LEVELSET_OCTREE<T>* levelset;

    DUALCONTOUR_OCTREE(LEVELSET_OCTREE<T>* levelset_input)
        :num_topology(0),levelset(levelset_input)
    {
        Dualcontour_Octree();
    }

//#####################################################################
    void Dualcontour_Octree();
    TRIANGULATED_SURFACE<T>* Get_Triangulated_Surface();
//#####################################################################
};
}
#endif
#endif
