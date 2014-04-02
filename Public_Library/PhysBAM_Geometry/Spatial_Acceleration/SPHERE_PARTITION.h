//#####################################################################
// Copyright 2002-2005, Ronald Fedkiw, Geoffrey Irving.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class SPHERE_PARTITION
//##################################################################### 
#ifndef __SPHERE_PARTITION__
#define __SPHERE_PARTITION__

#include <PhysBAM_Tools/Arrays/ARRAY.h>
#include <PhysBAM_Tools/Grids_Uniform/GRID.h>
#include <PhysBAM_Tools/Utilities/NONCOPYABLE.h>
#include <PhysBAM_Geometry/Basic_Geometry/BOX.h>
#include <PhysBAM_Geometry/Basic_Geometry/SPHERE.h>
namespace PhysBAM{

template<class T>
class SPHERE_PARTITION:public NONCOPYABLE
{
    typedef VECTOR<T,3> TV;
public:
    ARRAY<SPHERE<TV> > spheres; // list of all the spheres
    GRID<TV> grid;
    BOX<TV> box; // box containing the spheres
    ARRAY<ARRAY<int>*,VECTOR<int,3> > voxel_sphere_list;
    
    SPHERE_PARTITION(const int number_input)
        :spheres(number_input)
    {
        Set_Up_Grid(1,1,1);
    }

    virtual ~SPHERE_PARTITION()
    {voxel_sphere_list.Delete_Pointers_And_Clean_Memory();}

//#####################################################################
    TV Normal(const TV& location,const int aggregate=0) const;
    bool Inside(const TV& location,const T thickness_over_two=0) const;
    bool Outside(const TV& location,const T thickness_over_two=0) const;
    bool Boundary(const TV& location,const T thickness_over_two) const;
    void Set_Up_Grid(const int m,const int n,const int mn); // use this to initialize!
    void Find_Voxel(const TV& location,int& i_left,int& j_bottom,int& ij_front) const;
//#####################################################################
};   
}
#endif

