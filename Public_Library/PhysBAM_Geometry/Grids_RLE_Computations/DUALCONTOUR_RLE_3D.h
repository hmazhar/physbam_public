//#####################################################################
// Copyright 2005, Geoffrey Irving.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class DUALCONTOUR_RLE_3D
//##################################################################### 
#ifndef COMPILE_WITHOUT_RLE_SUPPORT
#ifndef __DUALCONTOUR_RLE_3D__
#define __DUALCONTOUR_RLE_3D__

#include <PhysBAM_Tools/Arrays/ARRAY.h>
#include <PhysBAM_Tools/Utilities/NONCOPYABLE.h>
#include <PhysBAM_Tools/Vectors/VECTOR.h>
namespace PhysBAM{

template<class T> class RLE_GRID_3D;
template<class T_GRID> class LEVELSET_RLE;
template<class T> class TRIANGULATED_SURFACE;

template<class T>
class DUALCONTOUR_RLE_3D:public NONCOPYABLE
{
    typedef VECTOR<T,3> TV;
private:
    LEVELSET_RLE<RLE_GRID_3D<T> >& levelset;
    const int number_of_ghost_cells;
    ARRAY<VECTOR<int,4> > quads;
    ARRAY<VECTOR<TV,2> > positions_and_normals;
public:

    DUALCONTOUR_RLE_3D(LEVELSET_RLE<RLE_GRID_3D<T> >& levelset_input,const int number_of_ghost_cells_input=3)
        :levelset(levelset_input),number_of_ghost_cells(number_of_ghost_cells_input)
    {}

    static TRIANGULATED_SURFACE<T>* Create_Triangulated_Surface_From_Levelset(LEVELSET_RLE<RLE_GRID_3D<T> >& levelset,const int number_of_ghost_cells=3)
    {DUALCONTOUR_RLE_3D<T> dualcontour(levelset,number_of_ghost_cells);dualcontour.Dualcontour();return dualcontour.Get_Triangulated_Surface();} 

//#####################################################################
    void Dualcontour();
    TRIANGULATED_SURFACE<T>* Get_Triangulated_Surface();
private:
    struct Generate_Quadrilateral_Indices{template<class T_FACE> static void Apply(const DUALCONTOUR_RLE_3D<T>& dualcontour,ARRAY<int>& face_quads,int& total_quads);};
//#####################################################################
};
}
#endif
#endif
