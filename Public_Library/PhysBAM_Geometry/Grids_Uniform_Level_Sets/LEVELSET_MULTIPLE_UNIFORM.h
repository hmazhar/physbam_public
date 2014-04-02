//#####################################################################
// Copyright 2005, Ron Fedkiw, Eran Guendelman, Frank Losasso, Tamar Shinar.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class LEVELSET_MULTIPLE_UNIFORM  
//##################################################################### 
#ifndef __LEVELSET_MULTIPLE_UNIFORM__
#define __LEVELSET_MULTIPLE_UNIFORM__ 

#include <PhysBAM_Geometry/Grids_Uniform_Level_Sets/LEVELSET_MULTIPLE.h>
namespace PhysBAM{

template<class T_GRID>
class LEVELSET_MULTIPLE_UNIFORM:public LEVELSET_MULTIPLE<T_GRID>
{
    typedef typename T_GRID::SCALAR T;
    typedef typename GRID_ARRAYS_POLICY<T_GRID>::ARRAYS_SCALAR T_ARRAYS_SCALAR;typedef typename T_GRID::CELL_ITERATOR CELL_ITERATOR;typedef typename LEVELSET_POLICY<T_GRID>::LEVELSET T_LEVELSET;
    typedef typename GRID_ARRAYS_POLICY<T_GRID>::FACE_ARRAYS T_FACE_ARRAYS_SCALAR;typedef typename T_FACE_ARRAYS_SCALAR::template REBIND<bool>::TYPE T_FACE_ARRAYS_BOOL;
    typedef typename T_GRID::FACE_ITERATOR FACE_ITERATOR;typedef typename GRID_ARRAYS_POLICY<T_GRID>::FLOOD_FILL T_FLOOD_FILL;typedef typename T_ARRAYS_SCALAR::template REBIND<int>::TYPE T_ARRAYS_INT;
public:
    using LEVELSET_MULTIPLE<T_GRID>::grid;using LEVELSET_MULTIPLE<T_GRID>::phis;using LEVELSET_MULTIPLE<T_GRID>::Set_Collision_Body_List;

    LEVELSET_MULTIPLE_UNIFORM(T_GRID& grid_input,ARRAY<T_ARRAYS_SCALAR>& phis_input,const bool use_external_levelsets_input=false);
    ~LEVELSET_MULTIPLE_UNIFORM();

//#####################################################################
    void Project_Levelset(const int number_of_ghost_cells=0);
    void Get_Single_Levelset(const ARRAY<bool>& positive_regions,T_LEVELSET& levelset,const bool flood_fill_for_bubbles);
//#####################################################################
};   
}
#endif
