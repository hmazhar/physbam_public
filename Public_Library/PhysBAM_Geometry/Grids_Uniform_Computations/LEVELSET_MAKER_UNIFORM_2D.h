//#####################################################################
// Copyright 2010.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class LEVELSET_MAKER_UNIFORM_2D
//##################################################################### 
#ifndef __LEVELSET_MAKER_UNIFORM_2D__
#define __LEVELSET_MAKER_UNIFORM_2D__

#include <PhysBAM_Tools/Grids_Uniform_Arrays/ARRAYS_ND.h>
namespace PhysBAM{

template<class TV> class GRID;
template<class T> class SEGMENTED_CURVE_2D;

template<class T>
class LEVELSET_MAKER_UNIFORM_2D
{
    typedef VECTOR<int,2> TV_INT;typedef VECTOR<T,2> TV;
public:
//#####################################################################
    static void Compute_Level_Set(SEGMENTED_CURVE_2D<T>& curve,GRID<TV>& grid,int ghost_cells,ARRAY<T,TV_INT>& phi);
private:
    static void Compute_Level_Set_Helper(const TV_INT& index,T next,ARRAY<TV_INT>& next_todo,ARRAY<T,TV_INT>& phi);
//#####################################################################
};
}
#endif
