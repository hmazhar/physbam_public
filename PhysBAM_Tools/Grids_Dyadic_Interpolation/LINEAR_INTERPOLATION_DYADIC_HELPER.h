//#####################################################################
// Copyright 2005, Frank Losasso.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class LINEAR_INTERPOLATION_DYADIC_HELPER
//#####################################################################
#if !COMPILE_WITHOUT_DYADIC_SUPPORT || COMPILE_WITH_BINTREE_SUPPORT
#ifndef __LINEAR_INTERPOLATION_DYADIC_HELPER__
#define __LINEAR_INTERPOLATION_DYADIC_HELPER__

#include <PhysBAM_Tools/Arrays/ARRAYS_FORWARD.h>
namespace PhysBAM{

template<class T_GRID>
class LINEAR_INTERPOLATION_DYADIC_HELPER
{
    typedef typename T_GRID::VECTOR_T TV;typedef typename TV::SCALAR T;typedef typename T_GRID::CELL T_CELL;
private:
    LINEAR_INTERPOLATION_DYADIC_HELPER(); // disallow construction
public:
    template<class T2> static void Interpolate_From_Cells_To_Faces(const T_GRID& grid,const ARRAY<T2>& cell_based,ARRAY<T2>& face_based)
    {ARRAY<T2> node_based(grid.number_of_nodes,false);Interpolate_From_Cells_To_Nodes(grid,cell_based,node_based);
    Interpolate_From_Cells_To_Faces(grid,cell_based,node_based,face_based);}

//#####################################################################
    template<class T2> static void Interpolate_From_Nodes_To_Cells(const T_GRID& grid,const ARRAY<T2>& node_based,ARRAY<T2>& cell_based);
    template<class T2> static void Interpolate_From_Cells_To_Faces(const T_GRID& grid,const ARRAY<T2>& cell_based,const ARRAY<T2>& node_based,ARRAY<T2>& face_based);
    template<class T2> static void Interpolate_From_Cells_To_Nodes(const T_GRID& grid,const ARRAY<T2>& cell_based,ARRAY<T2>& node_based);
    static void Interpolate_From_Faces_To_Nodes(const T_GRID& grid,const ARRAY<T>& face_based,ARRAY<TV>& node_based);    
    static void Interpolate_From_Faces_To_Nodes(const T_GRID& grid,const ARRAY<T>& face_based,ARRAY<T>& node_based);    
//#####################################################################
};
}
#endif
#endif
