//#####################################################################
// Copyright 2009, Geoffrey Irving, Nipun Kwatra.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#ifndef COMPILE_WITHOUT_RLE_SUPPORT
#include <PhysBAM_Tools/Grids_RLE_Computations/VORTICITY_RLE.h>
#include <PhysBAM_Tools/Grids_RLE_Interpolation/LINEAR_INTERPOLATION_RLE_HELPER.h>
using namespace PhysBAM;
//#####################################################################
// Function Cell_Centered_Face_Velocity_Difference
//#####################################################################
template<class T_GRID,class TV> static inline typename T_GRID::SCALAR
Cell_Centered_Face_Velocity_Difference(const T_GRID& grid,const int axis,const int difference_axis,const ARRAY<TV>& V_cell,const int neighbors[T_GRID::number_of_neighbors_per_cell])
{   
    typedef typename T_GRID::SCALAR T;
    T one_over_two_dh=(T).5*grid.uniform_grid.one_over_dX[difference_axis];
    int c1=neighbors[2*difference_axis-2],c2=neighbors[2*difference_axis-1];
    return one_over_two_dh*(V_cell(c2)[axis]-V_cell(c1)[axis]);
}   
//#####################################################################
// Function Vorticity
//#####################################################################
template<class T> static inline VECTOR<T,1>
Vorticity_Helper(const RLE_GRID_2D<T>& grid,const ARRAY<VECTOR<T,2> >& V_cell,const int neighbors[4])
{
    return VECTOR<T,1>(Cell_Centered_Face_Velocity_Difference<RLE_GRID_2D<T> >(grid,2,1,V_cell,neighbors)-Cell_Centered_Face_Velocity_Difference<RLE_GRID_2D<T> >(grid,1,2,V_cell,neighbors));
}   
template<class T> static inline VECTOR<T,3>
Vorticity_Helper(const RLE_GRID_3D<T>& grid,const ARRAY<VECTOR<T,3> >& V_cell,const int neighbors[6])
{   
    return VECTOR<T,3>(Cell_Centered_Face_Velocity_Difference<RLE_GRID_3D<T> >(grid,3,2,V_cell,neighbors)-Cell_Centered_Face_Velocity_Difference<RLE_GRID_3D<T> >(grid,2,3,V_cell,neighbors),
                        Cell_Centered_Face_Velocity_Difference<RLE_GRID_3D<T> >(grid,1,3,V_cell,neighbors)-Cell_Centered_Face_Velocity_Difference<RLE_GRID_3D<T> >(grid,3,1,V_cell,neighbors),
                        Cell_Centered_Face_Velocity_Difference<RLE_GRID_3D<T> >(grid,2,1,V_cell,neighbors)-Cell_Centered_Face_Velocity_Difference<RLE_GRID_3D<T> >(grid,1,2,V_cell,neighbors));
}
template<class TV> void VORTICITY_RLE<TV>::
Vorticity(const RLE_GRID& grid,const ARRAY<T>& V,ARRAY<T_SPIN>& vorticity)
{
    ARRAY<TV> V_cell(grid.number_of_cells,false);LINEAR_INTERPOLATION_RLE_HELPER<RLE_GRID>::Interpolate_From_Faces_To_Short_Cells(grid,V,V_cell);
    const ARRAY<VECTOR<int,RLE_GRID::number_of_neighbors_per_cell> >& short_cell_neighbors=grid.Short_Cell_Neighbors();
    for(int c=1;c<=grid.number_of_cells;c++){
        int neighbors[RLE_GRID::number_of_neighbors_per_cell];bool has_all_neighbors=true;
        for(int n=1;n<=RLE_GRID::number_of_neighbors_per_cell;n++){
            neighbors[n-1]=short_cell_neighbors(c)(n);
            if(!neighbors[n-1]){has_all_neighbors=false;break;}}
        vorticity(c)=has_all_neighbors?Vorticity_Helper(grid,V_cell,neighbors):T_SPIN();}
}
template class VORTICITY_RLE<VECTOR<float,2> >;
template class VORTICITY_RLE<VECTOR<float,3> >;
#ifndef COMPILE_WITHOUT_DOUBLE_SUPPORT
template class VORTICITY_RLE<VECTOR<double,2> >;
template class VORTICITY_RLE<VECTOR<double,3> >;
#endif
#endif
