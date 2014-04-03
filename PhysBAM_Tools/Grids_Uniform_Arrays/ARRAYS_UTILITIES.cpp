//#####################################################################
// Copyright 2006-2007, Nipun Kwatra, Frank Losasso, Jerry Talton.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#include <PhysBAM_Tools/Grids_Uniform/UNIFORM_GRID_ITERATOR_CELL.h>
#include <PhysBAM_Tools/Grids_Uniform/UNIFORM_GRID_ITERATOR_FACE.h>
#include <PhysBAM_Tools/Grids_Uniform/UNIFORM_GRID_ITERATOR_NODE.h>
#include <PhysBAM_Tools/Grids_Uniform_Arrays/ARRAYS_UTILITIES.h>
#include <PhysBAM_Tools/Grids_Uniform_Arrays/FACE_ARRAYS.h>
using namespace PhysBAM;
//#####################################################################
// Function Make_Ghost_Mask_From_Active_Mask
//#####################################################################
template<class T_GRID,class T2> void ARRAYS_UTILITIES<T_GRID,T2>::
Make_Ghost_Mask_From_Active_Mask(const T_GRID& grid,const T_ARRAYS_BOOL& input_mask,T_ARRAYS_BOOL& output_mask,const int stencil_width,const int number_of_ghost_cells)
{
    // TODO: this could be rewritten with linear complexity independent of the size of stencil_width by handling one dimension at a time
    assert(!grid.Is_MAC_Grid());T_ARRAYS_INT temp_mask(grid.Domain_Indices(number_of_ghost_cells+1));
    for(NODE_ITERATOR iterator(grid,number_of_ghost_cells);iterator.Valid();iterator.Next()){TV_INT node_index=iterator.Node_Index();temp_mask(node_index)=input_mask(node_index);}
    for(int i=1;i<=stencil_width;i++)for(NODE_ITERATOR iterator(grid,number_of_ghost_cells);iterator.Valid();iterator.Next()){TV_INT node_index=iterator.Node_Index();
        if(temp_mask(node_index)) continue;
        for(int neighbor=1;neighbor<=grid.number_of_neighbors_per_node;neighbor++){TV_INT neighbor_index=iterator.Node_Neighbor(neighbor);
            if(temp_mask(neighbor_index)==i){temp_mask(node_index)=i+1;break;}}}
    for(NODE_ITERATOR iterator(grid,number_of_ghost_cells);iterator.Valid();iterator.Next()){TV_INT node_index=iterator.Node_Index();
        output_mask(node_index)=temp_mask(node_index)>1;}
}
//#####################################################################
// Function Compute_Face_Data_From_Cell_Data
//#####################################################################
template<class T_GRID,class T2> void ARRAYS_UTILITIES<T_GRID,T2>::
Compute_Face_Data_From_Cell_Data(const T_GRID& face_grid,T_FACE_ARRAYS_T2& face_array,const T_ARRAYS_DIMENSION_T2& cell_array,const int number_of_ghost_cells)
{
    for(FACE_ITERATOR iterator(face_grid,number_of_ghost_cells);iterator.Valid();iterator.Next()){int axis=iterator.Axis();
        TV_INT face_index=iterator.Face_Index(),first_cell_index=iterator.First_Cell_Index(),second_cell_index=iterator.Second_Cell_Index();
        face_array.Component(axis)(face_index)=(cell_array(first_cell_index)+cell_array(second_cell_index))*(T).5;}
}
//#####################################################################
// Function Compute_Gradient_At_Faces_From_Cell_Data
//#####################################################################
template<class T_GRID,class T2> void ARRAYS_UTILITIES<T_GRID,T2>::
Compute_Gradient_At_Faces_From_Cell_Data(const T_GRID& face_grid,T_FACE_ARRAYS_T2& grad_face_array,const T_ARRAYS_DIMENSION_T2& cell_array,const int number_of_ghost_cells)
{
    TV one_over_dx=face_grid.one_over_dX;
    for(FACE_ITERATOR iterator(face_grid,number_of_ghost_cells);iterator.Valid();iterator.Next()){int axis=iterator.Axis();
        TV_INT face_index=iterator.Face_Index(),first_cell_index=iterator.First_Cell_Index(),second_cell_index=iterator.Second_Cell_Index();
        grad_face_array.Component(axis)(face_index)=one_over_dx[iterator.Axis()]*(cell_array(second_cell_index)-cell_array(first_cell_index));}
}
//#####################################################################
// Function Compute_Gradient_At_Cells_From_Face_Data
//#####################################################################
template<class T_GRID,class T2> void ARRAYS_UTILITIES<T_GRID,T2>::
Compute_Gradient_At_Cells_From_Face_Data(const T_GRID& face_grid,T_ARRAYS_DIMENSION_VECTOR_T2& grad_cell_array,const T_FACE_ARRAYS_T2& face_array,const int number_of_ghost_cells)
{
    TV one_over_dx=face_grid.one_over_dX;
    for(CELL_ITERATOR iterator(face_grid,number_of_ghost_cells);iterator.Valid();iterator.Next()){TV_INT cell_index=iterator.Cell_Index();
        for(int axis=1;axis<=T_GRID::dimension;axis++){TV_INT first_face_index=iterator.First_Face_Index(axis),second_face_index=iterator.Second_Face_Index(axis);
            grad_cell_array(cell_index)[axis]=one_over_dx[axis]*(face_array.Component(axis)(second_face_index)-face_array.Component(axis)(first_face_index));}}
}
//#####################################################################
// Function Compute_Divergence_At_Cells_From_Face_Data
//#####################################################################
template<class T_GRID,class T2> void ARRAYS_UTILITIES<T_GRID,T2>::
Compute_Divergence_At_Cells_From_Face_Data(const T_GRID& face_grid,T_ARRAYS_DIMENSION_T2& div_cell_array,const T_FACE_ARRAYS_T2& face_array,const int number_of_ghost_cells)
{
    TV one_over_dx=face_grid.one_over_dX;
    for(CELL_ITERATOR iterator(face_grid,number_of_ghost_cells);iterator.Valid();iterator.Next()){TV_INT cell_index=iterator.Cell_Index();
        div_cell_array(cell_index)=T2();
        for(int axis=1;axis<=T_GRID::dimension;axis++){TV_INT first_face_index=iterator.First_Face_Index(axis),second_face_index=iterator.Second_Face_Index(axis);
            div_cell_array(cell_index)+=one_over_dx[axis]*(face_array.Component(axis)(second_face_index)-face_array.Component(axis)(first_face_index));}}
}
//#####################################################################
#define INSTANTIATION_HELPER(T) \
    template class ARRAYS_UTILITIES<GRID<VECTOR<T,1> >,T>; \
    template class ARRAYS_UTILITIES<GRID<VECTOR<T,2> >,T>; \
    template class ARRAYS_UTILITIES<GRID<VECTOR<T,3> >,T>; \
    template class ARRAYS_UTILITIES<GRID<VECTOR<T,1> >,VECTOR<T,1> >; \
    template class ARRAYS_UTILITIES<GRID<VECTOR<T,1> >,VECTOR<T,2> >; \
    template class ARRAYS_UTILITIES<GRID<VECTOR<T,1> >,VECTOR<T,3> >; \
    template class ARRAYS_UTILITIES<GRID<VECTOR<T,1> >,VECTOR<T,4> >; \
    template class ARRAYS_UTILITIES<GRID<VECTOR<T,1> >,VECTOR<T,5> >; \
    template class ARRAYS_UTILITIES<GRID<VECTOR<T,1> >,VECTOR<T,6> >; \
    template class ARRAYS_UTILITIES<GRID<VECTOR<T,2> >,VECTOR<T,1> >; \
    template class ARRAYS_UTILITIES<GRID<VECTOR<T,2> >,VECTOR<T,2> >; \
    template class ARRAYS_UTILITIES<GRID<VECTOR<T,2> >,VECTOR<T,3> >; \
    template class ARRAYS_UTILITIES<GRID<VECTOR<T,2> >,VECTOR<T,4> >; \
    template class ARRAYS_UTILITIES<GRID<VECTOR<T,2> >,VECTOR<T,5> >; \
    template class ARRAYS_UTILITIES<GRID<VECTOR<T,2> >,VECTOR<T,6> >; \
    template class ARRAYS_UTILITIES<GRID<VECTOR<T,3> >,VECTOR<T,1> >; \
    template class ARRAYS_UTILITIES<GRID<VECTOR<T,3> >,VECTOR<T,2> >; \
    template class ARRAYS_UTILITIES<GRID<VECTOR<T,3> >,VECTOR<T,3> >; \
    template class ARRAYS_UTILITIES<GRID<VECTOR<T,3> >,VECTOR<T,4> >; \
    template class ARRAYS_UTILITIES<GRID<VECTOR<T,3> >,VECTOR<T,5> >; \
    template class ARRAYS_UTILITIES<GRID<VECTOR<T,3> >,VECTOR<T,6> >;
INSTANTIATION_HELPER(float)
#ifndef COMPILE_WITHOUT_DOUBLE_SUPPORT
INSTANTIATION_HELPER(double)
#endif
