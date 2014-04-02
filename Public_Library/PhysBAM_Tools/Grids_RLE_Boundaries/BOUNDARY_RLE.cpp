//#####################################################################
// Copyright 2002-2005, Ronald Fedkiw, Geoffrey Irving, Frank Losasso, Andrew Selle, Tamar Shinar.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class BOUNDARY_RLE
//#####################################################################
#ifndef COMPILE_WITHOUT_RLE_SUPPORT
#include <PhysBAM_Tools/Grids_RLE/RLE_GRID_ITERATOR_CELL_2D.h>
#include <PhysBAM_Tools/Grids_RLE/RLE_GRID_ITERATOR_CELL_3D.h>
#include <PhysBAM_Tools/Grids_RLE/RLE_GRID_SIMPLE_ITERATOR.h>
#include <PhysBAM_Tools/Grids_RLE/RLE_GRID_TRANSFER_ITERATOR.h>
#include <PhysBAM_Tools/Grids_RLE_Boundaries/BOUNDARY_RLE.h>
#include <PhysBAM_Tools/Matrices/SYMMETRIC_MATRIX_2X2.h>
#include <PhysBAM_Tools/Matrices/SYMMETRIC_MATRIX_3X3.h>
using namespace PhysBAM;
//#####################################################################
// Constructor
//#####################################################################
template<class T_GRID,class T2> BOUNDARY_RLE<T_GRID,T2>::
BOUNDARY_RLE()
{}
//#####################################################################
// Constructor
//#####################################################################
template<class T_GRID,class T2> BOUNDARY_RLE<T_GRID,T2>::
BOUNDARY_RLE(const TV_SIDES& constant_extrapolation)
{       
    Set_Constant_Extrapolation(constant_extrapolation);
}
//#####################################################################
// Destructor
//#####################################################################
template<class T_GRID,class T2> BOUNDARY_RLE<T_GRID,T2>::
~BOUNDARY_RLE()
{}
//#####################################################################
// Function Fill_Ghost_Cells
//#####################################################################
template<class T_GRID,class T2> void BOUNDARY_RLE<T_GRID,T2>::
Fill_Ghost_Cells(const T_GRID& grid,const typename BOUNDARY_RLE<T_GRID,T2>::T_ARRAYS_T2_BASE& u,typename BOUNDARY_RLE<T_GRID,T2>::T_ARRAYS_T2_BASE& u_ghost,const T dt,const T time,
    const int number_of_ghost_cells)
{
    ARRAY<T2>::Copy(u,u_ghost);
    Fill_Ghost_Cells_Helper<T2,typename BOUNDARY_RLE<T_GRID,T2>::T_ARRAYS_T2_BASE::ELEMENT,CELL_ITERATOR>(grid,u_ghost,fixed_boundary_value,lower_threshold,upper_threshold,number_of_ghost_cells);
}
//#####################################################################
// Function Fill_Ghost_Cells_Face
//#####################################################################
template<class T_GRID,class T2> struct Fill_Ghost_Cells_Face_Horizontal
{
    typedef typename T_GRID::VECTOR_T TV;typedef typename TV::SCALAR T;
    typedef typename GRID_ARRAYS_POLICY<T_GRID>::FACE_ARRAYS T_FACE_ARRAYS;typedef typename GRID_ARRAYS_POLICY<T_GRID>::ARRAYS_SCALAR T_ARRAYS_SCALAR;
    typedef typename T_GRID::VECTOR_HORIZONTAL_INT TV_HORIZONTAL_INT;typedef typename T_GRID::BOX_HORIZONTAL_INT T_BOX_HORIZONTAL_INT;
    typedef typename T_GRID::CELL_ITERATOR CELL_ITERATOR;
    typedef typename T_ARRAYS_SCALAR::template REBIND<T2>::TYPE T_ARRAYS_T2_BASE;
    typedef typename T_FACE_ARRAYS::template REBIND<T2>::TYPE T_FACE_ARRAYS_T2_BASE;
    typedef typename REBIND<typename BOUNDARY_POLICY<T_GRID>::BOUNDARY_SCALAR,T2>::TYPE T_BOUNDARY_T2;

    template<class T_FACE> static void
    Apply(const BOUNDARY_RLE<T_GRID,T2>& boundary,const T_GRID& grid,T_FACE_ARRAYS_T2_BASE& u_ghost,
        const typename T_GRID::VECTOR_T& fixed,const typename T_GRID::VECTOR_T& lower,
        const typename T_GRID::VECTOR_T& upper,const int number_of_ghost_cells)
    {
        int axis=T_FACE::Axis();
        boundary.template Fill_Ghost_Cells_Helper<typename T_GRID::SCALAR,typename T_FACE_ARRAYS_T2_BASE::ELEMENT,T_FACE>(grid,u_ghost,fixed[axis],lower[axis],upper[axis],number_of_ghost_cells);
    }
};
template<class T_GRID,class T2> void BOUNDARY_RLE<T_GRID,T2>::
Fill_Ghost_Cells_Face(const T_GRID& grid,const typename BOUNDARY_RLE<T_GRID,T2>::T_FACE_ARRAYS_T2_BASE& u,typename BOUNDARY_RLE<T_GRID,T2>::T_FACE_ARRAYS_T2_BASE& u_ghost,const T time,
    const int number_of_ghost_cells)
{
    ARRAY<typename BOUNDARY_RLE<T_GRID,T2>::T_FACE_ARRAYS_T2_BASE::ELEMENT>::Copy(u,u_ghost);
    //TV &fixed=Hack<TV,T2>::Convert(fixed_boundary_value),&lower=Hack<TV,T2>::Convert(lower_threshold),&upper=Hack<TV,T2>::Convert(upper_threshold);
    if(use_fixed_boundary||clamp_below||clamp_above) PHYSBAM_NOT_IMPLEMENTED();
    TV fixed,lower,upper;
    T_GRID::template Horizontal_Face_Loop<Fill_Ghost_Cells_Face_Horizontal<T_GRID,T2> >(*this,grid,u_ghost,fixed,lower,upper,number_of_ghost_cells);
    Fill_Ghost_Cells_Helper<T,typename BOUNDARY_RLE<T_GRID,T2>::T_FACE_ARRAYS_T2_BASE::ELEMENT,typename T_GRID::FACE_Y_ITERATOR>(grid,u_ghost,fixed[2],lower[2],upper[2],number_of_ghost_cells);
}
template<class T_GRID,class T2> template<class T3,class T4,class T_ITERATOR> void BOUNDARY_RLE<T_GRID,T2>::
Fill_Ghost_Cells_Helper(const T_GRID& grid,ARRAY_BASE<T4,ARRAY<T4> >& u_ghost,const T3& fixed_boundary_value,const T3& lower_threshold,const T3& upper_threshold,const int number_of_ghost_cells) const
{
    PHYSBAM_ASSERT(number_of_ghost_cells==grid.number_of_ghost_cells);
    ARRAY<T_BOX_HORIZONTAL_INT> regions;grid.Find_Ghost_Regions(regions,T_ITERATOR::Sentinels(),grid.number_of_ghost_cells);
    for(int r=1;r<=regions.m;r++)Fill_Single_Ghost_Region<T3,T4,T_ITERATOR>(grid,u_ghost,r,regions(r),fixed_boundary_value,lower_threshold,upper_threshold);
}
//#####################################################################
// Function Fill_Single_Ghost_Region
//#####################################################################
template<class T_GRID,class T2> template<class T3,class T4,class T_ITERATOR> void BOUNDARY_RLE<T_GRID,T2>::
Fill_Single_Ghost_Region(const T_GRID& grid,ARRAY_BASE<T4,ARRAY<T4> >& u_ghost,const int horizontal_side,const typename T_GRID::BOX_HORIZONTAL_INT& region,
    const T3& fixed_boundary_value,const T3& lower_threshold,const typename DISABLE_IF<IS_SAME<T3,T4>::value,T3>::TYPE& upper_threshold) const
{
    PHYSBAM_FATAL_ERROR();
}
template<class T_GRID,class T2> template<class T3,class T4,class T_ITERATOR> void BOUNDARY_RLE<T_GRID,T2>::
Fill_Single_Ghost_Region(const T_GRID& grid,ARRAY_BASE<T4,ARRAY<T4> >& u_ghost,const int horizontal_side,const typename T_GRID::BOX_HORIZONTAL_INT& region,
    const T3& fixed_boundary_value,const T3& lower_threshold,const typename ENABLE_IF<IS_SAME<T3,T4>::value,T3>::TYPE& upper_threshold) const
{
#if 0
    int horizontal_axis=(horizontal_side+1)/2;
    int boundary=horizontal_side&1?region.Maximum_Corner()[horizontal_axis]+1:region.Minimum_Corner()[horizontal_axis]-1;
    if(!use_fixed_boundary){
        for(RLE_GRID_SIMPLE_ITERATOR<T_GRID,T_ITERATOR> simple(grid,region);simple;simple++)u_ghost(simple.index)=T4();
        TV_HORIZONTAL_INT minimum=region.Minimum_Corner(),maximum=region.Maximum_Corner();
        TV_HORIZONTAL_INT boundary_min=minimum,boundary_max=maximum;boundary_min[horizontal_axis]=boundary;boundary_max[horizontal_axis]=boundary;
        T_BOX_HORIZONTAL_INT boundary_layer(boundary_min,boundary_max);
        for(int layer=minimum[horizontal_axis];layer<=maximum[horizontal_axis];layer++){
            TV_HORIZONTAL_INT ghost_min=minimum,ghost_max=maximum;ghost_min[horizontal_axis]=layer;ghost_max[horizontal_axis]=layer;
            T_BOX_HORIZONTAL_INT ghost_layer(ghost_min,ghost_max);
            for(RLE_GRID_TRANSFER_ITERATOR<T_GRID,T_ITERATOR> transfer(grid,boundary_layer,grid,ghost_layer);transfer;transfer++)
                Transfer(transfer,u_ghost,u_ghost);}}
    RLE_GRID_SIMPLE_ITERATOR<T_GRID,T_ITERATOR> simple(grid,region);
    if(use_fixed_boundary) for(;simple;simple++)u_ghost(simple.index)=fixed_boundary_value;
    else if(clamp_below&&clamp_above) for(;simple;simple++)u_ghost(simple.index)=clamp(u_ghost(simple.index),lower_threshold,upper_threshold);
    else if(clamp_below) for(;simple;simple++)u_ghost(simple.index)=clamp_min(u_ghost(simple.index),lower_threshold);
    else if(clamp_above) for(;simple;simple++)u_ghost(simple.index)=clamp_max(u_ghost(simple.index),upper_threshold);
#endif
}
//#####################################################################
#define INSTANTIATION_HELPER(T_GRID) \
    template class BOUNDARY_RLE<T_GRID,T_GRID::SCALAR>; \
    template class BOUNDARY_RLE<T_GRID,T_GRID::VECTOR_T>; \
    template class BOUNDARY_RLE<T_GRID,MATRIX_POLICY<T_GRID::VECTOR_T>::SYMMETRIC_MATRIX>;
INSTANTIATION_HELPER(RLE_GRID_2D<float>)
INSTANTIATION_HELPER(RLE_GRID_3D<float>)
#ifndef COMPILE_WITHOUT_DOUBLE_SUPPORT
INSTANTIATION_HELPER(RLE_GRID_2D<double>)
INSTANTIATION_HELPER(RLE_GRID_3D<double>)
#endif
#endif
