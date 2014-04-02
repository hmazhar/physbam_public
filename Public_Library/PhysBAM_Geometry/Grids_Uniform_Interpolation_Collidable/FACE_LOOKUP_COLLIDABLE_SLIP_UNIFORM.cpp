//#####################################################################
// Copyright 2009, Elliot English, Avi Robinson-Mosher.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#include <PhysBAM_Tools/Grids_Uniform/UNIFORM_GRID_ITERATOR_CELL.h>
#include <PhysBAM_Tools/Grids_Uniform_Arrays/FACE_ARRAYS_BINARY_UNIFORM.h>
#include <PhysBAM_Tools/Grids_Uniform_Interpolation/LINEAR_INTERPOLATION_MAC_2D_HELPER_DEFINITIONS.h>
#include <PhysBAM_Tools/Grids_Uniform_Interpolation/LINEAR_INTERPOLATION_MAC_3D_HELPER_DEFINITIONS.h>
#include <PhysBAM_Geometry/Grids_Uniform_Collisions/GRID_BASED_COLLISION_GEOMETRY_UNIFORM.h>
#include <PhysBAM_Geometry/Grids_Uniform_Interpolation_Collidable/FACE_LOOKUP_COLLIDABLE_SLIP_UNIFORM.h>
using namespace PhysBAM;
//#####################################################################
// Constructor
//#####################################################################
template<class T_GRID,class T_NESTED_LOOKUP> FACE_LOOKUP_COLLIDABLE_SLIP_UNIFORM<T_GRID,T_NESTED_LOOKUP>::
FACE_LOOKUP_COLLIDABLE_SLIP_UNIFORM(const T_NESTED_LOOKUP& nested_face_lookup_input,const T_GRID_BASED_COLLISION_GEOMETRY& body_list_input,const FACE_ARRAYS_BOOL* valid_value_mask_input)
    :nested_face_lookup(nested_face_lookup_input),body_list(body_list_input),valid_value_mask(valid_value_mask_input)
{
}
//#####################################################################
// Constructor
//#####################################################################
template<class T_GRID,class T_NESTED_LOOKUP> FACE_LOOKUP_COLLIDABLE_SLIP_UNIFORM<T_GRID,T_NESTED_LOOKUP>::LOOKUP::
LOOKUP(const FACE_LOOKUP_COLLIDABLE_SLIP_UNIFORM<T_GRID,T_NESTED_LOOKUP>& face_lookup_input,const typename T_NESTED_LOOKUP::LOOKUP& nested_lookup_input)
    :face_lookup(face_lookup_input),nested_lookup(nested_lookup_input),reference_point_set(false),reference_point_inside(false),both_cells_visible(true),using_reference_point(true)
{
}
//#####################################################################
// Constructor
//#####################################################################
template<class T_GRID,class T_NESTED_LOOKUP> FACE_LOOKUP_COLLIDABLE_SLIP_UNIFORM<T_GRID,T_NESTED_LOOKUP>::LOOKUP::
LOOKUP(const FACE_LOOKUP_COLLIDABLE_SLIP_UNIFORM<T_GRID,T_NESTED_LOOKUP>& face_lookup_input,const typename T_NESTED_LOOKUP::LOOKUP& nested_lookup_input,const TV_INT& cell_input,
    const bool both_cells_visible_input)
    :face_lookup(face_lookup_input),nested_lookup(nested_lookup_input),reference_point_set(false),reference_point_inside(false),cell(cell_input),both_cells_visible(both_cells_visible_input),using_reference_point(false)
{
}
//#####################################################################
// Function Set_Reference_Point
//#####################################################################
template<class T_GRID,class T_NESTED_LOOKUP> void FACE_LOOKUP_COLLIDABLE_SLIP_UNIFORM<T_GRID,T_NESTED_LOOKUP>::LOOKUP::
Set_Reference_Point(const TV& reference_point_input) const
{
    reference_point=reference_point_input;
    assert(!reference_point_set && !reference_point_inside);reference_point_set=true;
    COLLISION_GEOMETRY_ID body_id(0);
    int aggregate_id=0;
    found_valid_point=false;
    // TODO: Inside_Any_Body doesn't do quite the right thing for this body - so this is wrong, until we make the rigid body behave right simplicially
    if(face_lookup.body_list.Inside_Any_Simplex_Of_Any_Body(reference_point,body_id,aggregate_id)){
        reference_point_inside=true;
        object_velocity=face_lookup.body_list.Object_Velocity(body_id,aggregate_id,reference_point);}
    else{
        found_valid_point=true; // TODO: check this is correctly placed!
        BLOCK_UNIFORM<T_GRID> block(face_lookup.body_list.grid,reference_point,3);
        block.All_Cell_Indices(cells_in_block);}
    BLOCK_UNIFORM<T_GRID> block(face_lookup.body_list.grid,reference_point,3);
    TV_INT reference_point_cell=block.block_index-TV_INT::All_Ones_Vector();
    nested_lookup.Set_Reference_Point((TV)reference_point_cell);
}
//#####################################################################
// Function operator()
//#####################################################################
template<class T_GRID,class T_NESTED_LOOKUP> typename T_GRID::VECTOR_T::SCALAR FACE_LOOKUP_COLLIDABLE_SLIP_UNIFORM<T_GRID,T_NESTED_LOOKUP>::LOOKUP::
operator()(const int axis,const TV_INT& face) const
{
    if(reference_point_set && using_reference_point){
        TV location=face_lookup.body_list.grid.Face(axis,face);
        int side;
        if(location[axis]>=reference_point[axis])
            side=1;
        else
            side=2;
        T face_velocity=0;
        if(face_lookup.body_list.Face_Velocity(side,axis,face,cells_in_block,T_GRID::number_of_cells_per_block,reference_point,face_velocity)){
            found_valid_point=true;return face_velocity;}
        return nested_lookup(axis,face);}
    else{
        if(both_cells_visible)
            return nested_lookup(axis,face);
        else{ // "reflect" values to get correct averaging to face - is this correct for more general interpolation?
            if(face(axis)==cell(axis)) // use "right side" of face
                return nested_lookup(axis,cell);
            else if(face(axis)==cell(axis)+1) // use "left side" of face
                return nested_lookup(axis,cell+TV_INT::Axis_Vector(axis));}}
    return 0;
}
template class FACE_LOOKUP_COLLIDABLE_SLIP_UNIFORM<GRID<VECTOR<float,1> >,FACE_LOOKUP_UNIFORM<GRID<VECTOR<float,1> > > >;
template class FACE_LOOKUP_COLLIDABLE_SLIP_UNIFORM<GRID<VECTOR<float,2> >,FACE_LOOKUP_UNIFORM<GRID<VECTOR<float,2> > > >;
template class FACE_LOOKUP_COLLIDABLE_SLIP_UNIFORM<GRID<VECTOR<float,3> >,FACE_LOOKUP_UNIFORM<GRID<VECTOR<float,3> > > >;
template VECTOR<float,2> LINEAR_INTERPOLATION_MAC_2D_HELPER<GRID<VECTOR<float,2> > >::Extrema_Face_X_Transformed<BLOCK_UNIFORM<GRID<VECTOR<float,2> > >,FACE_LOOKUP_COLLIDABLE_SLIP_UNIFORM<GRID<VECTOR<float,2> >,FACE_LOOKUP_UNIFORM<GRID<VECTOR<float,2> > > >::LOOKUP>(BLOCK_UNIFORM<GRID<VECTOR<float,2> > > const&,FACE_LOOKUP_COLLIDABLE_SLIP_UNIFORM<GRID<VECTOR<float,2> >,FACE_LOOKUP_UNIFORM<GRID<VECTOR<float,2> > > >::LOOKUP const&,FACE_LOOKUP_COLLIDABLE_SLIP_UNIFORM<GRID<VECTOR<float,2> >,FACE_LOOKUP_UNIFORM<GRID<VECTOR<float,2> > > >::LOOKUP const&,VECTOR<float,2> const&);
template VECTOR<float,2> LINEAR_INTERPOLATION_MAC_2D_HELPER<GRID<VECTOR<float,2> > >::Extrema_Face_Y_Transformed<BLOCK_UNIFORM<GRID<VECTOR<float,2> > >,FACE_LOOKUP_COLLIDABLE_SLIP_UNIFORM<GRID<VECTOR<float,2> >,FACE_LOOKUP_UNIFORM<GRID<VECTOR<float,2> > > >::LOOKUP>(BLOCK_UNIFORM<GRID<VECTOR<float,2> > > const&,FACE_LOOKUP_COLLIDABLE_SLIP_UNIFORM<GRID<VECTOR<float,2> >,FACE_LOOKUP_UNIFORM<GRID<VECTOR<float,2> > > >::LOOKUP const&,FACE_LOOKUP_COLLIDABLE_SLIP_UNIFORM<GRID<VECTOR<float,2> >,FACE_LOOKUP_UNIFORM<GRID<VECTOR<float,2> > > >::LOOKUP const&,VECTOR<float,2> const&);
template VECTOR<float,2> LINEAR_INTERPOLATION_MAC_3D_HELPER<GRID<VECTOR<float,3> > >::Extrema_Face_X_Transformed<BLOCK_UNIFORM<GRID<VECTOR<float,3> > >,FACE_LOOKUP_COLLIDABLE_SLIP_UNIFORM<GRID<VECTOR<float,3> >,FACE_LOOKUP_UNIFORM<GRID<VECTOR<float,3> > > >::LOOKUP>(BLOCK_UNIFORM<GRID<VECTOR<float,3> > > const&,FACE_LOOKUP_COLLIDABLE_SLIP_UNIFORM<GRID<VECTOR<float,3> >,FACE_LOOKUP_UNIFORM<GRID<VECTOR<float,3> > > >::LOOKUP const&,FACE_LOOKUP_COLLIDABLE_SLIP_UNIFORM<GRID<VECTOR<float,3> >,FACE_LOOKUP_UNIFORM<GRID<VECTOR<float,3> > > >::LOOKUP const&,VECTOR<float,3> const&);
template VECTOR<float,2> LINEAR_INTERPOLATION_MAC_3D_HELPER<GRID<VECTOR<float,3> > >::Extrema_Face_Y_Transformed<BLOCK_UNIFORM<GRID<VECTOR<float,3> > >,FACE_LOOKUP_COLLIDABLE_SLIP_UNIFORM<GRID<VECTOR<float,3> >,FACE_LOOKUP_UNIFORM<GRID<VECTOR<float,3> > > >::LOOKUP>(BLOCK_UNIFORM<GRID<VECTOR<float,3> > > const&,FACE_LOOKUP_COLLIDABLE_SLIP_UNIFORM<GRID<VECTOR<float,3> >,FACE_LOOKUP_UNIFORM<GRID<VECTOR<float,3> > > >::LOOKUP const&,FACE_LOOKUP_COLLIDABLE_SLIP_UNIFORM<GRID<VECTOR<float,3> >,FACE_LOOKUP_UNIFORM<GRID<VECTOR<float,3> > > >::LOOKUP const&,VECTOR<float,3> const&);
template VECTOR<float,2> LINEAR_INTERPOLATION_MAC_3D_HELPER<GRID<VECTOR<float,3> > >::Extrema_Face_Z_Transformed<BLOCK_UNIFORM<GRID<VECTOR<float,3> > >,FACE_LOOKUP_COLLIDABLE_SLIP_UNIFORM<GRID<VECTOR<float,3> >,FACE_LOOKUP_UNIFORM<GRID<VECTOR<float,3> > > >::LOOKUP>(BLOCK_UNIFORM<GRID<VECTOR<float,3> > > const&,FACE_LOOKUP_COLLIDABLE_SLIP_UNIFORM<GRID<VECTOR<float,3> >,FACE_LOOKUP_UNIFORM<GRID<VECTOR<float,3> > > >::LOOKUP const&,FACE_LOOKUP_COLLIDABLE_SLIP_UNIFORM<GRID<VECTOR<float,3> >,FACE_LOOKUP_UNIFORM<GRID<VECTOR<float,3> > > >::LOOKUP const&,VECTOR<float,3> const&);
template float LINEAR_INTERPOLATION_MAC_2D_HELPER<GRID<VECTOR<float,2> > >::Interpolate_Face_X_Transformed<BLOCK_UNIFORM<GRID<VECTOR<float,2> > >,FACE_LOOKUP_COLLIDABLE_SLIP_UNIFORM<GRID<VECTOR<float,2> >,FACE_LOOKUP_UNIFORM<GRID<VECTOR<float,2> > > >::LOOKUP>(BLOCK_UNIFORM<GRID<VECTOR<float,2> > > const&,FACE_LOOKUP_COLLIDABLE_SLIP_UNIFORM<GRID<VECTOR<float,2> >,FACE_LOOKUP_UNIFORM<GRID<VECTOR<float,2> > > >::LOOKUP const&,VECTOR<float,2> const&);
template float LINEAR_INTERPOLATION_MAC_2D_HELPER<GRID<VECTOR<float,2> > >::Interpolate_Face_Y_Transformed<BLOCK_UNIFORM<GRID<VECTOR<float,2> > >,FACE_LOOKUP_COLLIDABLE_SLIP_UNIFORM<GRID<VECTOR<float,2> >,FACE_LOOKUP_UNIFORM<GRID<VECTOR<float,2> > > >::LOOKUP>(BLOCK_UNIFORM<GRID<VECTOR<float,2> > > const&,FACE_LOOKUP_COLLIDABLE_SLIP_UNIFORM<GRID<VECTOR<float,2> >,FACE_LOOKUP_UNIFORM<GRID<VECTOR<float,2> > > >::LOOKUP const&,VECTOR<float,2> const&);
template float LINEAR_INTERPOLATION_MAC_3D_HELPER<GRID<VECTOR<float,3> > >::Interpolate_Face_X_Transformed<BLOCK_UNIFORM<GRID<VECTOR<float,3> > >,FACE_LOOKUP_COLLIDABLE_SLIP_UNIFORM<GRID<VECTOR<float,3> >,FACE_LOOKUP_UNIFORM<GRID<VECTOR<float,3> > > >::LOOKUP>(BLOCK_UNIFORM<GRID<VECTOR<float,3> > > const&,FACE_LOOKUP_COLLIDABLE_SLIP_UNIFORM<GRID<VECTOR<float,3> >,FACE_LOOKUP_UNIFORM<GRID<VECTOR<float,3> > > >::LOOKUP const&,VECTOR<float,3> const&);
template float LINEAR_INTERPOLATION_MAC_3D_HELPER<GRID<VECTOR<float,3> > >::Interpolate_Face_Y_Transformed<BLOCK_UNIFORM<GRID<VECTOR<float,3> > >,FACE_LOOKUP_COLLIDABLE_SLIP_UNIFORM<GRID<VECTOR<float,3> >,FACE_LOOKUP_UNIFORM<GRID<VECTOR<float,3> > > >::LOOKUP>(BLOCK_UNIFORM<GRID<VECTOR<float,3> > > const&,FACE_LOOKUP_COLLIDABLE_SLIP_UNIFORM<GRID<VECTOR<float,3> >,FACE_LOOKUP_UNIFORM<GRID<VECTOR<float,3> > > >::LOOKUP const&,VECTOR<float,3> const&);
template float LINEAR_INTERPOLATION_MAC_3D_HELPER<GRID<VECTOR<float,3> > >::Interpolate_Face_Z_Transformed<BLOCK_UNIFORM<GRID<VECTOR<float,3> > >,FACE_LOOKUP_COLLIDABLE_SLIP_UNIFORM<GRID<VECTOR<float,3> >,FACE_LOOKUP_UNIFORM<GRID<VECTOR<float,3> > > >::LOOKUP>(BLOCK_UNIFORM<GRID<VECTOR<float,3> > > const&,FACE_LOOKUP_COLLIDABLE_SLIP_UNIFORM<GRID<VECTOR<float,3> >,FACE_LOOKUP_UNIFORM<GRID<VECTOR<float,3> > > >::LOOKUP const&,VECTOR<float,3> const&);
template ARRAY<PAIR<FACE_INDEX<2>,float> > LINEAR_INTERPOLATION_MAC_2D_HELPER<GRID<VECTOR<float,2> > >::Interpolate_Face_X_Transformed_Weights<BLOCK_UNIFORM<GRID<VECTOR<float,2> > >,FACE_LOOKUP_COLLIDABLE_SLIP_UNIFORM<GRID<VECTOR<float,2> >,FACE_LOOKUP_UNIFORM<GRID<VECTOR<float,2> > > >::LOOKUP>(BLOCK_UNIFORM<GRID<VECTOR<float,2> > > const&,FACE_LOOKUP_COLLIDABLE_SLIP_UNIFORM<GRID<VECTOR<float,2> >,FACE_LOOKUP_UNIFORM<GRID<VECTOR<float,2> > > >::LOOKUP const&,VECTOR<float,2> const&);
template ARRAY<PAIR<FACE_INDEX<2>,float> > LINEAR_INTERPOLATION_MAC_2D_HELPER<GRID<VECTOR<float,2> > >::Interpolate_Face_Y_Transformed_Weights<BLOCK_UNIFORM<GRID<VECTOR<float,2> > >,FACE_LOOKUP_COLLIDABLE_SLIP_UNIFORM<GRID<VECTOR<float,2> >,FACE_LOOKUP_UNIFORM<GRID<VECTOR<float,2> > > >::LOOKUP>(BLOCK_UNIFORM<GRID<VECTOR<float,2> > > const&,FACE_LOOKUP_COLLIDABLE_SLIP_UNIFORM<GRID<VECTOR<float,2> >,FACE_LOOKUP_UNIFORM<GRID<VECTOR<float,2> > > >::LOOKUP const&,VECTOR<float,2> const&);
template ARRAY<PAIR<FACE_INDEX<3>,float> > LINEAR_INTERPOLATION_MAC_3D_HELPER<GRID<VECTOR<float,3> > >::Interpolate_Face_X_Transformed_Weights<BLOCK_UNIFORM<GRID<VECTOR<float,3> > >,FACE_LOOKUP_COLLIDABLE_SLIP_UNIFORM<GRID<VECTOR<float,3> >,FACE_LOOKUP_UNIFORM<GRID<VECTOR<float,3> > > >::LOOKUP>(BLOCK_UNIFORM<GRID<VECTOR<float,3> > > const&,FACE_LOOKUP_COLLIDABLE_SLIP_UNIFORM<GRID<VECTOR<float,3> >,FACE_LOOKUP_UNIFORM<GRID<VECTOR<float,3> > > >::LOOKUP const&,VECTOR<float,3> const&);
template ARRAY<PAIR<FACE_INDEX<3>,float> > LINEAR_INTERPOLATION_MAC_3D_HELPER<GRID<VECTOR<float,3> > >::Interpolate_Face_Y_Transformed_Weights<BLOCK_UNIFORM<GRID<VECTOR<float,3> > >,FACE_LOOKUP_COLLIDABLE_SLIP_UNIFORM<GRID<VECTOR<float,3> >,FACE_LOOKUP_UNIFORM<GRID<VECTOR<float,3> > > >::LOOKUP>(BLOCK_UNIFORM<GRID<VECTOR<float,3> > > const&,FACE_LOOKUP_COLLIDABLE_SLIP_UNIFORM<GRID<VECTOR<float,3> >,FACE_LOOKUP_UNIFORM<GRID<VECTOR<float,3> > > >::LOOKUP const&,VECTOR<float,3> const&);
template ARRAY<PAIR<FACE_INDEX<3>,float> > LINEAR_INTERPOLATION_MAC_3D_HELPER<GRID<VECTOR<float,3> > >::Interpolate_Face_Z_Transformed_Weights<BLOCK_UNIFORM<GRID<VECTOR<float,3> > >,FACE_LOOKUP_COLLIDABLE_SLIP_UNIFORM<GRID<VECTOR<float,3> >,FACE_LOOKUP_UNIFORM<GRID<VECTOR<float,3> > > >::LOOKUP>(BLOCK_UNIFORM<GRID<VECTOR<float,3> > > const&,FACE_LOOKUP_COLLIDABLE_SLIP_UNIFORM<GRID<VECTOR<float,3> >,FACE_LOOKUP_UNIFORM<GRID<VECTOR<float,3> > > >::LOOKUP const&,VECTOR<float,3> const&);
#ifndef COMPILE_WITHOUT_DOUBLE_SUPPORT
template class FACE_LOOKUP_COLLIDABLE_SLIP_UNIFORM<GRID<VECTOR<double,1> >,FACE_LOOKUP_UNIFORM<GRID<VECTOR<double,1> > > >;
template class FACE_LOOKUP_COLLIDABLE_SLIP_UNIFORM<GRID<VECTOR<double,2> >,FACE_LOOKUP_UNIFORM<GRID<VECTOR<double,2> > > >;
template class FACE_LOOKUP_COLLIDABLE_SLIP_UNIFORM<GRID<VECTOR<double,3> >,FACE_LOOKUP_UNIFORM<GRID<VECTOR<double,3> > > >;
template VECTOR<double,2> LINEAR_INTERPOLATION_MAC_2D_HELPER<GRID<VECTOR<double,2> > >::Extrema_Face_X_Transformed<BLOCK_UNIFORM<GRID<VECTOR<double,2> > >,FACE_LOOKUP_COLLIDABLE_SLIP_UNIFORM<GRID<VECTOR<double,2> >,FACE_LOOKUP_UNIFORM<GRID<VECTOR<double,2> > > >::LOOKUP>(BLOCK_UNIFORM<GRID<VECTOR<double,2> > > const&,FACE_LOOKUP_COLLIDABLE_SLIP_UNIFORM<GRID<VECTOR<double,2> >,FACE_LOOKUP_UNIFORM<GRID<VECTOR<double,2> > > >::LOOKUP const&,FACE_LOOKUP_COLLIDABLE_SLIP_UNIFORM<GRID<VECTOR<double,2> >,FACE_LOOKUP_UNIFORM<GRID<VECTOR<double,2> > > >::LOOKUP const&,VECTOR<double,2> const&);
template VECTOR<double,2> LINEAR_INTERPOLATION_MAC_2D_HELPER<GRID<VECTOR<double,2> > >::Extrema_Face_Y_Transformed<BLOCK_UNIFORM<GRID<VECTOR<double,2> > >,FACE_LOOKUP_COLLIDABLE_SLIP_UNIFORM<GRID<VECTOR<double,2> >,FACE_LOOKUP_UNIFORM<GRID<VECTOR<double,2> > > >::LOOKUP>(BLOCK_UNIFORM<GRID<VECTOR<double,2> > > const&,FACE_LOOKUP_COLLIDABLE_SLIP_UNIFORM<GRID<VECTOR<double,2> >,FACE_LOOKUP_UNIFORM<GRID<VECTOR<double,2> > > >::LOOKUP const&,FACE_LOOKUP_COLLIDABLE_SLIP_UNIFORM<GRID<VECTOR<double,2> >,FACE_LOOKUP_UNIFORM<GRID<VECTOR<double,2> > > >::LOOKUP const&,VECTOR<double,2> const&);
template VECTOR<double,2> LINEAR_INTERPOLATION_MAC_3D_HELPER<GRID<VECTOR<double,3> > >::Extrema_Face_X_Transformed<BLOCK_UNIFORM<GRID<VECTOR<double,3> > >,FACE_LOOKUP_COLLIDABLE_SLIP_UNIFORM<GRID<VECTOR<double,3> >,FACE_LOOKUP_UNIFORM<GRID<VECTOR<double,3> > > >::LOOKUP>(BLOCK_UNIFORM<GRID<VECTOR<double,3> > > const&,FACE_LOOKUP_COLLIDABLE_SLIP_UNIFORM<GRID<VECTOR<double,3> >,FACE_LOOKUP_UNIFORM<GRID<VECTOR<double,3> > > >::LOOKUP const&,FACE_LOOKUP_COLLIDABLE_SLIP_UNIFORM<GRID<VECTOR<double,3> >,FACE_LOOKUP_UNIFORM<GRID<VECTOR<double,3> > > >::LOOKUP const&,VECTOR<double,3> const&);
template VECTOR<double,2> LINEAR_INTERPOLATION_MAC_3D_HELPER<GRID<VECTOR<double,3> > >::Extrema_Face_Y_Transformed<BLOCK_UNIFORM<GRID<VECTOR<double,3> > >,FACE_LOOKUP_COLLIDABLE_SLIP_UNIFORM<GRID<VECTOR<double,3> >,FACE_LOOKUP_UNIFORM<GRID<VECTOR<double,3> > > >::LOOKUP>(BLOCK_UNIFORM<GRID<VECTOR<double,3> > > const&,FACE_LOOKUP_COLLIDABLE_SLIP_UNIFORM<GRID<VECTOR<double,3> >,FACE_LOOKUP_UNIFORM<GRID<VECTOR<double,3> > > >::LOOKUP const&,FACE_LOOKUP_COLLIDABLE_SLIP_UNIFORM<GRID<VECTOR<double,3> >,FACE_LOOKUP_UNIFORM<GRID<VECTOR<double,3> > > >::LOOKUP const&,VECTOR<double,3> const&);
template VECTOR<double,2> LINEAR_INTERPOLATION_MAC_3D_HELPER<GRID<VECTOR<double,3> > >::Extrema_Face_Z_Transformed<BLOCK_UNIFORM<GRID<VECTOR<double,3> > >,FACE_LOOKUP_COLLIDABLE_SLIP_UNIFORM<GRID<VECTOR<double,3> >,FACE_LOOKUP_UNIFORM<GRID<VECTOR<double,3> > > >::LOOKUP>(BLOCK_UNIFORM<GRID<VECTOR<double,3> > > const&,FACE_LOOKUP_COLLIDABLE_SLIP_UNIFORM<GRID<VECTOR<double,3> >,FACE_LOOKUP_UNIFORM<GRID<VECTOR<double,3> > > >::LOOKUP const&,FACE_LOOKUP_COLLIDABLE_SLIP_UNIFORM<GRID<VECTOR<double,3> >,FACE_LOOKUP_UNIFORM<GRID<VECTOR<double,3> > > >::LOOKUP const&,VECTOR<double,3> const&);
template double LINEAR_INTERPOLATION_MAC_2D_HELPER<GRID<VECTOR<double,2> > >::Interpolate_Face_X_Transformed<BLOCK_UNIFORM<GRID<VECTOR<double,2> > >,FACE_LOOKUP_COLLIDABLE_SLIP_UNIFORM<GRID<VECTOR<double,2> >,FACE_LOOKUP_UNIFORM<GRID<VECTOR<double,2> > > >::LOOKUP>(BLOCK_UNIFORM<GRID<VECTOR<double,2> > > const&,FACE_LOOKUP_COLLIDABLE_SLIP_UNIFORM<GRID<VECTOR<double,2> >,FACE_LOOKUP_UNIFORM<GRID<VECTOR<double,2> > > >::LOOKUP const&,VECTOR<double,2> const&);
template double LINEAR_INTERPOLATION_MAC_2D_HELPER<GRID<VECTOR<double,2> > >::Interpolate_Face_Y_Transformed<BLOCK_UNIFORM<GRID<VECTOR<double,2> > >,FACE_LOOKUP_COLLIDABLE_SLIP_UNIFORM<GRID<VECTOR<double,2> >,FACE_LOOKUP_UNIFORM<GRID<VECTOR<double,2> > > >::LOOKUP>(BLOCK_UNIFORM<GRID<VECTOR<double,2> > > const&,FACE_LOOKUP_COLLIDABLE_SLIP_UNIFORM<GRID<VECTOR<double,2> >,FACE_LOOKUP_UNIFORM<GRID<VECTOR<double,2> > > >::LOOKUP const&,VECTOR<double,2> const&);
template double LINEAR_INTERPOLATION_MAC_3D_HELPER<GRID<VECTOR<double,3> > >::Interpolate_Face_X_Transformed<BLOCK_UNIFORM<GRID<VECTOR<double,3> > >,FACE_LOOKUP_COLLIDABLE_SLIP_UNIFORM<GRID<VECTOR<double,3> >,FACE_LOOKUP_UNIFORM<GRID<VECTOR<double,3> > > >::LOOKUP>(BLOCK_UNIFORM<GRID<VECTOR<double,3> > > const&,FACE_LOOKUP_COLLIDABLE_SLIP_UNIFORM<GRID<VECTOR<double,3> >,FACE_LOOKUP_UNIFORM<GRID<VECTOR<double,3> > > >::LOOKUP const&,VECTOR<double,3> const&);
template double LINEAR_INTERPOLATION_MAC_3D_HELPER<GRID<VECTOR<double,3> > >::Interpolate_Face_Y_Transformed<BLOCK_UNIFORM<GRID<VECTOR<double,3> > >,FACE_LOOKUP_COLLIDABLE_SLIP_UNIFORM<GRID<VECTOR<double,3> >,FACE_LOOKUP_UNIFORM<GRID<VECTOR<double,3> > > >::LOOKUP>(BLOCK_UNIFORM<GRID<VECTOR<double,3> > > const&,FACE_LOOKUP_COLLIDABLE_SLIP_UNIFORM<GRID<VECTOR<double,3> >,FACE_LOOKUP_UNIFORM<GRID<VECTOR<double,3> > > >::LOOKUP const&,VECTOR<double,3> const&);
template double LINEAR_INTERPOLATION_MAC_3D_HELPER<GRID<VECTOR<double,3> > >::Interpolate_Face_Z_Transformed<BLOCK_UNIFORM<GRID<VECTOR<double,3> > >,FACE_LOOKUP_COLLIDABLE_SLIP_UNIFORM<GRID<VECTOR<double,3> >,FACE_LOOKUP_UNIFORM<GRID<VECTOR<double,3> > > >::LOOKUP>(BLOCK_UNIFORM<GRID<VECTOR<double,3> > > const&,FACE_LOOKUP_COLLIDABLE_SLIP_UNIFORM<GRID<VECTOR<double,3> >,FACE_LOOKUP_UNIFORM<GRID<VECTOR<double,3> > > >::LOOKUP const&,VECTOR<double,3> const&);
template ARRAY<PAIR<FACE_INDEX<2>,double> > LINEAR_INTERPOLATION_MAC_2D_HELPER<GRID<VECTOR<double,2> > >::Interpolate_Face_X_Transformed_Weights<BLOCK_UNIFORM<GRID<VECTOR<double,2> > >,FACE_LOOKUP_COLLIDABLE_SLIP_UNIFORM<GRID<VECTOR<double,2> >,FACE_LOOKUP_UNIFORM<GRID<VECTOR<double,2> > > >::LOOKUP>(BLOCK_UNIFORM<GRID<VECTOR<double,2> > > const&,FACE_LOOKUP_COLLIDABLE_SLIP_UNIFORM<GRID<VECTOR<double,2> >,FACE_LOOKUP_UNIFORM<GRID<VECTOR<double,2> > > >::LOOKUP const&,VECTOR<double,2> const&);
template ARRAY<PAIR<FACE_INDEX<2>,double> > LINEAR_INTERPOLATION_MAC_2D_HELPER<GRID<VECTOR<double,2> > >::Interpolate_Face_Y_Transformed_Weights<BLOCK_UNIFORM<GRID<VECTOR<double,2> > >,FACE_LOOKUP_COLLIDABLE_SLIP_UNIFORM<GRID<VECTOR<double,2> >,FACE_LOOKUP_UNIFORM<GRID<VECTOR<double,2> > > >::LOOKUP>(BLOCK_UNIFORM<GRID<VECTOR<double,2> > > const&,FACE_LOOKUP_COLLIDABLE_SLIP_UNIFORM<GRID<VECTOR<double,2> >,FACE_LOOKUP_UNIFORM<GRID<VECTOR<double,2> > > >::LOOKUP const&,VECTOR<double,2> const&);
template ARRAY<PAIR<FACE_INDEX<3>,double> > LINEAR_INTERPOLATION_MAC_3D_HELPER<GRID<VECTOR<double,3> > >::Interpolate_Face_X_Transformed_Weights<BLOCK_UNIFORM<GRID<VECTOR<double,3> > >,FACE_LOOKUP_COLLIDABLE_SLIP_UNIFORM<GRID<VECTOR<double,3> >,FACE_LOOKUP_UNIFORM<GRID<VECTOR<double,3> > > >::LOOKUP>(BLOCK_UNIFORM<GRID<VECTOR<double,3> > > const&,FACE_LOOKUP_COLLIDABLE_SLIP_UNIFORM<GRID<VECTOR<double,3> >,FACE_LOOKUP_UNIFORM<GRID<VECTOR<double,3> > > >::LOOKUP const&,VECTOR<double,3> const&);
template ARRAY<PAIR<FACE_INDEX<3>,double> > LINEAR_INTERPOLATION_MAC_3D_HELPER<GRID<VECTOR<double,3> > >::Interpolate_Face_Y_Transformed_Weights<BLOCK_UNIFORM<GRID<VECTOR<double,3> > >,FACE_LOOKUP_COLLIDABLE_SLIP_UNIFORM<GRID<VECTOR<double,3> >,FACE_LOOKUP_UNIFORM<GRID<VECTOR<double,3> > > >::LOOKUP>(BLOCK_UNIFORM<GRID<VECTOR<double,3> > > const&,FACE_LOOKUP_COLLIDABLE_SLIP_UNIFORM<GRID<VECTOR<double,3> >,FACE_LOOKUP_UNIFORM<GRID<VECTOR<double,3> > > >::LOOKUP const&,VECTOR<double,3> const&);
template ARRAY<PAIR<FACE_INDEX<3>,double> > LINEAR_INTERPOLATION_MAC_3D_HELPER<GRID<VECTOR<double,3> > >::Interpolate_Face_Z_Transformed_Weights<BLOCK_UNIFORM<GRID<VECTOR<double,3> > >,FACE_LOOKUP_COLLIDABLE_SLIP_UNIFORM<GRID<VECTOR<double,3> >,FACE_LOOKUP_UNIFORM<GRID<VECTOR<double,3> > > >::LOOKUP>(BLOCK_UNIFORM<GRID<VECTOR<double,3> > > const&,FACE_LOOKUP_COLLIDABLE_SLIP_UNIFORM<GRID<VECTOR<double,3> >,FACE_LOOKUP_UNIFORM<GRID<VECTOR<double,3> > > >::LOOKUP const&,VECTOR<double,3> const&);
#endif
