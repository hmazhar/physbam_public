//#####################################################################
// Copyright 2004-2008, Ron Fedkiw, Eran Guendelman, Geoffrey Irving, Frank Losasso, Andrew Selle.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class LEVELSET_MAKER
//##################################################################### 
#ifndef __LEVELSET_MAKER__
#define __LEVELSET_MAKER__

#include <PhysBAM_Tools/Arrays/ARRAY.h>
#include <PhysBAM_Tools/Arrays/ARRAYS_FORWARD.h>
#include <PhysBAM_Tools/Utilities/NONCOPYABLE.h>
#include <PhysBAM_Tools/Vectors/VECTOR_3D.h>
namespace PhysBAM{

template<class T>
class LEVELSET_MAKER:public NONCOPYABLE
{
protected:
    typedef VECTOR<int,3> TV_INT;
public:
    bool verbose;
    bool compute_signed_distance_function,compute_unsigned_distance_function,compute_heaviside_function;
    bool use_fmm;
    T fmm_one_sided_band_width; // in number of cells (if zero then compute phi everywhere)
    bool extrapolate_velocity;
    T velocity_extrapolation_one_sided_band_width; // in number of cells
    bool only_boundary_region_is_outside;
    bool keep_only_largest_inside_region;
    bool flip_sign_if_corners_are_inside;
    T surface_thickness_or_zero; // if zero, use grid.dx/100
    T surface_padding_for_flood_fill_or_negative; // if negative, use grid.dx/100 (0 means no padding)
    ARRAY<VECTOR<int,3> > initialized_indices;
    ARRAY<int> initialized_indices_octree;
    bool write_debug_data;
    bool write_debug_path;
    VECTOR<int,3> path_start_node,path_end_node;
    bool use_orthogonal_vote;
    int positive_boundary_band;
    T remove_degenerate_triangles_area_threshold;
    T phi_offset;
   
    LEVELSET_MAKER()
    {
        Verbose_Mode(false);
        Write_Debug_Data(false);
        Write_Debug_Path(false);
        Compute_Signed_Distance_Function();
        Use_Fast_Marching_Method();
        Extrapolate_Velocity(false);
        Only_Boundary_Region_Is_Outside(false);
        Keep_Only_Largest_Inside_Region(false);
        Flip_Sign_If_Corners_Are_Inside(false);
        Set_Surface_Thickness();
        Set_Surface_Padding_For_Flood_Fill();
        Use_Orthogonal_Vote(false);
        Set_Positive_Boundary_Band(0);
        Use_Remove_Degenerate_Triangles();
        Set_Phi_Offset();
    }

    void Verbose_Mode(const bool verbose_input=true)
    {verbose=verbose_input;}

    void Write_Debug_Data(const bool write_debug_data_input=true)
    {write_debug_data=write_debug_data_input;}

    void Write_Debug_Path(const bool write_debug_path_input=true,const VECTOR<int,3>& path_start_node_input=TV_INT(),const VECTOR<int,3>& path_end_node_input=TV_INT())
    {write_debug_path=write_debug_path_input;if(write_debug_path){path_start_node=path_start_node_input;path_end_node=path_end_node_input;}}

    void Compute_Signed_Distance_Function()
    {compute_signed_distance_function=true;compute_unsigned_distance_function=false;compute_heaviside_function=false;}

    void Compute_Unsigned_Distance_Function()
    {compute_signed_distance_function=false;compute_unsigned_distance_function=true;compute_heaviside_function=false;}

    void Compute_Heaviside_Function()
    {compute_signed_distance_function=false;compute_unsigned_distance_function=false;compute_heaviside_function=true;}

    // fmm_one_sided_band_width_input==0 indicates an infinite band width (run fmm everywhere)
    void Use_Fast_Marching_Method(const bool use_fmm_input=true,const T fmm_one_sided_band_width_input=0)
    {use_fmm=use_fmm_input;if(use_fmm)fmm_one_sided_band_width=fmm_one_sided_band_width_input;}

    void Extrapolate_Velocity(const bool extrapolate_velocity_input=true,const T velocity_extrapolation_one_sided_band_width_input=3)
    {extrapolate_velocity=extrapolate_velocity_input;if(extrapolate_velocity)velocity_extrapolation_one_sided_band_width=velocity_extrapolation_one_sided_band_width_input;}

    void Only_Boundary_Region_Is_Outside(const bool only_boundary_region_is_outside_input=true)
    {only_boundary_region_is_outside=only_boundary_region_is_outside_input;}

    void Keep_Only_Largest_Inside_Region(const bool keep=true)
    {keep_only_largest_inside_region=keep;}

    void Flip_Sign_If_Corners_Are_Inside(const bool flip=true)
    {flip_sign_if_corners_are_inside=flip;}

    void Set_Surface_Thickness(const T surface_thickness_input=0) // zero means grid.dx/100
    {surface_thickness_or_zero=surface_thickness_input;}

    void Set_Surface_Padding_For_Flood_Fill(const T surface_padding_for_flood_fill_input=-1) // negative means grid.dx/100
    {surface_padding_for_flood_fill_or_negative=surface_padding_for_flood_fill_input;}

    void Use_Orthogonal_Vote(const bool use_orthogonal_vote_input=true)
    {use_orthogonal_vote=use_orthogonal_vote_input;}

    void Set_Positive_Boundary_Band(const int positive_boundary_band_input=0)
    {positive_boundary_band=positive_boundary_band_input;}

    void Use_Remove_Degenerate_Triangles(const bool remove_degenerate_triangles_input=true,const T area_threshold_input=1e-15)
    {if(!remove_degenerate_triangles_input) remove_degenerate_triangles_area_threshold=0;
    else remove_degenerate_triangles_area_threshold=area_threshold_input;}

    void Set_Phi_Offset(const T phi_offset_input=0)
    {phi_offset=phi_offset_input;}

protected:
    template<class T_GRID> T Surface_Thickness_Over_Two(const T_GRID& grid) const
    {return (T).5*(surface_thickness_or_zero?surface_thickness_or_zero:grid.Minimum_Edge_Length()/100);}

    template<class T_GRID> T Surface_Padding_For_Flood_Fill(const T_GRID& grid) const
    {return surface_padding_for_flood_fill_or_negative>=0?surface_padding_for_flood_fill_or_negative:grid.Minimum_Edge_Length()/100;}
};
}
#endif
