//#####################################################################
// Copyright 2004-2008, Ron Fedkiw, Eran Guendelman, Geoffrey Irving, Frank Losasso, Andrew Selle.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class LEVELSET_MAKER_DYADIC
//##################################################################### 
#ifndef COMPILE_WITHOUT_DYADIC_SUPPORT
#ifndef __LEVELSET_MAKER_DYADIC__
#define __LEVELSET_MAKER_DYADIC__

#include <PhysBAM_Geometry/Level_Sets/LEVELSET_MAKER.h>
#include <PhysBAM_Geometry/Topology_Based_Geometry/TOPOLOGY_BASED_GEOMETRY_FORWARD.h>
namespace PhysBAM{
template<class TV> class GRID;
template<class T> class OCTREE_GRID;

template<class T>
class LEVELSET_MAKER_DYADIC:public LEVELSET_MAKER<T>
{
    typedef LEVELSET_MAKER<T> BASE;
    typedef VECTOR<int,3> TV_INT;

    using BASE::verbose;using BASE::compute_signed_distance_function;using BASE::compute_unsigned_distance_function;using BASE::compute_heaviside_function;
    using BASE::use_fmm;using BASE::fmm_one_sided_band_width;using BASE::extrapolate_velocity;using BASE::velocity_extrapolation_one_sided_band_width;
    using BASE::only_boundary_region_is_outside;using BASE::keep_only_largest_inside_region;using BASE::flip_sign_if_corners_are_inside;
    using BASE::surface_thickness_or_zero;using BASE::surface_padding_for_flood_fill_or_negative;using BASE::initialized_indices;
    using BASE::initialized_indices_octree;using BASE::write_debug_data;using BASE::write_debug_path;using BASE::path_start_node;using BASE::path_end_node;
    using BASE::use_orthogonal_vote;using BASE::positive_boundary_band;using BASE::remove_degenerate_triangles_area_threshold;using BASE::phi_offset;
public:
//#####################################################################
    bool Compute_Level_Set(TRIANGULATED_SURFACE<T>& trianglated_surface,OCTREE_GRID<T>& grid,ARRAY<T>& phi,ARRAY<T>* phi_nodes=0,ARRAY<VECTOR<T,3> >* velocity=0,
        const T coarsen_bandwidth=0,int (*depth_function)(int,const RANGE<VECTOR<T,3> >&,void*)=0,void* depth_data=0);
//#####################################################################
};
}
#endif
#endif
