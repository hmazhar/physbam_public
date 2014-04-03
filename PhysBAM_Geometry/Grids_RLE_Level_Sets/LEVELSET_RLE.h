//#####################################################################
// Copyright 2005, Eran Guendelman, Geoffrey Irving, Jerry Talton.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class LEVELSET_RLE  
//#####################################################################
#ifndef COMPILE_WITHOUT_RLE_SUPPORT
#ifndef __LEVELSET_RLE__
#define __LEVELSET_RLE__

#include <PhysBAM_Tools/Grids_RLE/RLE_GRID_2D.h>
#include <PhysBAM_Tools/Grids_RLE/RLE_GRID_3D.h>
#include <PhysBAM_Tools/Grids_RLE_Arrays/GRID_ARRAYS_POLICY_RLE.h>
#include <PhysBAM_Tools/Grids_RLE_Boundaries/BOUNDARY_POLICY_RLE.h>
#include <PhysBAM_Tools/Grids_RLE_Interpolation/INTERPOLATION_POLICY_RLE.h>
#include <PhysBAM_Tools/Grids_RLE_Interpolation/LINEAR_INTERPOLATION_RLE_HELPER.h>
#include <PhysBAM_Geometry/Grids_RLE_Collisions/GRID_BASED_COLLISION_GEOMETRY_COLLECTION_POLICY_RLE.h>
#include <PhysBAM_Geometry/Grids_RLE_Interpolation_Collidable/INTERPOLATION_COLLIDABLE_POLICY_RLE.h>
#include <PhysBAM_Geometry/Level_Sets/LEVELSET.h>
namespace PhysBAM{

template<class T_GRID>
class LEVELSET_RLE:public LEVELSET<typename T_GRID::SCALAR,T_GRID>
{
public:
    typedef typename T_GRID::SCALAR T;typedef typename T_GRID::VECTOR_T TV;typedef typename T_GRID::VECTOR_INT TV_INT;typedef typename T_GRID::BLOCK T_BLOCK;
    typedef typename T_GRID::LINEAR_INTERPOLATION_MAC_HELPER T_LINEAR_INTERPOLATION_MAC_HELPER;typedef typename INTERPOLATION_POLICY<T_GRID>::FACE_LOOKUP T_FACE_LOOKUP;typedef typename INTERPOLATION_COLLIDABLE_POLICY<T_GRID>::AVERAGING T_AVERAGING;
    typedef typename T_GRID::BLOCK_ITERATOR BLOCK_ITERATOR;typedef typename T_GRID::CELL_ITERATOR CELL_ITERATOR;
    typedef typename T_GRID::ARRAYS_HORIZONTAL::template REBIND<int>::TYPE T_ARRAYS_HORIZONTAL_INT;

    typedef LEVELSET<T,T_GRID> BASE;
    using BASE::small_number;using BASE::boundary;using BASE::max_time_step;
    using BASE::curvature_motion;using BASE::sigma;
    using BASE::collision_body_list;using BASE::valid_mask_current;using BASE::collidable_phi_replacement_value;using BASE::secondary_interpolation;

    T_GRID& grid;
    ARRAY<T>& phi;
    ARRAY<TV>* normals;
    ARRAY<T>* curvature;
    ARRAY<T> *block_minimum,*block_maximum;
    T half_bandwidth;

    LEVELSET_RLE(T_GRID& grid_input,ARRAY<T>& phi_input);
    virtual ~LEVELSET_RLE();

    void Initialize_Levelset_Grid_Values()
    {LEVELSET<T,T_GRID>::Initialize_Levelset_Grid_Values(grid);}

    void Set_Band_Width(const T number_of_cells)
    {half_bandwidth=(T).5*number_of_cells*grid.uniform_grid.dX.Max();}

    T Phi(const TV& X) const
    {return Phi(phi,X);}

    T Phi(const ARRAY<T>& phi,const TV& X) const
    {T_BLOCK block(grid,X);return block?T_LINEAR_INTERPOLATION_MAC_HELPER::Interpolate_Cell(block,phi,X):grid.positive_bandwidth;}

    T Phi(const T_BLOCK& block,const TV& X) const
    {return T_LINEAR_INTERPOLATION_MAC_HELPER::Interpolate_Cell(block,phi,X);}

    T Clamped_Phi(const TV& X) const
    {return Phi(grid.Clamp(X));}

    T Phi_Far(const TV& X) const
    {return Phi_Far(T_BLOCK(grid,X),X);}

    T Phi_Far(const T_BLOCK& block,const TV& X) const
    {if(block) return Phi(block,X);
    CELL_ITERATOR cell(grid,X);return phi(cell.Cell());}

    T Extended_Phi(const TV& X) const
    {TV clamped_X=grid.Clamp(X);return Phi_Far(clamped_X)+(clamped_X-X).Magnitude();}

    TV Normal(const T_BLOCK& block,const TV& X) const
    {assert(normals);return T_LINEAR_INTERPOLATION_MAC_HELPER::Interpolate_Cell(block,*normals,X).Normalized();}

//#####################################################################
    void Clean_Memory();
    T Collision_Aware_Phi(const T_BLOCK& block,const TV& location) const;
    TV Normal(const TV& X) const;
    TV Extended_Normal(const TV& X) const;
    T Phi_Secondary(const TV& X) const;
    T CFL(const ARRAY<T>& V) const;
    T Curvature_CFL() const;
    void Compute_Normals(const T time=0);
    void Compute_Normals(ARRAY<T>& phi_ghost);
    void Compute_Curvature(const T time=0);
    void Compute_Curvature(const ARRAY<T>& phi_ghost);
public:
    void Compute_Block_Minimum_And_Maximum(const bool recompute_if_exists=true);
    void Curvature_Motion(const T dt,const T time);
    void Curvature_Motion(const T dt,const ARRAY<T>& phi_ghost);
public:
    void Fast_Marching_Method(const T time);
    T Approximate_Negative_Size(const T interface_thickness=3,const T time=0) const;
    T Approximate_Positive_Size(const T interface_thickness=3,const T time=0) const;
    void Transfer_Phi(const T_GRID& new_grid);
    void Rebuild_Grid(const int extra_cells,const T_ARRAYS_HORIZONTAL_INT* ground_j);
//#####################################################################
};
}
#endif
#endif
