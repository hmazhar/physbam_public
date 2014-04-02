//#####################################################################
// Copyright 2002-2007, Zhaosheng Bao, Robert Bridson, Ronald Fedkiw, Eran Guendelman, Geoffrey Irving, Sergey Koltakov, Neil Molino, Igor Neverov, Andrew Selle, Jerry Talton, Rachel Weinstein.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class LEVELSET_UNIFORM
//#####################################################################
#ifndef __LEVELSET_UNIFORM__
#define __LEVELSET_UNIFORM__

#include <PhysBAM_Tools/Grids_Uniform_Boundaries/BOUNDARY_POLICY_UNIFORM.h>
#include <PhysBAM_Tools/Grids_Uniform_Interpolation/FACE_LOOKUP_UNIFORM.h>
#include <PhysBAM_Tools/Grids_Uniform_Interpolation/LINEAR_INTERPOLATION_UNIFORM.h>
#include <PhysBAM_Tools/Parallel_Computation/THREAD_QUEUE.h>
#include <PhysBAM_Geometry/Grids_Uniform_Collisions/GRID_BASED_COLLISION_GEOMETRY_COLLECTION_POLICY_UNIFORM.h>
#include <PhysBAM_Geometry/Grids_Uniform_Interpolation_Collidable/INTERPOLATION_COLLIDABLE_POLICY_UNIFORM.h>
#include <PhysBAM_Geometry/Level_Sets/LEVELSET.h>
namespace PhysBAM{

template<class TV> class VOF_ADVECTION;

template<class T_GRID_input>
class LEVELSET_UNIFORM:public LEVELSET<typename T_GRID_input::SCALAR,T_GRID_input>
{
    typedef typename T_GRID_input::VECTOR_T TV;typedef typename T_GRID_input::SCALAR T;
    typedef typename T_GRID_input::VECTOR_INT TV_INT;typedef typename GRID_ARRAYS_POLICY<T_GRID_input>::ARRAYS_SCALAR T_ARRAYS_SCALAR;
    typedef typename T_ARRAYS_SCALAR::template REBIND<TV>::TYPE T_ARRAYS_VECTOR;typedef typename T_ARRAYS_SCALAR::template REBIND<bool>::TYPE T_ARRAYS_BOOL;
    typedef typename REBIND<T_ARRAYS_SCALAR,RANGE<VECTOR<T,1> > >::TYPE T_ARRAYS_INTERVAL;
    typedef typename GRID_ARRAYS_POLICY<T_GRID_input>::FACE_ARRAYS T_FACE_ARRAYS_SCALAR;typedef typename T_GRID_input::NODE_ITERATOR NODE_ITERATOR;typedef typename T_GRID_input::CELL_ITERATOR CELL_ITERATOR;
    typedef typename INTERPOLATION_POLICY<T_GRID_input>::INTERPOLATION_SCALAR T_INTERPOLATION_SCALAR;
public:
    typedef T_GRID_input T_GRID;
    typedef LEVELSET<T,T_GRID> BASE;
    using BASE::interpolation;using BASE::secondary_interpolation;using BASE::curvature_interpolation;using BASE::curvature_motion;using BASE::sigma;using BASE::max_time_step;
    using BASE::collision_body_list;using BASE::collidable_phi_replacement_value;using BASE::valid_mask_current;using BASE::boundary;

    T_GRID& grid;
    T_ARRAYS_SCALAR& phi;
    T_ARRAYS_VECTOR* normals;
    T_ARRAYS_SCALAR *curvature;
    T_ARRAYS_INTERVAL *cell_range;
    THREAD_QUEUE *thread_queue;
    int number_of_ghost_cells;
protected:

public:
    LEVELSET_UNIFORM(T_GRID& grid_input,T_ARRAYS_SCALAR& phi_input,const int number_of_ghost_cells_input);
    ~LEVELSET_UNIFORM();

    void Initialize_Levelset_Grid_Values()
    {LEVELSET<T,T_GRID>::Initialize_Levelset_Grid_Values(grid);}

    T Phi(const TV& location) const
    {return interpolation->Clamped_To_Array(grid,phi,location);}

    T Phi_Secondary(const TV& location) const
    {return secondary_interpolation->Clamped_To_Array(grid,phi,location);}

    T Extended_Phi(const TV& location) const
    {TV clamped_location(grid.Clamp(location));
    T phi_value=interpolation->Clamped_To_Array(grid,phi,clamped_location);
    T magnitude_squared=(clamped_location-location).Magnitude_Squared();
    if(magnitude_squared) phi_value=sqrt(magnitude_squared+sqr(max((T)0,phi_value)));
    return phi_value;}

    T Curvature(const TV& location) const // later add finite difference for curvature, like in normal above
    {assert(curvature);return curvature_interpolation->Clamped_To_Array(grid,*curvature,location);}

    bool Lazy_Inside(const TV& clamped_X,const T contour_value=0) const
    {assert(cell_range);TV_INT index=T_INTERPOLATION_SCALAR::Clamped_Index_End_Minus_One(grid,phi,clamped_X);
    const RANGE<VECTOR<T,1> >& range=(*cell_range)(index);
    if(range.min_corner.x>contour_value) return false;
    else if(range.max_corner.x<=contour_value) return true;
    return interpolation->From_Base_Node(grid,phi,clamped_X,index)<=contour_value;}

    bool Lazy_Inside_And_Value(const TV& clamped_X,T& phi_value,const T contour_value=0) const
    {assert(cell_range);TV_INT index=T_INTERPOLATION_SCALAR::Clamped_Index_End_Minus_One(grid,phi,clamped_X);
    if((*cell_range)(index).min_corner.x>contour_value) return false;
    phi_value=interpolation->From_Base_Node(grid,phi,clamped_X,index);
    return phi_value<=contour_value;}

    bool Lazy_Inside_Extended_Levelset(const TV& unclamped_X,const T contour_value=0) const
    {if(grid.domain.Lazy_Inside(unclamped_X)) return Lazy_Inside(unclamped_X,contour_value);
    if(contour_value<=0) return false;
    T sqr_contour_value=sqr(contour_value);
    TV clamped_X(grid.Clamp(unclamped_X));
    T magnitude_squared=(unclamped_X-clamped_X).Magnitude_Squared();
    if(magnitude_squared>sqr_contour_value) return false;
    TV_INT index=T_INTERPOLATION_SCALAR::Clamped_Index_End_Minus_One(grid,phi,clamped_X);
    T phi_value=interpolation->From_Base_Node(grid,phi,clamped_X,index);
    return magnitude_squared+sqr(max((T)0,phi_value))<=sqr_contour_value;}

    // phi_value only guaranteed to be set if function returns true
    bool Lazy_Inside_Extended_Levelset_And_Value(const TV& unclamped_X,T& phi_value,const T contour_value=0) const
    {assert(cell_range);TV clamped_X(grid.Clamp(unclamped_X));
    T magnitude_squared=(unclamped_X-clamped_X).Magnitude_Squared();
    if(contour_value <= 0){if(magnitude_squared > 0) return false;}
    else if(magnitude_squared > sqr(contour_value)) return false;
    TV_INT index=T_INTERPOLATION_SCALAR::Clamped_Index_End_Minus_One(grid,phi,clamped_X);
    if((*cell_range)(index).min_corner.x>contour_value) return false;
    phi_value=interpolation->From_Base_Node(grid,phi,clamped_X,index);
    if(magnitude_squared>0) phi_value=sqrt(magnitude_squared+sqr(max((T)0,phi_value)));
    return phi_value<=contour_value;}

    bool Lazy_Outside(const TV& clamped_X,const T contour_value=0) const
    {return !Lazy_Inside(clamped_X,contour_value);}

    bool Lazy_Outside_Extended_Levelset(const TV& unclamped_X,const T contour_value=0) const
    {return !Lazy_Inside_Extended_Levelset(unclamped_X,contour_value);}

    // phi_value only guaranteed to be set if function returns true
    bool Lazy_Outside_Extended_Levelset_And_Value(const TV& unclamped_X,T& phi_value,const T contour_value=0) const
    {assert(cell_range);TV clamped_X(grid.Clamp(unclamped_X));
    T magnitude_squared=(unclamped_X-clamped_X).Magnitude_Squared();
    TV_INT index=T_INTERPOLATION_SCALAR::Clamped_Index_End_Minus_One(grid,phi,clamped_X);
    if(magnitude_squared==0 && (*cell_range)(index).max_corner.x<=contour_value) return false;
    phi_value=interpolation->From_Base_Node(grid,phi,clamped_X,index);
    if(magnitude_squared>0) phi_value=sqrt(magnitude_squared+sqr(max((T)0,phi_value)));
    return phi_value>contour_value;}

//#####################################################################
    T Collision_Aware_Phi(const TV& location) const;
    T CFL(const T_FACE_ARRAYS_SCALAR& face_velocities) const;
    T CFL(const T_ARRAYS_VECTOR& velocity) const;
    TV Iterative_Find_Interface(TV left,TV right,const int iterations=3) const;
    void Compute_Gradient(T_ARRAYS_VECTOR& gradient,const T time=0) const;
//#####################################################################
};
}
#endif
