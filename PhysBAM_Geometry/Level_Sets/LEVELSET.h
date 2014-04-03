//#####################################################################
// Copyright 2002-2006, Doug Enright, Ronald Fedkiw, Frederic Gibou, Geoffrey Irving, Frank Losasso, Neil Molino, Avi Robinson-Mosher, Tamar Shinar, Jerry Talton.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class LEVELSET
//#####################################################################
#ifndef __LEVELSET__
#define __LEVELSET__

#include <PhysBAM_Tools/Arrays/ARRAYS_FORWARD.h>
#include <PhysBAM_Tools/Math_Tools/constants.h>
#include <PhysBAM_Tools/Utilities/NONCOPYABLE.h>
#include <PhysBAM_Tools/Utilities/PHYSBAM_OVERRIDE.h>
#include <PhysBAM_Tools/Vectors/SCALAR_POLICY.h>
#include <PhysBAM_Geometry/Interpolation_Collidable/INTERPOLATION_COLLIDABLE_POLICY.h>
#include <cassert>
#include <cfloat>
namespace PhysBAM{

template<class T_GRID> struct COLLISION_GEOMETRY_COLLECTION_POLICY;
template<class T_GRID> struct INTERPOLATION_POLICY;
template<class T_GRID> class LEVELSET_CALLBACKS;
template<class T_GRID> class LEVELSET_NORMAL_COMPUTATION;
template<class T_GRID> struct BOUNDARY_POLICY;
template<class T_GRID> struct GRID_ARRAYS_POLICY;
template<class TV> class GRID;

template<class T,class T_GRID=GRID<VECTOR<T,1> > >
class LEVELSET:public NONCOPYABLE
{
    STATIC_ASSERT((IS_SCALAR<T>::value));
    typedef typename T_GRID::VECTOR_T TV;typedef typename GRID_ARRAYS_POLICY<T_GRID>::ARRAYS_SCALAR T_ARRAYS_SCALAR;typedef typename T_ARRAYS_SCALAR::template REBIND<bool>::TYPE T_ARRAYS_BOOL;
    typedef typename COLLISION_GEOMETRY_COLLECTION_POLICY<T_GRID>::GRID_BASED_COLLISION_GEOMETRY T_GRID_BASED_COLLISION_GEOMETRY;
    typedef typename INTERPOLATION_POLICY<T_GRID>::LINEAR_INTERPOLATION_SCALAR T_LINEAR_INTERPOLATION_SCALAR;
    typedef typename INTERPOLATION_COLLIDABLE_POLICY<T_GRID>::LINEAR_INTERPOLATION_COLLIDABLE_CELL_SCALAR T_LINEAR_INTERPOLATION_COLLIDABLE_CELL_SCALAR;
    typedef typename INTERPOLATION_POLICY<T_GRID>::INTERPOLATION_SCALAR T_INTERPOLATION_SCALAR;
    typedef typename REBIND<T_INTERPOLATION_SCALAR,TV>::TYPE T_INTERPOLATION_VECTOR;
    typedef typename INTERPOLATION_POLICY<T_GRID>::FACE_LOOKUP T_FACE_LOOKUP;typedef typename INTERPOLATION_COLLIDABLE_POLICY<T_GRID>::FACE_LOOKUP_COLLIDABLE T_FACE_LOOKUP_COLLIDABLE;
    typedef typename REBIND<typename GRID_ARRAYS_POLICY<T_GRID>::FACE_ARRAYS,bool>::TYPE T_FACE_ARRAYS_BOOL;typedef typename BOUNDARY_POLICY<T_GRID>::BOUNDARY_SCALAR T_BOUNDARY_SCALAR;
    typedef typename INTERPOLATION_COLLIDABLE_POLICY<T_GRID>::FACE_LOOKUP_COLLIDABLE_SLIP T_FACE_LOOKUP_COLLIDABLE_SLIP;
public:
    T small_number;
    T_INTERPOLATION_SCALAR *interpolation,*curvature_interpolation,*secondary_interpolation;
    T_INTERPOLATION_VECTOR* normal_interpolation;
    T max_time_step;
    bool curvature_motion;
    T sigma;
public:
    bool refine_fmm_initialization_with_iterative_solver;
    int fmm_initialization_iterations;
    T fmm_initialization_iterative_tolerance;
    T fmm_initialization_iterative_drift_fraction;

    T_BOUNDARY_SCALAR* boundary;
    LEVELSET_CALLBACKS<T_GRID>* levelset_callbacks;
    T_GRID_BASED_COLLISION_GEOMETRY* collision_body_list;
    const T_FACE_ARRAYS_BOOL* face_velocities_valid_mask_current;
    bool collision_aware_signed_distance,clamp_phi_with_collision_bodies;
//protected:
    T_BOUNDARY_SCALAR& boundary_default;
    static T_LINEAR_INTERPOLATION_SCALAR interpolation_default;
    T_INTERPOLATION_SCALAR *collision_aware_interpolation_plus,*collision_aware_interpolation_minus,*collision_unaware_interpolation;
    static typename REBIND<T_LINEAR_INTERPOLATION_SCALAR,TV>::TYPE normal_interpolation_default;
    LEVELSET_NORMAL_COMPUTATION<T_GRID>* custom_normal_computation;
    T collidable_phi_replacement_value;
public:
    T_ARRAYS_BOOL valid_mask_current;
    T_ARRAYS_BOOL valid_mask_next;
protected:

    LEVELSET()
        :levelset_callbacks(0),collision_body_list(0),face_velocities_valid_mask_current(0),clamp_phi_with_collision_bodies(true),boundary_default(*new T_BOUNDARY_SCALAR),
        collision_aware_interpolation_plus(0),collision_aware_interpolation_minus(0),collision_unaware_interpolation(0),custom_normal_computation(0),
        collidable_phi_replacement_value((T)1e-5)
    {
        Set_Small_Number();
        Set_Max_Time_Step();
        curvature_motion=false; // default is no curvature motion
        Initialize_FMM_Initialization_Iterative_Solver();
        boundary=&boundary_default;
        interpolation=&interpolation_default;
        curvature_interpolation=&interpolation_default;
        normal_interpolation=&normal_interpolation_default;
        secondary_interpolation=0;
    }

    virtual ~LEVELSET()
    {
        assert(!collision_unaware_interpolation);delete collision_aware_interpolation_plus;delete collision_aware_interpolation_minus;
        delete &boundary_default;
    }

public:
    void Set_Small_Number(const T small_number_input=1e-8)
    {small_number=small_number_input;}

    void Set_Max_Time_Step(const T max_time_step_input=1e8)
    {max_time_step=max_time_step_input;}

    void Set_Curvature_Motion(const T sigma_input=1) // times by the grid spacing
    {curvature_motion=true;sigma=sigma_input;assert(sigma >= 0);}

    void Initialize_FMM_Initialization_Iterative_Solver(const bool refine_fmm_initialization_with_iterative_solver_input=true,const int fmm_initialization_iterations_input=10,
        const T fmm_initialization_iterative_tolerance_input=1e-2,const T fmm_initialization_iterative_drift_fraction_input=.1)
    {refine_fmm_initialization_with_iterative_solver=refine_fmm_initialization_with_iterative_solver_input;fmm_initialization_iterations=fmm_initialization_iterations_input;
    fmm_initialization_iterative_tolerance=fmm_initialization_iterative_tolerance_input;fmm_initialization_iterative_drift_fraction=fmm_initialization_iterative_drift_fraction_input;}

    virtual void Initialize_Levelset_Grid_Values()=0; // should call version below

protected:
    void Initialize_Levelset_Grid_Values(const T_GRID& grid_input)
    {}
public:

    void Initialize_Valid_Masks(const T_GRID& grid_input)
    {
        valid_mask_current.Resize(grid_input.Cell_Indices(3),true,true,true);
        valid_mask_next.Resize(grid_input.Cell_Indices(3),false);
    }

    void Set_Custom_Boundary(T_BOUNDARY_SCALAR& boundary_input)
    {boundary=&boundary_input;}

    void Set_Custom_Interpolation(T_INTERPOLATION_SCALAR& interpolation_input)
    {interpolation=&interpolation_input;}

    void Set_Custom_Secondary_Interpolation(T_INTERPOLATION_SCALAR& interpolation_input)
    {secondary_interpolation=&interpolation_input;}

    void Set_Custom_Normal_Interpolation(T_INTERPOLATION_VECTOR& normal_interpolation_input)
    {normal_interpolation=&normal_interpolation_input;}

    void Set_Custom_Normal_Computation(LEVELSET_NORMAL_COMPUTATION<T_GRID>* normal_computation)
    {custom_normal_computation=normal_computation;}

    void Set_Custom_Curvature_Interpolation(T_INTERPOLATION_SCALAR& curvature_interpolation_input)
    {curvature_interpolation=&curvature_interpolation_input;}

    void Enable_Collision_Aware_Interpolation(const int sign)
    {PHYSBAM_ASSERT(!collision_unaware_interpolation);collision_unaware_interpolation=interpolation;
    interpolation=sign==1?collision_aware_interpolation_plus:collision_aware_interpolation_minus;}

    void Disable_Collision_Aware_Interpolation()
    {interpolation=collision_unaware_interpolation;collision_unaware_interpolation=0;}

    void Set_Levelset_Callbacks(LEVELSET_CALLBACKS<T_GRID>& levelset_callbacks_input)
    {levelset_callbacks=&levelset_callbacks_input;}

    void Set_Collision_Body_List(T_GRID_BASED_COLLISION_GEOMETRY& collision_body_list_input,const bool set_secondary_interpolation=false)
    {
        collision_body_list=&collision_body_list_input;
        delete collision_aware_interpolation_plus;delete collision_aware_interpolation_minus;
        collision_aware_interpolation_plus=new T_LINEAR_INTERPOLATION_SCALAR;
        collision_aware_interpolation_minus=new T_LINEAR_INTERPOLATION_COLLIDABLE_CELL_SCALAR(*collision_body_list,&valid_mask_current,collidable_phi_replacement_value);
        if(set_secondary_interpolation) secondary_interpolation=collision_aware_interpolation_minus;
    }

    void Set_Face_Velocities_Valid_Mask(const T_FACE_ARRAYS_BOOL* face_velocities_valid_mask_current_input)
    {
        face_velocities_valid_mask_current=face_velocities_valid_mask_current_input;
    }
};
}
#endif
