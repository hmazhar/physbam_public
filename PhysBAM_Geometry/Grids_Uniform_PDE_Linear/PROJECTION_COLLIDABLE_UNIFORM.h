//#####################################################################
// Copyright 2002-2010, Ronald Fedkiw, Eran Guendelman, Michael Lentine, Frank Losasso, Andrew Selle, Tamar Shinar, Jerry Talton.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class PROJECTION_COLLIDABLE_UNIFORM  
//#####################################################################
#ifndef __PROJECTION_COLLIDABLE_UNIFORM__
#define __PROJECTION_COLLIDABLE_UNIFORM__

#include <PhysBAM_Tools/Data_Structures/TRIPLE.h>
#include <PhysBAM_Tools/Grids_Uniform_PDE_Linear/PROJECTION_UNIFORM.h>
#include <PhysBAM_Geometry/Grids_Uniform_PDE_Linear/POISSON_COLLIDABLE_UNIFORM.h>
#include <PhysBAM_Geometry/Interpolation_Collidable/INTERPOLATION_COLLIDABLE_POLICY.h>
namespace PhysBAM{

template<class T_GRID> class DETONATION_SHOCK_DYNAMICS;

template<class T_GRID>
class PROJECTION_COLLIDABLE_UNIFORM:public PROJECTION_UNIFORM<T_GRID>
{
    typedef typename T_GRID::VECTOR_T TV;typedef typename T_GRID::SCALAR T;
    typedef typename GRID_ARRAYS_POLICY<T_GRID>::ARRAYS_SCALAR T_ARRAYS_SCALAR;typedef typename GRID_ARRAYS_POLICY<T_GRID>::ARRAYS_BASE T_ARRAYS_BASE;typedef typename LEVELSET_POLICY<T_GRID>::LEVELSET T_LEVELSET;
    typedef typename GRID_ARRAYS_POLICY<T_GRID>::FACE_ARRAYS T_FACE_ARRAYS_SCALAR;
    typedef typename T_GRID::FACE_ITERATOR FACE_ITERATOR;typedef typename T_ARRAYS_SCALAR::template REBIND<bool>::TYPE T_ARRAYS_BOOL;
    typedef typename T_FACE_ARRAYS_SCALAR::template REBIND<bool>::TYPE T_FACE_ARRAYS_BOOL;typedef typename T_GRID::VECTOR_INT TV_INT;
public:
    typedef PROJECTION_UNIFORM<T_GRID> BASE;
    using BASE::p_grid;using BASE::p;using BASE::elliptic_solver;using BASE::poisson;using BASE::laplace;using BASE::Zero_Out_Neumann_Pocket_Velocities;
    
    LAPLACE_COLLIDABLE<T_GRID>* collidable_solver;
    LAPLACE_COLLIDABLE_UNIFORM<T_GRID>* laplace_collidable; 
    POISSON_COLLIDABLE_UNIFORM<T_GRID>* poisson_collidable;

    PROJECTION_COLLIDABLE_UNIFORM(const T_GRID& mac_grid,const bool multiphase,const bool use_poisson,const bool use_variable_beta,THREAD_QUEUE* thread_queue=0);
    PROJECTION_COLLIDABLE_UNIFORM(const T_GRID& mac_grid,T_LEVELSET& levelset_input);
    virtual ~PROJECTION_COLLIDABLE_UNIFORM();

//#####################################################################
    virtual void Initialize_Grid(const T_GRID& mac_grid) PHYSBAM_OVERRIDE;
    virtual void Apply_Pressure(T_FACE_ARRAYS_SCALAR& face_velocities,const T dt,const T time,bool scale_by_dt=false) PHYSBAM_OVERRIDE;
//#####################################################################
};
}
#endif
