//#####################################################################
// Copyright 2002-2010, Ronald Fedkiw, Eran Guendelman, Michael Lentine, Frank Losasso, Avi Robinson-Mosher, Andrew Selle, Tamar Shinar, Jerry Talton.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class PROJECTION_UNIFORM  
//#####################################################################
#ifndef __PROJECTION_UNIFORM__
#define __PROJECTION_UNIFORM__

#include <PhysBAM_Tools/Data_Structures/TRIPLE.h>
#include <PhysBAM_Tools/Grids_PDE_Linear/PROJECTION.h>
#include <PhysBAM_Tools/Grids_Uniform_PDE_Linear/POISSON_UNIFORM.h>
namespace PhysBAM{

template<class T_GRID>
class PROJECTION_UNIFORM:public PROJECTION<typename T_GRID::SCALAR>
{
    typedef typename T_GRID::VECTOR_T TV;typedef typename T_GRID::SCALAR T;
    typedef typename GRID_ARRAYS_POLICY<T_GRID>::ARRAYS_SCALAR T_ARRAYS_SCALAR;typedef typename GRID_ARRAYS_POLICY<T_GRID>::ARRAYS_BASE T_ARRAYS_BASE;
    typedef typename GRID_ARRAYS_POLICY<T_GRID>::FACE_ARRAYS T_FACE_ARRAYS_SCALAR;typedef ARRAY<T,SIDED_FACE_INDEX<T_GRID::dimension> > T_FACE_ARRAYS_SLIP_SCALAR;
    typedef typename T_GRID::FACE_ITERATOR FACE_ITERATOR;typedef typename T_GRID::CELL_ITERATOR CELL_ITERATOR;typedef typename T_ARRAYS_SCALAR::template REBIND<bool>::TYPE T_ARRAYS_BOOL;
    typedef typename T_FACE_ARRAYS_SCALAR::template REBIND<bool>::TYPE T_FACE_ARRAYS_BOOL;typedef typename T_GRID::VECTOR_INT TV_INT;
    typedef typename INTERPOLATION_POLICY<T_GRID>::FACE_LOOKUP T_FACE_LOOKUP;
public:
    typedef PROJECTION<T> BASE;
    using PROJECTION<T>::use_non_zero_divergence;
    
    T_GRID p_grid; // p_grid is a cell centered MAC grid
    T_ARRAYS_SCALAR p;
    T_ARRAYS_SCALAR p_save_for_projection;
    T_FACE_ARRAYS_SCALAR face_velocities_save_for_projection;
    LAPLACE_UNIFORM<T_GRID>* elliptic_solver;
    LAPLACE_UNIFORM<T_GRID>* laplace; 
    POISSON_UNIFORM<T_GRID>* poisson;     
    T_ARRAYS_SCALAR divergence; // use this to set up a non-zero divergence
    bool use_divergence_multiplier;
    T_ARRAYS_SCALAR divergence_multiplier;
    THREAD_QUEUE* thread_queue;

    PROJECTION_UNIFORM(const T_GRID& mac_grid,const bool use_variable_beta=false,const bool use_poisson=false,THREAD_QUEUE* thread_queue=0);
protected:
    PROJECTION_UNIFORM(THREAD_QUEUE* thread_queue_input=0);
public:
    virtual ~PROJECTION_UNIFORM();

    void Use_Non_Zero_Divergence(const bool use_non_zero_divergence_input=true)
    {use_non_zero_divergence=use_non_zero_divergence_input;
    if(use_non_zero_divergence) divergence.Resize(p_grid.Domain_Indices());else divergence.Clean_Memory();}

    void Use_Divergence_Multiplier(const bool use_divergence_multiplier_input=true)
    {use_divergence_multiplier=use_divergence_multiplier_input;
    if(use_divergence_multiplier) divergence_multiplier.Resize(p_grid.Domain_Indices(3));else divergence_multiplier.Clean_Memory();}

//#####################################################################
    virtual void Initialize_Grid(const T_GRID& mac_grid);
    virtual void Make_Divergence_Free(T_FACE_ARRAYS_SCALAR& face_velocities,const T dt,const T time);
    virtual void Calculate_Kinetic_Energy_Error(T_FACE_ARRAYS_SCALAR& face_velocities,ARRAY<TV,TV_INT>& kinetic_energy_error) {PHYSBAM_FATAL_ERROR();}
    void Zero_Out_Neumann_Pocket_Velocities(T_FACE_ARRAYS_SCALAR& face_velocities);
    virtual void Apply_Pressure(T_FACE_ARRAYS_SCALAR& face_velocities,const T dt,const T time,bool scale_by_dt=false);
    void Enforce_Velocity_Compatibility(T_FACE_ARRAYS_SCALAR& face_velocities);
    void Set_Up_For_Projection(T_FACE_ARRAYS_SCALAR& face_velocities);
    void Restore_After_Projection(T_FACE_ARRAYS_SCALAR& face_velocities);
    void Exchange_Pressures_For_Projection();
    void Compute_Divergence(const T_FACE_LOOKUP& face_lookup,LAPLACE_UNIFORM<T_GRID>* solver);
    void Compute_Divergence_Threaded(RANGE<TV_INT>& domain,const T_FACE_LOOKUP& face_lookup,LAPLACE_UNIFORM<T_GRID>* solver);
//#####################################################################
};
}
#endif
