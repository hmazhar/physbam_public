//#####################################################################
// Copyright 2002-2007, Ronald Fedkiw, Jon Gretarsson, Eran Guendelman, Frank Losasso, Andrew Selle, Tamar Shinar.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class POISSON_UNIFORM
//#####################################################################
#ifndef __POISSON_UNIFORM__
#define __POISSON_UNIFORM__

#include <PhysBAM_Tools/Grids_PDE_Linear/POISSON.h>
#include <PhysBAM_Tools/Grids_Uniform_PDE_Linear/LAPLACE_UNIFORM.h>
namespace PhysBAM{

template<class T> class SPARSE_MATRIX_FLAT_NXN;

template<class T_GRID>
class POISSON_UNIFORM:public POISSON<typename T_GRID::SCALAR>,public LAPLACE_UNIFORM<T_GRID>
{
    typedef typename T_GRID::VECTOR_T TV;typedef typename TV::SCALAR T;
    typedef typename T_GRID::VECTOR_INT TV_INT;typedef typename T_GRID::CELL_ITERATOR CELL_ITERATOR;
    typedef typename GRID_ARRAYS_POLICY<T_GRID>::ARRAYS_SCALAR T_ARRAYS_SCALAR;typedef typename T_ARRAYS_SCALAR::template REBIND<int>::TYPE T_ARRAYS_INT;
    typedef typename GRID_ARRAYS_POLICY<T_GRID>::FACE_ARRAYS T_FACE_ARRAYS_SCALAR;typedef typename T_GRID::FACE_ITERATOR FACE_ITERATOR;
public:
    using LAPLACE_UNIFORM<T_GRID>::grid;using POISSON<T>::use_variable_beta;using POISSON<T>::beta_given_on_faces;
    using POISSON<T>::use_weighted_divergence;using POISSON<T>::multiphase;
    using LAPLACE_UNIFORM<T_GRID>::psi_N;using LAPLACE_UNIFORM<T_GRID>::periodic_boundary;
    using LAPLACE_UNIFORM<T_GRID>::filled_region_colors;
    using LAPLACE_UNIFORM<T_GRID>::filled_region_touches_dirichlet;using LAPLACE_UNIFORM<T_GRID>::solve_neumann_regions;

    T_FACE_ARRAYS_SCALAR beta_face;
    T_ARRAYS_SCALAR variable_beta;
    T_FACE_ARRAYS_SCALAR divergence_face_weights;
public:

    POISSON_UNIFORM(const T_GRID& grid_input,T_ARRAYS_SCALAR& u_input,const bool initialize_grid,const bool multiphase_input,const bool enforce_compatibility_input);
    virtual ~POISSON_UNIFORM();

    void Set_Variable_beta(const bool beta_given_on_faces_input=false)
    {use_variable_beta=true;beta_given_on_faces=beta_given_on_faces_input;
    if(!beta_given_on_faces) variable_beta.Resize(grid.Domain_Indices(1));else variable_beta.Clean_Memory();}

    void Use_Weighted_Divergence(bool use_face_velocity_weights_input=true)
    {use_weighted_divergence=use_face_velocity_weights_input;
    if(use_weighted_divergence) divergence_face_weights.Resize(grid,1);}

//#####################################################################
    void Initialize_Grid(const T_GRID& grid_input);
    void Find_Variable_beta();
    void Find_A_Part_Two(RANGE<TV_INT>& domain,ARRAY<SPARSE_MATRIX_FLAT_NXN<T> >& A_array,ARRAY<VECTOR_ND<T> >& b_array,T_ARRAYS_INT& cell_index_to_matrix_index) PHYSBAM_OVERRIDE;
//#####################################################################
};
}
#endif
