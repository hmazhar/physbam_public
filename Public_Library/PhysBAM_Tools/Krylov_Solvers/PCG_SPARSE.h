//#####################################################################
// Copyright 2004-2008, Ron Fedkiw, Geoffrey Irving, Michael Lentine, Frank Losasso, Andrew Selle.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class PCG_SPARSE
//#####################################################################
//
// Solve Ax=b with Preconditioned Conjugate Gradient, Modified Incomplete Cholesky preconditioner.
// In the Neumann case, the input A, b, & x need to sum to 0 for compatibility.
// Input A as a sparse (symmetric positive definite) matrix & b as a vector. Input x as a vector. 
// Input tolerance for stopping iteration when the maximum residual is less than or equal to the tolerance, i.e. max|r|=max|b-Ax| <= tolerance. Otherwise stops at max iterations.
//               
//#####################################################################
#ifndef __PCG_SPARSE__
#define __PCG_SPARSE__

#include <PhysBAM_Tools/Krylov_Solvers/KRYLOV_SOLVER.h>
#include <PhysBAM_Tools/Utilities/NONCOPYABLE.h>
namespace PhysBAM{

template<class T> class VECTOR_ND;
template<class T> class SPARSE_MATRIX_FLAT_NXN;

template<class T>
class PCG_SPARSE:public NONCOPYABLE
{
public:
    bool enforce_compatibility; // for Neumann boundary conditions
    bool remove_null_space_solution_component; 
    bool show_residual,show_results;
    bool incomplete_cholesky; // true when using this preconditioner or the modified version below
    bool modified_incomplete_cholesky; // true when using this preconditioner
    T modified_incomplete_cholesky_coefficient;
    int maximum_iterations;
    int cg_restart_iterations;
    T preconditioner_zero_tolerance,preconditioner_zero_replacement;
    KRYLOV_SOLVER_TYPE evolution_solver_type;

    PCG_SPARSE()
        :cg_restart_iterations(0),preconditioner_zero_tolerance((T)1e-8),preconditioner_zero_replacement((T)1e-8),evolution_solver_type(krylov_solver_cg)
    {
        Enforce_Compatibility(false);
        Remove_Null_Space_Solution_Component(false);
        Show_Residuals(false);Show_Results(false);
        Use_Modified_Incomplete_Cholesky();
        Set_Maximum_Iterations();
    }

    virtual ~PCG_SPARSE()
    {}

    void Enforce_Compatibility(const bool enforce_compatibility_input=true) //  A, b, & x inputs need to sum to 0 for compatibility
    {enforce_compatibility=enforce_compatibility_input;}

    void Remove_Null_Space_Solution_Component(const bool remove_null_space_solution_component_input=true)
    {remove_null_space_solution_component=remove_null_space_solution_component_input;}
    
    void Show_Residuals(const bool show_residual_input=true)
    {show_residual=show_residual_input;}

    void Show_Results(const bool show_results_input=true)
    {show_results=show_results_input;}

    void Do_Not_Show_Results()
    {show_results=false;}

    void Set_Maximum_Iterations(const int maximum_iterations_input=0) // initially set to zero - needs to be set somewhere else!
    {maximum_iterations=maximum_iterations_input;}

    void Use_Conjugate_Gradient() // no preconditioner
    {incomplete_cholesky=false;modified_incomplete_cholesky=false;}

    void Use_Incomplete_Cholesky()
    {incomplete_cholesky=true;modified_incomplete_cholesky=false;}
    
    // use .97 for octrees and .99 for uniform grids!
    void Use_Modified_Incomplete_Cholesky(const T modified_incomplete_cholesky_coefficient_input=(T).97)
    {incomplete_cholesky=true;modified_incomplete_cholesky=true;modified_incomplete_cholesky_coefficient=modified_incomplete_cholesky_coefficient_input;} // note that both are true
    
//#####################################################################
    virtual void Solve(SPARSE_MATRIX_FLAT_NXN<T>& A_matrix,VECTOR_ND<T>& x,VECTOR_ND<T>& b,VECTOR_ND<T>& q,VECTOR_ND<T>& s,VECTOR_ND<T>& r,VECTOR_ND<T>& k,VECTOR_ND<T>& z,const T tolerance=1e-7,const bool recompute_preconditioner=true);
//#####################################################################
};
}
#endif
