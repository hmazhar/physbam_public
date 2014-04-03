//#####################################################################
// Copyright 2003-2008, Ron Fedkiw, Frederic Gibou, Geoffrey Irving, Michael Lentine, Frank Losasso, Craig Schroeder, Andrew Selle, Tamar Shinar.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#include <PhysBAM_Tools/Krylov_Solvers/CONJUGATE_GRADIENT.h>
#include <PhysBAM_Tools/Krylov_Solvers/CONJUGATE_RESIDUAL.h>
#include <PhysBAM_Tools/Krylov_Solvers/PCG_SPARSE.h>
#include <PhysBAM_Tools/Krylov_Solvers/PCG_SPARSE_SYSTEM.h>
#include <PhysBAM_Tools/Krylov_Solvers/SYMMQMR.h>
#include <PhysBAM_Tools/Matrices/SPARSE_MATRIX_FLAT_NXN.h>
using namespace PhysBAM;
//#####################################################################
// Function Solve  
//#####################################################################
template<class T> void PCG_SPARSE<T>::
Solve(SPARSE_MATRIX_FLAT_NXN<T>& A_matrix,VECTOR_ND<T>& x,VECTOR_ND<T>& b,VECTOR_ND<T>& q,VECTOR_ND<T>& s,VECTOR_ND<T>& r,VECTOR_ND<T>& k,VECTOR_ND<T>& z,const T tolerance,const bool recompute_preconditioner)
{
    int desired_iterations=A_matrix.n-enforce_compatibility;if(maximum_iterations) desired_iterations=maximum_iterations;

    CONJUGATE_GRADIENT<T> cg;
    CONJUGATE_RESIDUAL<T> cr;
    SYMMQMR<T> symmqmr;
    KRYLOV_SOLVER<T>* solver=0;
    const char* solver_name=0;
    if(evolution_solver_type==krylov_solver_cg){solver=&cg;solver_name="CG";}
    else if(evolution_solver_type==krylov_solver_cr){solver=&cr;solver_name="CONJUGATE_RESIDUAL";}
    else if(evolution_solver_type==krylov_solver_symmqmr){solver=&symmqmr;solver_name="SYMMQMR";}
    solver->print_diagnostics=show_results;solver->print_residuals=show_residual;
    solver->restart_iterations=cg_restart_iterations;
    PCG_SPARSE_SYSTEM<T> system(*this,A_matrix);

    if(incomplete_cholesky){
        system.use_preconditioner=true;
        system.preconditioner_commutes_with_projection=false;
        if(recompute_preconditioner || !A_matrix.C)
            A_matrix.Construct_Incomplete_Cholesky_Factorization(modified_incomplete_cholesky,modified_incomplete_cholesky_coefficient,
                preconditioner_zero_tolerance,preconditioner_zero_replacement);}

    KRYLOV_VECTOR_WRAPPER<T,VECTOR_ND<T>&> kx(x),kb(b),kq(q),ks(s),kr(r),kk(k),kz(z);
    PHYSBAM_ASSERT(&kx.v==&x && &kb.v==&b && &kq.v==&q && &ks.v==&s && &kr.v==&r && &kk.v==&k && &kz.v==&z);
    solver->Solve(system,kx,kb,kq,ks,kr,kk,kz,tolerance,0,desired_iterations);
}
//#####################################################################
template class PCG_SPARSE<float>;
#ifndef COMPILE_WITHOUT_DOUBLE_SUPPORT
template class PCG_SPARSE<double>;
#endif
