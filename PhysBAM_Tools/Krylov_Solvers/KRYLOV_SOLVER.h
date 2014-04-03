//#####################################################################
// Copyright 2008, Craig Schroeder.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class KRYLOV_SOLVER
//#####################################################################
#ifndef __KRYLOV_SOLVER__
#define __KRYLOV_SOLVER__
#include <PhysBAM_Tools/Krylov_Solvers/KRYLOV_SYSTEM_BASE.h>
#include <PhysBAM_Tools/Log/DEBUG_UTILITIES.h>
namespace PhysBAM{
namespace KRYLOV{template<class T_MATRIX> struct IS_MATRIX;}

enum KRYLOV_SOLVER_TYPE {krylov_solver_cg,krylov_solver_cr,krylov_solver_symmqmr};

template<class T>
class KRYLOV_SOLVER
{
public:
    bool print_diagnostics,print_residuals;
    T nullspace_tolerance; // don't attempt to invert eigenvalues approximately less than nullspace_tolerance*max_eigenvalue
    int* iterations_used;
    T residual_magnitude_squared,nullspace_measure; // extra convergence information
    int restart_iterations;

    KRYLOV_SOLVER();
    virtual ~KRYLOV_SOLVER();

     // version without preconditioning
    bool Solve(const KRYLOV_SYSTEM_BASE<T>& system,KRYLOV_VECTOR_BASE<T>& x,const KRYLOV_VECTOR_BASE<T>& b,KRYLOV_VECTOR_BASE<T>& p,KRYLOV_VECTOR_BASE<T>& ap,KRYLOV_VECTOR_BASE<T>& ar,
        KRYLOV_VECTOR_BASE<T>& r,const T tolerance,const int min_iterations,const int max_iterations)
    {KRYLOV_VECTOR_BASE<T>* z=0;PHYSBAM_ASSERT(!system.use_preconditioner);return Solve(system,x,b,p,ap,ar,r,*z,tolerance,min_iterations,max_iterations);}

//#####################################################################
    virtual bool Solve(const KRYLOV_SYSTEM_BASE<T>& system,KRYLOV_VECTOR_BASE<T>& x,const KRYLOV_VECTOR_BASE<T>& b,KRYLOV_VECTOR_BASE<T>& q,KRYLOV_VECTOR_BASE<T>& s,
        KRYLOV_VECTOR_BASE<T>& r,KRYLOV_VECTOR_BASE<T>& k,KRYLOV_VECTOR_BASE<T>& z,const T tolerance,const int min_iterations,const int max_iterations)=0;
//#####################################################################
};
}
#endif
