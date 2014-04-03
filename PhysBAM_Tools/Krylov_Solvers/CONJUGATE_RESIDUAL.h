//#####################################################################
// Copyright 2007-2008, Geoffrey Irving, Avi Robinson-Mosher, Craig Schroeder, Tamar Shinar.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class CONJUGATE_RESIDUAL
//#####################################################################
#ifndef __CONJUGATE_RESIDUAL__
#define __CONJUGATE_RESIDUAL__

#include <PhysBAM_Tools/Krylov_Solvers/KRYLOV_SOLVER.h>
namespace PhysBAM{

// for details see Choi, Sou-Cheng, "Iterative methods for singular linear equations and least-squares problems", p. 27 (Table 2.12, Algorithm CR)
// (the variables ar,ap,rho here correspond to z,w,mu there)

// cost per iteration = 1 matrix multiply/project, 3 inner products, 1 convergence norm, 4 saxpy's
// approximate total flops = 2v + 14n (TODO: this is 5n more than it should be according to the above reference)

template<class T>
class CONJUGATE_RESIDUAL:public KRYLOV_SOLVER<T>
{
    typedef KRYLOV_SOLVER<T> BASE;
public:
    using BASE::restart_iterations;using BASE::residual_magnitude_squared;using BASE::iterations_used;using BASE::print_diagnostics;using BASE::print_residuals;
    using BASE::nullspace_measure;using BASE::nullspace_tolerance;using BASE::Solve;

    CONJUGATE_RESIDUAL()
    {}

    virtual ~CONJUGATE_RESIDUAL();

//#####################################################################
    bool Solve(const KRYLOV_SYSTEM_BASE<T>& system,KRYLOV_VECTOR_BASE<T>& x,const KRYLOV_VECTOR_BASE<T>& b,KRYLOV_VECTOR_BASE<T>& p,KRYLOV_VECTOR_BASE<T>& ap,KRYLOV_VECTOR_BASE<T>& ar,
        KRYLOV_VECTOR_BASE<T>& r,KRYLOV_VECTOR_BASE<T>& z,const T tolerance,const int min_iterations,const int max_iterations) PHYSBAM_OVERRIDE;
//#####################################################################
};
}
#endif
