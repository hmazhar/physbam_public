//#####################################################################
// Copyright 2006-2008, Geoffrey Irving.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class LANCZOS_ITERATION
//#####################################################################
#ifndef __LANCZOS_ITERATION__
#define __LANCZOS_ITERATION__

#include <PhysBAM_Tools/Krylov_Solvers/KRYLOV_SYSTEM_BASE.h>
#include <PhysBAM_Tools/Matrices/BANDED_SYMMETRIC_MATRIX.h>
namespace PhysBAM{

// See Golub and Van Loan, 9.2.1, p. 480 for details.
// This version uses an extra vector q for simplicity.
// Lanczos sounds something like Lansosz (in Polish)

// cost per iteration = 1 matrix multiply/project, 2 inner products, 1 convergence norm, 1 vector add, 2 saxpy's
// approximate total flops = 2v + 9n

template<class T>
class LANCZOS_ITERATION
{
public:
    bool print_diagnostics,print_residuals;

    BANDED_SYMMETRIC_MATRIX<T,3> tridiagonal;

    LANCZOS_ITERATION()
        :print_diagnostics(true),print_residuals(false)
    {}

    ~LANCZOS_ITERATION();

//#####################################################################
    bool Tridiagonalize(const KRYLOV_SYSTEM_BASE<T>& system,KRYLOV_VECTOR_BASE<T>& w,KRYLOV_VECTOR_BASE<T>& v,KRYLOV_VECTOR_BASE<T>& q,const T tolerance,const int max_iterations);
    void Test_Symmetry(const KRYLOV_SYSTEM_BASE<T>& system,KRYLOV_VECTOR_BASE<T>& w,KRYLOV_VECTOR_BASE<T>& v,KRYLOV_VECTOR_BASE<T>& q,const int iterations);
    static void Print_Spectral_Information(const KRYLOV_SYSTEM_BASE<T>& system,KRYLOV_VECTOR_BASE<T>& w,KRYLOV_VECTOR_BASE<T>& v,KRYLOV_VECTOR_BASE<T>& q,const T tolerance,const int max_iterations);
//#####################################################################
};
}
#endif
