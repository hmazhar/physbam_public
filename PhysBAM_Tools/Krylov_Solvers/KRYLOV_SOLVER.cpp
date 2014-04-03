//#####################################################################
// Copyright 2008, Craig Schroeder.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#include <PhysBAM_Tools/Krylov_Solvers/KRYLOV_SOLVER.h>
using namespace PhysBAM;
//#####################################################################
// Constructor
//#####################################################################
template<class T> KRYLOV_SOLVER<T>::
KRYLOV_SOLVER()
    :print_diagnostics(true),print_residuals(false),nullspace_tolerance((T)1e-5),iterations_used(0),residual_magnitude_squared(0),nullspace_measure(0),restart_iterations(0)
{
}
//#####################################################################
// Destructor
//#####################################################################
template<class T> KRYLOV_SOLVER<T>::
~KRYLOV_SOLVER()
{
}
template class KRYLOV_SOLVER<float>;
#ifndef COMPILE_WITHOUT_DOUBLE_SUPPORT
template class KRYLOV_SOLVER<double>;
#endif
