//#####################################################################
// Copyright 2004-2008, Ron Fedkiw, Geoffrey Irving, Michael Lentine, Craig Schroeder, Andrew Selle, Tamar Shinar, Eftychios Sifakis.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class IMPLICIT_SOLVE_PARAMETERS
//#####################################################################
#include <PhysBAM_Tools/Krylov_Solvers/IMPLICIT_SOLVE_PARAMETERS.h>
#include <PhysBAM_Tools/Vectors/VECTOR.h>
#include <climits>
namespace PhysBAM{
//#####################################################################
// Constructor
//#####################################################################
template<class TV> IMPLICIT_SOLVE_PARAMETERS<TV>::
IMPLICIT_SOLVE_PARAMETERS()
    :cg_iterations(200),cg_restart_iterations(0),cg_tolerance((T)1e-2),spectral_analysis(false),lanczos_iterations(20),cg_projection_iterations(5),
    throw_exception_on_backward_euler_failure(true),evolution_solver_type(krylov_solver_cg),project_nullspace_frequency(1),use_half_fully_implicit(false),
    test_system(false),print_matrix(false)
{
}
//#####################################################################
// Destructor
//#####################################################################
template<class TV> IMPLICIT_SOLVE_PARAMETERS<TV>::
~IMPLICIT_SOLVE_PARAMETERS()
{
}
//#####################################################################
template class IMPLICIT_SOLVE_PARAMETERS<VECTOR<float,1> >;
template class IMPLICIT_SOLVE_PARAMETERS<VECTOR<float,2> >;
template class IMPLICIT_SOLVE_PARAMETERS<VECTOR<float,3> >;
#ifndef COMPILE_WITHOUT_DOUBLE_SUPPORT
template class IMPLICIT_SOLVE_PARAMETERS<VECTOR<double,1> >;
template class IMPLICIT_SOLVE_PARAMETERS<VECTOR<double,2> >;
template class IMPLICIT_SOLVE_PARAMETERS<VECTOR<double,3> >;
#endif
}
