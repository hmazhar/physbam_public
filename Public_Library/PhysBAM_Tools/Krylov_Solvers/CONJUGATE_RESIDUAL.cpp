//#####################################################################
// Copyright 2007-2008, Geoffrey Irving, Avi Robinson-Mosher, Craig Schroeder, Tamar Shinar.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class CONJUGATE_RESIDUAL
//#####################################################################
#include <PhysBAM_Tools/Krylov_Solvers/CONJUGATE_RESIDUAL.h>
#include <PhysBAM_Tools/Krylov_Solvers/KRYLOV_VECTOR_BASE.h>
#include <PhysBAM_Tools/Log/LOG.h>
#include <cfloat>
#include <cmath>
#include <limits>
using namespace PhysBAM;
using ::std::abs;
//#####################################################################
// Destructor
//#####################################################################
template<class T> CONJUGATE_RESIDUAL<T>::
~CONJUGATE_RESIDUAL()
{}
//#####################################################################
// Function Solve
//#####################################################################
template<class T> bool CONJUGATE_RESIDUAL<T>::
Solve(const KRYLOV_SYSTEM_BASE<T>& system,KRYLOV_VECTOR_BASE<T>& x,const KRYLOV_VECTOR_BASE<T>& b,KRYLOV_VECTOR_BASE<T>& p,KRYLOV_VECTOR_BASE<T>& ap,KRYLOV_VECTOR_BASE<T>& ar,
    KRYLOV_VECTOR_BASE<T>& r,KRYLOV_VECTOR_BASE<T>& z,const T tolerance,const int min_iterations,const int max_iterations)
{
    // TODO: there used to a note here claiming this routine assumed x=0, but it seems to work fine in that case as long as x doesn't have a component in the nullspace
    static const T small_number=std::numeric_limits<T>::epsilon();
    system.Set_Boundary_Conditions(x);

    T rho_old=0;T convergence_norm=0;
    int iterations;for(iterations=0;;iterations++){
        bool restart=!iterations || (restart_iterations && iterations%restart_iterations==0);
        if(restart){
#ifndef COMPILE_WITHOUT_READ_WRITE_SUPPORT
            if(print_residuals) LOG::cout<<"restarting conjugate_residual"<<std::endl;
#endif
            r=b;system.Multiply(x,p);r-=p;system.Project(r);}
        // stopping conditions
        convergence_norm=system.Convergence_Norm(r);
#ifndef COMPILE_WITHOUT_READ_WRITE_SUPPORT
        if(print_residuals) LOG::cout<<convergence_norm<<std::endl;
#endif
        residual_magnitude_squared=(T)system.Inner_Product(r,r);
        nullspace_measure=(residual_magnitude_squared>small_number*small_number*100)?abs(rho_old/residual_magnitude_squared):0;
        if((convergence_norm<=tolerance || (iterations && nullspace_measure<=nullspace_tolerance)) &&
            (iterations>=min_iterations || convergence_norm<small_number)){ // TODO: get the stopping criterion right
            if(print_diagnostics) LOG::Stat("conjugate_residual iterations",iterations);if(iterations_used) *iterations_used=iterations;return true;}
        if(iterations==max_iterations) break;
        // actual iteration
        const KRYLOV_VECTOR_BASE<T>& mr=system.Precondition(r,z);
        system.Multiply(mr,ar);
        system.Project(ar);
        T rho=(T)system.Inner_Product(mr,ar);
        if(!rho) break;
        if(restart){p=mr;ap=ar;}
        else{T beta=rho/rho_old;p.Copy(beta,p,mr);ap.Copy(beta,ap,ar);}
        const KRYLOV_VECTOR_BASE<T>& map=system.Precondition(ap,z);
        T denom=(T)system.Inner_Product(map,ap);
        if(!denom) break;
        T alpha=rho/denom;
        x.Copy(alpha,p,x);
        r.Copy(-alpha,ap,r);
        rho_old=rho;}

#ifndef COMPILE_WITHOUT_READ_WRITE_SUPPORT
    if(print_diagnostics) LOG::Stat("conjugate_residual iterations",iterations);
    if(print_diagnostics) LOG::cout<<"conjugate_residual not converged after "<<max_iterations<<" iterations, error = "<<convergence_norm<<std::endl;
#endif
    if(iterations_used) *iterations_used=iterations;
    return false;
}
//#####################################################################
template class CONJUGATE_RESIDUAL<float>;
#ifndef COMPILE_WITHOUT_DOUBLE_SUPPORT
template class CONJUGATE_RESIDUAL<double>;
#endif
