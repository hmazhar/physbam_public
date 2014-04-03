//#####################################################################
// Copyright 2007-2008, Geoffrey Irving, Avi Robinson-Mosher, Craig Schroeder, Tamar Shinar.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class SYMMQMR
//#####################################################################
#include <PhysBAM_Tools/Krylov_Solvers/KRYLOV_VECTOR_BASE.h>
#include <PhysBAM_Tools/Krylov_Solvers/SYMMQMR.h>
#include <PhysBAM_Tools/Log/LOG.h>
#include <PhysBAM_Tools/Math_Tools/sqr.h>
#include <cfloat>
#include <cmath>
#include <limits>
using namespace PhysBAM;
using ::std::abs;
//#####################################################################
// Destructor
//#####################################################################
template<class T> SYMMQMR<T>::
~SYMMQMR()
{}
//#####################################################################
// Function Solve
//#####################################################################
template<class T> bool SYMMQMR<T>::
Solve(const KRYLOV_SYSTEM_BASE<T>& system,KRYLOV_VECTOR_BASE<T>& x,const KRYLOV_VECTOR_BASE<T>& b,KRYLOV_VECTOR_BASE<T>& p,KRYLOV_VECTOR_BASE<T>& ap,KRYLOV_VECTOR_BASE<T>& d,
    KRYLOV_VECTOR_BASE<T>& r,KRYLOV_VECTOR_BASE<T>& z,const T tolerance,const int min_iterations,const int max_iterations)
{
    static const T small_number=std::numeric_limits<T>::epsilon();
    system.Set_Boundary_Conditions(x);

    T rho_old=0;T convergence_norm=0;T tau_old=0;T nu_old=0;
    int iterations;for(iterations=0;;iterations++){
        bool restart=!iterations || (restart_iterations && iterations%restart_iterations==0);
        if(restart){
#ifndef COMPILE_WITHOUT_READ_WRITE_SUPPORT
            if(print_residuals) LOG::cout<<"restarting symmqmr"<<std::endl;
#endif
            r=b;system.Multiply(x,p);r-=p;system.Project(r);
            tau_old=sqrt((T)system.Inner_Product(r,r));
            p=system.Precondition(r,z);
            nu_old=(T)0;
            rho_old=(T)system.Inner_Product(r,p);}

        // actual iteration
        system.Multiply(p,ap);
        system.Project(ap);
        T theta=(T)system.Inner_Product(p,ap);
        if(!theta){if(print_diagnostics) LOG::Stat("symmqmr iterations",iterations);if(iterations_used) *iterations_used=iterations;return true;}

        T alpha=rho_old/theta;
        r.Copy(-alpha,ap,r);
        T nu=sqrt((T)system.Inner_Product(r,r))/tau_old;
        T c=1/sqrt(1+sqr(nu));
        T tau=tau_old*nu*c;
        if(!restart){d*=sqr(c)*sqr(nu_old);d.Copy(sqr(c)*alpha,p,d);}
        else d.Copy(sqr(c)*alpha,p);
        x.Copy((T)1,d,x);

        // stopping conditions
        convergence_norm=(T)system.Convergence_Norm(r);
#ifndef COMPILE_WITHOUT_READ_WRITE_SUPPORT
        if(print_residuals) LOG::cout<<convergence_norm<<std::endl;
#endif
        // since this handles indefinite matrices, using nullspace criterion from conjugate residual
        residual_magnitude_squared=(T)system.Inner_Product(r,r);
        nullspace_measure=residual_magnitude_squared?abs(rho_old/residual_magnitude_squared):0;
        if((!rho_old || convergence_norm<=tolerance || (iterations && nullspace_measure<=nullspace_tolerance)) &&
            (iterations>=min_iterations || convergence_norm<small_number)){ // TODO: get the stopping criterion right
            if(print_diagnostics) LOG::Stat("symmqmr iterations",iterations);if(iterations_used) *iterations_used=iterations;return true;}
        if(iterations==max_iterations) break;

        const KRYLOV_VECTOR_BASE<T>& mr=system.Precondition(r,z);
        T rho=(T)system.Inner_Product(r,mr);
        T beta=rho/rho_old;
        p.Copy(beta,p,mr);
        rho_old=rho;tau_old=tau;nu_old=nu;}

#ifndef COMPILE_WITHOUT_READ_WRITE_SUPPORT
    if(print_diagnostics) LOG::Stat("symmqmr iterations",iterations);if(iterations_used) *iterations_used=iterations;
    if(print_diagnostics) LOG::cout<<"symmqmr not converged after "<<max_iterations<<" iterations, error = "<<convergence_norm<<std::endl;
#endif
    return false;
}
//#####################################################################
template class SYMMQMR<float>;
#ifndef COMPILE_WITHOUT_DOUBLE_SUPPORT
template class SYMMQMR<double>;
#endif
