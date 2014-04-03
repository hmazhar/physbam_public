//#####################################################################
// Copyright 2006-2008, Geoffrey Irving.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class LANCZOS_ITERATION
//#####################################################################
#include <PhysBAM_Tools/Krylov_Solvers/KRYLOV_VECTOR_BASE.h>
#include <PhysBAM_Tools/Krylov_Solvers/LANCZOS_ITERATION.h>
#include <PhysBAM_Tools/Log/LOG.h>
#include <cfloat>
#include <limits>
using namespace PhysBAM;
//#####################################################################
// Destructor
//#####################################################################
template<class T> LANCZOS_ITERATION<T>::
~LANCZOS_ITERATION()
{}
//#####################################################################
// Function Tridiagonalize
//#####################################################################
template<class T> bool LANCZOS_ITERATION<T>::
Tridiagonalize(const KRYLOV_SYSTEM_BASE<T>& system,KRYLOV_VECTOR_BASE<T>& w,KRYLOV_VECTOR_BASE<T>& v,KRYLOV_VECTOR_BASE<T>& q,const T tolerance,const int max_iterations)
{
    PHYSBAM_ASSERT(!system.use_preconditioner);

    tridiagonal.A.Remove_All();

    system.Project(w);
    T magnitude_squared=(T)system.Inner_Product(w,w);
    if(!magnitude_squared){
#ifndef COMPILE_WITHOUT_READ_WRITE_SUPPORT
        if(print_diagnostics) LOG::cout<<"lanczos iterations = "<<0<<std::endl;
#endif
        return true;}
    w*=1/sqrt(magnitude_squared);

    T convergence_norm=0;
    for(int iterations=1;;iterations++){
        system.Multiply(w,q);
        if(iterations==1) v=q;else v+=q;
        system.Project(v);
        T alpha=(T)system.Inner_Product(w,v);
        v.Copy(-alpha,w,v);
        T beta=sqrt((T)system.Inner_Product(v,v));
        tridiagonal.A.Append(VECTOR<T,2>(alpha,beta));
        convergence_norm=system.Convergence_Norm(v);
#ifndef COMPILE_WITHOUT_READ_WRITE_SUPPORT
        if(print_residuals) LOG::cout<<convergence_norm<<std::endl;
#endif
        if(convergence_norm<=tolerance){
#ifndef COMPILE_WITHOUT_READ_WRITE_SUPPORT
            if(print_diagnostics) LOG::cout<<"lanczos iterations = "<<iterations<<std::endl;
#endif
            return true;}
        if(iterations==max_iterations) break;
        q=w;
        w.Copy(1/beta,v);
        v.Copy(-beta,q);}

#ifndef COMPILE_WITHOUT_READ_WRITE_SUPPORT
    if(print_diagnostics) LOG::cout<<"lanczos not converged after "<<max_iterations<<" iterations, error = "<<convergence_norm<<std::endl;
#endif
    return false;
}
//#####################################################################
// Function Test_Symmetry
//#####################################################################
template<class T> void LANCZOS_ITERATION<T>::
Test_Symmetry(const KRYLOV_SYSTEM_BASE<T>& system,KRYLOV_VECTOR_BASE<T>& w,KRYLOV_VECTOR_BASE<T>& v,KRYLOV_VECTOR_BASE<T>& q,const int iterations)
{
    PHYSBAM_ASSERT(!system.use_preconditioner);

    VECTOR<KRYLOV_VECTOR_BASE<T>*,3> V(&w,&v,&q);
    system.Project(*V.x);
    system.Multiply(*V.x,*V.y);
    system.Project(*V.y);

    for(int i=1;i<=iterations;i++){
        double Ax_dot_Ax=system.Inner_Product(*V.y,*V.y);
        if(!Ax_dot_Ax){
#ifndef COMPILE_WITHOUT_READ_WRITE_SUPPORT
            LOG::cout<<"Ax = 0, exiting"<<std::endl;
#endif
            break;}
        system.Multiply(*V.y,*V.z);
        system.Project(*V.z);
#ifndef COMPILE_WITHOUT_READ_WRITE_SUPPORT
        double x_dot_AAx=system.Inner_Product(*V.x,*V.z);
        LOG::cout<<"Ax_dot_Ax = "<<Ax_dot_Ax<<", x_dot_AAx = "<<x_dot_AAx<<", relative_error = "<<abs(Ax_dot_Ax-x_dot_AAx)/Ax_dot_Ax<<std::endl;
#endif
        V=Vector(V.y,V.z,V.x);
        // rescale to avoid overflow
        T scale=1/sqrt(Ax_dot_Ax);
        *V.x*=scale;
        *V.y*=scale;}
}
//#####################################################################
// Function Print_Spectral_Information
//#####################################################################
template<class T> void LANCZOS_ITERATION<T>::
Print_Spectral_Information(const KRYLOV_SYSTEM_BASE<T>& system,KRYLOV_VECTOR_BASE<T>& w,KRYLOV_VECTOR_BASE<T>& v,KRYLOV_VECTOR_BASE<T>& q,const T tolerance,const int max_iterations)
{
    LANCZOS_ITERATION lanczos;lanczos.Tridiagonalize(system,w,v,q,tolerance,max_iterations);
    lanczos.tridiagonal.Print_Spectral_Information();
}
//#####################################################################
template class LANCZOS_ITERATION<float>;
#ifndef COMPILE_WITHOUT_DOUBLE_SUPPORT
template class LANCZOS_ITERATION<double>;
#endif
