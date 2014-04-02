//#####################################################################
// Copyright 2008, Geoffrey Irving.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#include <PhysBAM_Tools/Krylov_Solvers/KRYLOV_SYSTEM_BASE.h>
#include <PhysBAM_Tools/Krylov_Solvers/KRYLOV_VECTOR_BASE.h>
#include <PhysBAM_Tools/Log/DEBUG_UTILITIES.h>
#include <PhysBAM_Tools/Log/LOG.h>
#include <PhysBAM_Tools/Math_Tools/maxabs.h>
#include <PhysBAM_Tools/Random_Numbers/RANDOM_NUMBERS.h>
#include <PhysBAM_Tools/Utilities/DEBUG_CAST.h>
using namespace PhysBAM;
template<class T> KRYLOV_SYSTEM_BASE<T>::
KRYLOV_SYSTEM_BASE(const bool use_preconditioner,const bool preconditioner_commutes_with_projection)
    :use_preconditioner(use_preconditioner),preconditioner_commutes_with_projection(preconditioner_commutes_with_projection)
{
}
template<class T> KRYLOV_SYSTEM_BASE<T>::
~KRYLOV_SYSTEM_BASE()
{
}
template<class T> void KRYLOV_SYSTEM_BASE<T>::
Test_System(KRYLOV_VECTOR_BASE<T>& x,KRYLOV_VECTOR_BASE<T>& y,KRYLOV_VECTOR_BASE<T>& z) const
{
    T tolerance=(T)1e-5;
    RANDOM_NUMBERS<T> random;
    double a,b,r;
    bool pass=true;
    int n=x.Raw_Size();
    for(int i=1;i<=n;i++){
        x.Raw_Get(i)=random.Get_Uniform_Number(-1,1);
        y.Raw_Get(i)=random.Get_Uniform_Number(-1,1);}

    a=Inner_Product(x,y);
    b=Inner_Product(y,x);
    r=abs(a-b)/maxabs(1e-30,a,b);
    if(r>tolerance) {pass=false;LOG::cout<<"Inner Product Symmetry Test: "<<a<<"  vs  "<<b<<"  relative  "<<r<<std::endl;}

    Multiply(x,z);
    a=Inner_Product(y,z);
    Multiply(y,z);
    b=Inner_Product(x,z);
    r=abs(a-b)/maxabs(1e-30,a,b);
    if(r>tolerance) {pass=false;LOG::cout<<"System Symmetry Test: "<<a<<"  vs  "<<b<<"  relative  "<<abs(a-b)/maxabs(1e-30,a,b)<<std::endl;}

    z=x;
    Project(z);
    a=Inner_Product(y,z);
    z=y;
    Project(z);
    b=Inner_Product(x,z);
    r=abs(a-b)/maxabs(1e-30,a,b);
    if(r>tolerance) {pass=false;LOG::cout<<"Projection Symmetry Test: "<<a<<"  vs  "<<b<<"  relative  "<<abs(a-b)/maxabs(1e-30,a,b)<<std::endl;}

    z=x;
    Project(z);
    a=Inner_Product(z,z);
    Project(z);
    b=Inner_Product(z,z);
    r=abs(a-b)/maxabs(1e-30,a,b);
    if(r>tolerance) {pass=false;LOG::cout<<"Projection Idempotence Test: "<<a<<"  vs  "<<b<<"  relative  "<<abs(a-b)/maxabs(1e-30,a,b)<<std::endl;}

    z=x;
    Project_Nullspace(z);
    a=Inner_Product(y,z);
    z=y;
    Project_Nullspace(z);
    b=Inner_Product(x,z);
    r=abs(a-b)/maxabs(1e-30,a,b);
    if(r>tolerance) {pass=false;LOG::cout<<"Project Nullspace Symmetry Test: "<<a<<"  vs  "<<b<<"  relative  "<<abs(a-b)/maxabs(1e-30,a,b)<<std::endl;}

    z=x;
    Project_Nullspace(z);
    a=Inner_Product(z,z);
    Project_Nullspace(z);
    b=Inner_Product(z,z);
    r=abs(a-b)/maxabs(1e-30,a,b);
    if(r>tolerance) {pass=false;LOG::cout<<"Project Nullspace Idempotence Test: "<<a<<"  vs  "<<b<<"  relative  "<<abs(a-b)/maxabs(1e-30,a,b)<<std::endl;}

    z*=(T)0;
    a=Inner_Product(y,Precondition(x,z));
    z*=(T)0;
    b=Inner_Product(x,Precondition(y,z));
    r=abs(a-b)/maxabs(1e-30,a,b);
    if(r>tolerance) {pass=false;LOG::cout<<"Preconditioner Symmetry Test: "<<a<<"  vs  "<<b<<"  relative  "<<abs(a-b)/maxabs(1e-30,a,b)<<std::endl;}

    LOG::cout<<"Krylov System Test Result: "<<(pass?"PASS":"FAIL")<<std::endl;
}
template<class T> void KRYLOV_SYSTEM_BASE<T>::
Apply_Preconditioner(const KRYLOV_VECTOR_BASE<T>& r,KRYLOV_VECTOR_BASE<T>& z) const
{
    PHYSBAM_FUNCTION_IS_NOT_DEFINED();
}
template<class T> void KRYLOV_SYSTEM_BASE<T>::
Project_Nullspace(KRYLOV_VECTOR_BASE<T>& x) const
{
    PHYSBAM_FUNCTION_IS_NOT_DEFINED();
}
template<class T> const KRYLOV_VECTOR_BASE<T>& KRYLOV_SYSTEM_BASE<T>::
Precondition(const KRYLOV_VECTOR_BASE<T>& r,KRYLOV_VECTOR_BASE<T>& z) const
{
    if(!use_preconditioner) return r;
    Apply_Preconditioner(r,z);
    if(!preconditioner_commutes_with_projection){Project(z);Project_Nullspace(z);}
    return z;
}
template class KRYLOV_SYSTEM_BASE<float>;
#ifndef COMPILE_WITHOUT_DOUBLE_SUPPORT
template class KRYLOV_SYSTEM_BASE<double>;
#endif
