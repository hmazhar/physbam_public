//#####################################################################
// Copyright 2009, Geoffrey Irving, Tamar Shinar.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#include <PhysBAM_Tools/Krylov_Solvers/KRYLOV_VECTOR_WRAPPER.h>
#include <PhysBAM_Tools/Krylov_Solvers/MATRIX_SYSTEM.h>
#include <PhysBAM_Tools/Matrices/SPARSE_MATRIX_FLAT_MXN.h>
using namespace PhysBAM;
using namespace KRYLOV;
//#####################################################################
// Constructor
//#####################################################################
template<class T_MATRIX,class T,class VECTOR_T,class T_MATRIX_PRECON> MATRIX_SYSTEM<T_MATRIX,T,VECTOR_T,T_MATRIX_PRECON>::
MATRIX_SYSTEM(const T_MATRIX& A_input)
    :KRYLOV_SYSTEM_BASE<T>(false,true),A(A_input)
{
}
//#####################################################################
// Destructor
//#####################################################################
template<class T_MATRIX,class T,class VECTOR_T,class T_MATRIX_PRECON> MATRIX_SYSTEM<T_MATRIX,T,VECTOR_T,T_MATRIX_PRECON>::
~MATRIX_SYSTEM()
{
}
//#####################################################################
// Function Multiply
//#####################################################################
template<class T_MATRIX,class T,class VECTOR_T,class T_MATRIX_PRECON> void MATRIX_SYSTEM<T_MATRIX,T,VECTOR_T,T_MATRIX_PRECON>::
Multiply(const KRYLOV_VECTOR_BASE<T>& x,KRYLOV_VECTOR_BASE<T>& result) const PHYSBAM_OVERRIDE
{
    const VECTOR_T& vx=dynamic_cast<const VECTOR_T&>(x);
    VECTOR_T& vresult=dynamic_cast<VECTOR_T&>(result);
    A.Times(vx.v,vresult.v);
}
//#####################################################################
// Function Inner_Product
//#####################################################################
template<class T_MATRIX,class T,class VECTOR_T,class T_MATRIX_PRECON> double MATRIX_SYSTEM<T_MATRIX,T,VECTOR_T,T_MATRIX_PRECON>::
Inner_Product(const KRYLOV_VECTOR_BASE<T>& x,const KRYLOV_VECTOR_BASE<T>& y) const PHYSBAM_OVERRIDE
{
    const VECTOR_T& vx=dynamic_cast<const VECTOR_T&>(x),vy=dynamic_cast<const VECTOR_T&>(y);
    return vx.v.Dot_Product_Double_Precision(vx.v,vy.v);
}
//#####################################################################
// Function Convergence_Norm
//#####################################################################
template<class T_MATRIX,class T,class VECTOR_T,class T_MATRIX_PRECON> T MATRIX_SYSTEM<T_MATRIX,T,VECTOR_T,T_MATRIX_PRECON>::
Convergence_Norm(const KRYLOV_VECTOR_BASE<T>& x) const PHYSBAM_OVERRIDE
{
    const VECTOR_T& vx=dynamic_cast<const VECTOR_T&>(x);
    return vx.v.Maximum_Magnitude();
}
//#####################################################################
// Function Project
//#####################################################################
template<class T_MATRIX,class T,class VECTOR_T,class T_MATRIX_PRECON> void MATRIX_SYSTEM<T_MATRIX,T,VECTOR_T,T_MATRIX_PRECON>::
Project(KRYLOV_VECTOR_BASE<T>& x) const PHYSBAM_OVERRIDE
{
}
//#####################################################################
// Function Set_Boundary_Conditions
//#####################################################################
template<class T_MATRIX,class T,class VECTOR_T,class T_MATRIX_PRECON> void MATRIX_SYSTEM<T_MATRIX,T,VECTOR_T,T_MATRIX_PRECON>::
Set_Boundary_Conditions(KRYLOV_VECTOR_BASE<T>& x) const PHYSBAM_OVERRIDE
{
}
//#####################################################################
// Function Project_Nullspace
//#####################################################################
template<class T_MATRIX,class T,class VECTOR_T,class T_MATRIX_PRECON> void MATRIX_SYSTEM<T_MATRIX,T,VECTOR_T,T_MATRIX_PRECON>::
Project_Nullspace(KRYLOV_VECTOR_BASE<T>& x) const PHYSBAM_OVERRIDE
{
} 
//#####################################################################
// Function Initialize_Preconditioner
//#####################################################################
template<class T_MATRIX,class T,class VECTOR_T,class T_MATRIX_PRECON> void MATRIX_SYSTEM<T_MATRIX,T,VECTOR_T,T_MATRIX_PRECON>::
Set_Preconditioner(const T_MATRIX_PRECON& preconditioner,VECTOR_T& vector)
{
    this->use_preconditioner=true;
    P=&preconditioner;
    temp_vector=&vector;
}
//#####################################################################
// Function Apply_Preconditioner
//#####################################################################
template<class T_MATRIX,class T,class VECTOR_T,class T_MATRIX_PRECON> void MATRIX_SYSTEM<T_MATRIX,T,VECTOR_T,T_MATRIX_PRECON>::
Apply_Preconditioner(const KRYLOV_VECTOR_BASE<T>& r,KRYLOV_VECTOR_BASE<T>& z) const
{
    const VECTOR_T& vr=dynamic_cast<const VECTOR_T&>(r);
    VECTOR_T& vz=dynamic_cast<VECTOR_T&>(z);
    P->Solve_Forward_Substitution(vr.v,temp_vector->v,true);
    P->Solve_Backward_Substitution(temp_vector->v,vz.v,false,true);
}
template class MATRIX_SYSTEM<SPARSE_MATRIX_FLAT_NXN<float>,float,KRYLOV_VECTOR_WRAPPER<float,VECTOR_ND<float> > >;
#ifndef COMPILE_WITHOUT_DOUBLE_SUPPORT
template class MATRIX_SYSTEM<SPARSE_MATRIX_FLAT_NXN<double>,double,KRYLOV_VECTOR_WRAPPER<double,VECTOR_ND<double> > >;
#endif
