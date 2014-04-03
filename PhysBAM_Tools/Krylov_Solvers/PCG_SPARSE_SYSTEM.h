//#####################################################################
// Copyright 2003-2008, Ron Fedkiw, Frederic Gibou, Geoffrey Irving, Frank Losasso, Craig Schroeder, Andrew Selle, Tamar Shinar.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#ifndef __PCG_SPARSE_SYSTEM__
#define __PCG_SPARSE_SYSTEM__
#include <PhysBAM_Tools/Krylov_Solvers/KRYLOV_SYSTEM_BASE.h>
#include <PhysBAM_Tools/Krylov_Solvers/KRYLOV_VECTOR_WRAPPER.h>
#include <PhysBAM_Tools/Vectors/Dot_Product.h>
namespace PhysBAM{

template<class T> class PCG_SPARSE;
template<class T> class SPARSE_MATRIX_FLAT_NXN;
template<class T> class VECTOR_ND;

//#####################################################################
// Class PCG_SPARSE_SYSTEM
//#####################################################################
template<class T>
class PCG_SPARSE_SYSTEM:public KRYLOV_SYSTEM_BASE<T>
{
    typedef VECTOR_ND<T> VECTOR_T;typedef SPARSE_MATRIX_FLAT_NXN<T> T_MATRIX;typedef KRYLOV_VECTOR_WRAPPER<T,VECTOR_T&> KRYLOV_VECTOR_T;

    const T_MATRIX& A;
    const bool enforce_compatibility,remove_null_space_solution_component;
    const double one_over_n;
    mutable VECTOR_T temp;
public:

    PCG_SPARSE_SYSTEM(const PCG_SPARSE<T>& pcg,const T_MATRIX& A_matrix)
        :KRYLOV_SYSTEM_BASE<T>(false,true),A(A_matrix),enforce_compatibility(pcg.enforce_compatibility),
        remove_null_space_solution_component(pcg.remove_null_space_solution_component),one_over_n(A.n?(double)1/A.n:0.0)
    {}

    void Multiply(const KRYLOV_VECTOR_BASE<T>& bx,KRYLOV_VECTOR_BASE<T>& bresult) const PHYSBAM_OVERRIDE
    {const KRYLOV_VECTOR_T& x=debug_cast<const KRYLOV_VECTOR_T&>(bx);KRYLOV_VECTOR_T& result=debug_cast<KRYLOV_VECTOR_T&>(bresult);result.v.Resize(A.n);A.Times(x.v,result.v);}

    double Inner_Product(const KRYLOV_VECTOR_BASE<T>& bx,const KRYLOV_VECTOR_BASE<T>& by) const PHYSBAM_OVERRIDE
    {const KRYLOV_VECTOR_T& x=debug_cast<const KRYLOV_VECTOR_T&>(bx),&y=debug_cast<const KRYLOV_VECTOR_T&>(by);return Dot_Product_Double_Precision(x.v,y.v);}

    T Convergence_Norm(const KRYLOV_VECTOR_BASE<T>& bx) const PHYSBAM_OVERRIDE
    {const KRYLOV_VECTOR_T& x=debug_cast<const KRYLOV_VECTOR_T&>(bx);return x.v.Maximum_Magnitude();}

    void Project_Nullspace(KRYLOV_VECTOR_BASE<T>& bx) const PHYSBAM_OVERRIDE
    {KRYLOV_VECTOR_T& x=debug_cast<KRYLOV_VECTOR_T&>(bx);if(enforce_compatibility) x.v-=(T)(x.v.Sum_Double_Precision()*one_over_n);}

    void Set_Boundary_Conditions(KRYLOV_VECTOR_BASE<T>& x) const PHYSBAM_OVERRIDE
    {if(remove_null_space_solution_component) Project_Nullspace(x);}

    void Apply_Preconditioner(const KRYLOV_VECTOR_BASE<T>& br,KRYLOV_VECTOR_BASE<T>& bz) const PHYSBAM_OVERRIDE // solve Mz=r
    {const KRYLOV_VECTOR_T& r=debug_cast<const KRYLOV_VECTOR_T&>(br);KRYLOV_VECTOR_T& z=debug_cast<KRYLOV_VECTOR_T&>(bz);
    temp.Resize(A.n);z.v.Resize(A.n);
    A.C->Solve_Forward_Substitution(r.v,temp,true); // diagonal should be treated as the identity
    A.C->Solve_Backward_Substitution(temp,z.v,false,true);} // diagonal is inverted to save on divides

    void Project(KRYLOV_VECTOR_BASE<T>& x) const PHYSBAM_OVERRIDE {}
//#####################################################################
};
}
#endif
