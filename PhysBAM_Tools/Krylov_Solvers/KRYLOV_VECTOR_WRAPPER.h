//#####################################################################
// Copyright 2009, Craig Schroeder.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class KRYLOV_VECTOR_WRAPPER
//#####################################################################
#ifndef __KRYLOV_VECTOR_WRAPPER__
#define __KRYLOV_VECTOR_WRAPPER__

#include <PhysBAM_Tools/Krylov_Solvers/KRYLOV_VECTOR_BASE.h>
#include <PhysBAM_Tools/Utilities/PHYSBAM_OVERRIDE.h>
namespace PhysBAM{

template<class T,class TV>
class KRYLOV_VECTOR_WRAPPER: public KRYLOV_VECTOR_BASE<T>
{
public:
    TV v;

    KRYLOV_VECTOR_WRAPPER();
    KRYLOV_VECTOR_WRAPPER(TV vector);
    template<class VECTOR,class INDICES> KRYLOV_VECTOR_WRAPPER(VECTOR& vector,INDICES& index);
    virtual ~KRYLOV_VECTOR_WRAPPER();

    KRYLOV_VECTOR_BASE<T>& operator+=(const KRYLOV_VECTOR_BASE<T>& bv) PHYSBAM_OVERRIDE;
    KRYLOV_VECTOR_BASE<T>& operator-=(const KRYLOV_VECTOR_BASE<T>& bv) PHYSBAM_OVERRIDE;
    KRYLOV_VECTOR_BASE<T>& operator*=(const T a) PHYSBAM_OVERRIDE;
    void Copy(const T c,const KRYLOV_VECTOR_BASE<T>& bv) PHYSBAM_OVERRIDE;
    void Copy(const T c,const KRYLOV_VECTOR_BASE<T>& bv1,const KRYLOV_VECTOR_BASE<T>& bv2) PHYSBAM_OVERRIDE;
    int Raw_Size() const PHYSBAM_OVERRIDE;
    T& Raw_Get(int i) PHYSBAM_OVERRIDE;
//#####################################################################
};
}
#endif
