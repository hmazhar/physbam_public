//#####################################################################
// Copyright 2009, Craig Schroeder.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#include <PhysBAM_Tools/Krylov_Solvers/KRYLOV_VECTOR_BASE.h>
#include <PhysBAM_Tools/Log/DEBUG_UTILITIES.h>
using namespace PhysBAM;
//#####################################################################
// Constructor
//#####################################################################
template<class T> KRYLOV_VECTOR_BASE<T>::
KRYLOV_VECTOR_BASE()
{
}
//#####################################################################
// Destructor
//#####################################################################
template<class T> KRYLOV_VECTOR_BASE<T>::
~KRYLOV_VECTOR_BASE()
{
}
//#####################################################################
// Function Raw_Get
//#####################################################################
template<class T> const T& KRYLOV_VECTOR_BASE<T>::
Raw_Get(int i) const
{
    return const_cast<KRYLOV_VECTOR_BASE<T>*>(this)->Raw_Get(i);
}
template class KRYLOV_VECTOR_BASE<float>;
#ifndef COMPILE_WITHOUT_DOUBLE_SUPPORT
template class KRYLOV_VECTOR_BASE<double>;
#endif
