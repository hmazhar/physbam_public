//#####################################################################
// Copyright 2008, Eftychios Sifakis, Jeffrey Hellrung.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#include <PhysBAM_Tools/Adaptive_Arithmetic/EXACT_FLOAT.h>
#include <PhysBAM_Tools/Adaptive_Arithmetic/EXPANSION_ARITHMETIC.h>
using namespace PhysBAM;
//#####################################################################
// operator +
//#####################################################################
template<class T>
EXACT_FLOAT<T> EXACT_FLOAT<T>::
operator+(const EXACT_FLOAT& exact_float) const
{
    EXACT_FLOAT result;
    EXPANSION_ARITHMETIC<T>::Add_Expansions(expansion,exact_float.expansion,result.expansion);
    return result;
}
//#####################################################################
// operator *
//#####################################################################
template<class T>
EXACT_FLOAT<T> EXACT_FLOAT<T>::
operator*(const EXACT_FLOAT& exact_float) const
{
    EXACT_FLOAT result;
    for(int i=1;i<=exact_float.expansion.m;i++) result+=(*this)*exact_float.expansion(i);
    return result;
}
//#####################################################################
// operator *
//#####################################################################
template<class T>
EXACT_FLOAT<T> EXACT_FLOAT<T>::
operator*(T base_float) const
{
    EXACT_FLOAT result;
    EXPANSION_ARITHMETIC<T>::Scale_Expansion(expansion,base_float,result.expansion);
    return result;
}
//#####################################################################
// Function Compress
//#####################################################################
template<class T> inline const EXACT_FLOAT<T>& EXACT_FLOAT<T>::
Compress() const
{
    EXPANSION_ARITHMETIC<T>::Compress_Expansion(const_cast<ARRAY<T>&>(expansion));
    return *this;
}
//#####################################################################
template class EXACT_FLOAT<float>;
template class EXACT_FLOAT<double>;
