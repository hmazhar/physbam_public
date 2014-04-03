//#####################################################################
// Copyright 2009.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Function MAGNITUDE 
//#####################################################################
#ifndef __MAGNITUDE__
#define __MAGNITUDE__

#include <PhysBAM_Tools/Math_Tools/max.h>
#include <PhysBAM_Tools/Utilities/TYPE_UTILITIES.h>
#include <PhysBAM_Tools/Vectors/VECTOR_BASE.h>
namespace PhysBAM{
template<class T,class T_ARRAY,class ID> class ARRAY_BASE;

namespace ARRAYS_COMPUTATIONS
{
    template<class T,class T_VECTOR>
    inline T Magnitude_Squared(const VECTOR_BASE<T,T_VECTOR>& v)
    {return v.Derived().Magnitude_Squared();}

    inline int Magnitude_Squared(const int a)
    {return a*a;}

    inline float Magnitude_Squared(const float a)
    {return a*a;}

    inline double Magnitude_Squared(const double a)
    {return a*a;}

    template<class T,class T_ARRAY,class ID>
    typename T_ARRAY::SCALAR Magnitude_Squared(const ARRAY_BASE<T,T_ARRAY,ID>& a)
    {const T_ARRAY& self=a.Derived();
    typename T_ARRAY::SCALAR result(0);ID m=self.Size();for(ID i(1);i<=m;i++) result+=Magnitude_Squared(self(i));return result;}

    template<class T,class T_ARRAY,class ID> typename ENABLE_IF<IS_SCALAR<T>::value,T>::TYPE
    Maximum_Magnitude(const ARRAY_BASE<T,T_ARRAY,ID>& a)
    {T result=(T)0;for(int i=1;i<=a.Size();i++) result=PhysBAM::max(result,abs(a(i)));return result;}

    template<class T,class T_ARRAY,class ID> typename T::SCALAR
    Maximum_Magnitude(const ARRAY_BASE<T,T_ARRAY,ID>& a)
    {typename T::SCALAR result(0);for(int i=1;i<=a.Size();i++) result=PhysBAM::max(result,Magnitude_Squared(a(i)));return sqrt(result);}

    template<class T,class T_ARRAY,class ID> typename ENABLE_IF<IS_SCALAR<T>::value,ID>::TYPE
    Arg_Maximum_Magnitude(const ARRAY_BASE<T,T_ARRAY,ID>& a)
    {const T_ARRAY& self=a.Derived();ID m=self.Size();
    T maximum=-1;ID argmax=ID();
    for(ID i(1);i<=m;i++){T current=abs(self(i));if(maximum<current){maximum=current;argmax=i;}}
    return argmax;}

    template<class T,class T_ARRAY,class ID> typename DISABLE_IF<IS_SCALAR<T>::value,ID>::TYPE
    Arg_Maximum_Magnitude(const ARRAY_BASE<T,T_ARRAY,ID>& a)
    {const T_ARRAY& self=a.Derived();ID m=self.Size();
    typename T::SCALAR maximum=-1;ID argmax=ID();
    for(ID i(1);i<=m;i++){
        typename T::SCALAR current=self(i).Magnitude_Squared();
        if(maximum<current){maximum=current;argmax=i;}}
    return argmax;}
}
}
#endif
