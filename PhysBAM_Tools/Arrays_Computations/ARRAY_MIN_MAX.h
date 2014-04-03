//#####################################################################
// Copyright 2009, Craig Schroeder, Jonathan Su.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class ARRAY_MIN_MAX
//#####################################################################
#ifndef __ARRAY_MIN_MAX__
#define __ARRAY_MIN_MAX__

#include <PhysBAM_Tools/Math_Tools/max.h>
#include <PhysBAM_Tools/Math_Tools/maxabs.h>
#include <PhysBAM_Tools/Math_Tools/maxmag.h>
#include <PhysBAM_Tools/Math_Tools/min.h>
#include <PhysBAM_Tools/Math_Tools/minmag.h>
namespace PhysBAM{
template<class T,class T_ARRAY,class ID>  class ARRAY_BASE;

namespace ARRAYS_COMPUTATIONS
{
    template<class T,class T_ARRAY,class ID>
    T Max(const ARRAY_BASE<T,T_ARRAY,ID>& a)
    {const T_ARRAY& self=a.Derived();T result=self(ID(1));ID m=self.Size();for(ID i(2);i<=m;i++) result=PhysBAM::max(result,self(i));return result;}

    template<class T,class T_ARRAY,class ID>
    T Maxabs(const ARRAY_BASE<T,T_ARRAY,ID>& a)
    {const T_ARRAY& self=a.Derived();T result=T();;ID m=self.Size();for(ID i(1);i<=m;i++) result=PhysBAM::max(result,abs(self(i)));return result;}

    template<class T,class T_ARRAY,class ID>
    T Maxmag(const ARRAY_BASE<T,T_ARRAY,ID>& a)
    {const T_ARRAY& self=a.Derived();T result=T();ID m=self.Size();for(ID i(1);i<=m;i++) result=PhysBAM::maxmag(result,self(i));return result;}

    template<class T,class T_ARRAY,class ID>
    ID Argmax(const ARRAY_BASE<T,T_ARRAY,ID>& a)
    {const T_ARRAY& self=a.Derived();ID result(1),m=self.Size();for(ID i(2);i<=m;i++) if(self(i)>self(result)) result=i;return result;}

    template<class T,class T_ARRAY,class ID>
    T Min(const ARRAY_BASE<T,T_ARRAY,ID>& a)
    {const T_ARRAY& self=a.Derived();T result=self(ID(1));ID m=self.Size();for(ID i(2);i<=m;i++) result=PhysBAM::min(result,self(i));return result;}

    template<class T,class T_ARRAY,class ID>
    T Minmag(const ARRAY_BASE<T,T_ARRAY,ID>& a)
    {const T_ARRAY& self=a.Derived();T result=self(ID(1));ID m=self.Size();for(ID i(2);i<=m;i++) result=PhysBAM::minmag(result,self(i));return result;}

    template<class T,class T_ARRAY,class ID>
    ID Argmin(const ARRAY_BASE<T,T_ARRAY,ID>& a)
    {const T_ARRAY& self=a.Derived();ID result(1),m=self.Size();for(ID i(2);i<=m;i++) if(self(i)<self(result)) result=i;return result;}

    template<class T,class T_ARRAY,class ID>
    T Componentwise_Maxabs(const ARRAY_BASE<T,T_ARRAY,ID>& a)
    {const T_ARRAY& self=a.Derived();T result=T();ID m=self.Size();for(ID i(1);i<=m;i++) result=T::Componentwise_Max(result,abs(self(i)));return result;}
}
}
#endif
