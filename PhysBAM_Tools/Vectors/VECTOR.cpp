//#####################################################################
// Copyright 2002-2007, Silvia Salinas-Blemker, Robert Bridson, Ronald Fedkiw, Eran Guendelman, Geoffrey Irving, Neil Molino, Igor Neverov, Craig Schroeder, Joseph Teran.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class VECTOR
//#####################################################################
#include <PhysBAM_Tools/Math_Tools/clamp.h>
#include <PhysBAM_Tools/Math_Tools/max.h>
#include <PhysBAM_Tools/Math_Tools/min.h>
#include <PhysBAM_Tools/Math_Tools/sqr.h>
#include <PhysBAM_Tools/Vectors/VECTOR.h>
#include <cmath>
#include <iostream>
using namespace PhysBAM;
using ::std::abs;
using ::std::floor;
using ::std::ceil;
using ::std::sqrt;
using ::std::exp;
using ::std::sin;
using ::std::cos;
using ::std::pow;
//#####################################################################
// Function Dominant_Axis
//#####################################################################
template<class T,int d> int VECTOR<T,d>::
Dominant_Axis() const
{
    int axis=1;
    T abs_max=abs(array[0]);
    for(int i=1;i<d;i++){T abs_i=abs(array[i]);if(abs_max<abs_i){abs_max=abs_i;axis=i;}}
    return axis;
}
//#####################################################################
// Function abs
//#####################################################################
template<class T,int d> VECTOR<T,d> PhysBAM::
abs(const VECTOR<T,d>& v)
{
    VECTOR<T,d> r;
    for(int i=0;i<d;i++) r.array[i]=abs(v.array[i]);
    return r;
}
//#####################################################################
// Function floor
//#####################################################################
template<class T,int d> VECTOR<T,d> PhysBAM::
floor(const VECTOR<T,d>& v)
{
    VECTOR<T,d> r;
    for(int i=0;i<d;i++) r.array[i]=floor(v.array[i]);
    return r;
}
//#####################################################################
// Function ceil
//#####################################################################
template<class T,int d> VECTOR<T,d> PhysBAM::
ceil(const VECTOR<T,d>& v)
{
    VECTOR<T,d> r;
    for(int i=0;i<d;i++) r.array[i]=ceil(v.array[i]);
    return r;
}
//#####################################################################
// Function rint
//#####################################################################
template<class T,int d> VECTOR<T,d> PhysBAM::
rint(const VECTOR<T,d>& v)
{
    VECTOR<T,d> r;
    for(int i=0;i<d;i++) r.array[i]=rint(v.array[i]);
    return r;
}
//#####################################################################
// Function exp
//#####################################################################
template<class T,int d> VECTOR<T,d> PhysBAM::
exp(const VECTOR<T,d>& v)
{
    VECTOR<T,d> r;
    for(int i=0;i<d;i++) r.array[i]=exp(v.array[i]);
    return r;
}
//#####################################################################
// Function sin
//#####################################################################
template<class T,int d> VECTOR<T,d> PhysBAM::
sin(const VECTOR<T,d>& v)
{
    VECTOR<T,d> r;
    for(int i=0;i<d;i++) r.array[i]=sin(v.array[i]);
    return r;
}
//#####################################################################
// Function cos
//#####################################################################
template<class T,int d> VECTOR<T,d> PhysBAM::
cos(const VECTOR<T,d>& v)
{
    VECTOR<T,d> r;
    for(int i=0;i<d;i++) r.array[i]=cos(v.array[i]);
    return r;
}
//#####################################################################
// Function sqrt
//#####################################################################
template<class T,int d> VECTOR<T,d> PhysBAM::
sqrt(const VECTOR<T,d>& v)
{
    VECTOR<T,d> r;
    for(int i=0;i<d;i++) r.array[i]=sqrt(v.array[i]);
    return r;
}
//#####################################################################
// Function clamp
//#####################################################################
template<class T,int d> VECTOR<T,d> PhysBAM::
clamp(const VECTOR<T,d>& v,const VECTOR<T,d>& vmin,const VECTOR<T,d>& vmax)
{
    VECTOR<T,d> r;
    for(int i=0;i<d;i++) r.array[i]=clamp(v.array[i],vmin.array[i],vmax.array[i]);
    return r;
}
//#####################################################################
// Function clamp
//#####################################################################
template<class T,int d> VECTOR<T,d> PhysBAM::
clamp(const VECTOR<T,d>& v,const T& min,const T& max)
{
    VECTOR<T,d> r;
    for(int i=0;i<d;i++) r.array[i]=clamp(v.array[i],min,max);
    return r;
}
//#####################################################################
// Function clamp_min
//#####################################################################
template<class T,int d> VECTOR<T,d> PhysBAM::
clamp_min(const VECTOR<T,d>& v,const VECTOR<T,d>& vmin)
{
    VECTOR<T,d> r;
    for(int i=0;i<d;i++) r.array[i]=clamp_min(v.array[i],vmin.array[i]);
    return r;
}
//#####################################################################
// Function clamp_min
//#####################################################################
template<class T,int d> VECTOR<T,d> PhysBAM::
clamp_min(const VECTOR<T,d>& v,const T& min)
{
    VECTOR<T,d> r;
    for(int i=0;i<d;i++) r.array[i]=clamp_min(v.array[i],min);
    return r;
}
//#####################################################################
// Function clamp_max
//#####################################################################
template<class T,int d> VECTOR<T,d> PhysBAM::
clamp_max(const VECTOR<T,d>& v,const VECTOR<T,d>& vmax)
{
    VECTOR<T,d> r;
    for(int i=0;i<d;i++) r.array[i]=clamp_max(v.array[i],vmax.array[i]);
    return r;
}
//#####################################################################
// Function clamp_max
//#####################################################################
template<class T,int d> VECTOR<T,d> PhysBAM::
clamp_max(const VECTOR<T,d>& v,const T& max)
{
    VECTOR<T,d> r;
    for(int i=0;i<d;i++) r.array[i]=clamp_max(v.array[i],max);
    return r;
}
//#####################################################################
// Function in_bounds
//#####################################################################
template<class T,int d> VECTOR<T,d> PhysBAM::
in_bounds(const VECTOR<T,d>& v,const VECTOR<T,d>& vmin,const VECTOR<T,d>& vmax)
{
    for(int i=0;i<d;i++) if(!in_bounds(v.array[i],vmin.array[i],vmax.array[i])) return false;
    return true;
}
//#####################################################################
template int VECTOR<float,4>::Dominant_Axis() const;
template int VECTOR<double,4>::Dominant_Axis() const;
template VECTOR<double,4> PhysBAM::clamp<double,4>(VECTOR<double,4> const&,VECTOR<double,4> const&,VECTOR<double,4> const&);
template VECTOR<double,4> PhysBAM::clamp_max<double,4>(VECTOR<double,4> const&,VECTOR<double,4> const&);
template VECTOR<double,4> PhysBAM::clamp_min<double,4>(VECTOR<double,4> const&,VECTOR<double,4> const&);
template VECTOR<double,5> PhysBAM::clamp<double,5>(VECTOR<double,5> const&,VECTOR<double,5> const&,VECTOR<double,5> const&);
template VECTOR<double,5> PhysBAM::clamp_max<double,5>(VECTOR<double,5> const&,VECTOR<double,5> const&);
template VECTOR<double,5> PhysBAM::clamp_min<double,5>(VECTOR<double,5> const&,VECTOR<double,5> const&);
template VECTOR<double,6> PhysBAM::clamp<double,6>(VECTOR<double,6> const&,VECTOR<double,6> const&,VECTOR<double,6> const&);
template VECTOR<double,6> PhysBAM::clamp_max<double,6>(VECTOR<double,6> const&,VECTOR<double,6> const&);
template VECTOR<double,6> PhysBAM::clamp_min<double,6>(VECTOR<double,6> const&,VECTOR<double,6> const&);
template VECTOR<float,4> PhysBAM::clamp<float,4>(VECTOR<float,4> const&,VECTOR<float,4> const&,VECTOR<float,4> const&);
template VECTOR<float,4> PhysBAM::clamp_max<float,4>(VECTOR<float,4> const&,VECTOR<float,4> const&);
template VECTOR<float,4> PhysBAM::clamp_min<float,4>(VECTOR<float,4> const&,VECTOR<float,4> const&);
template VECTOR<float,5> PhysBAM::clamp<float,5>(VECTOR<float,5> const&,VECTOR<float,5> const&,VECTOR<float,5> const&);
template VECTOR<float,5> PhysBAM::clamp_max<float,5>(VECTOR<float,5> const&,VECTOR<float,5> const&);
template VECTOR<float,5> PhysBAM::clamp_min<float,5>(VECTOR<float,5> const&,VECTOR<float,5> const&);
template VECTOR<float,6> PhysBAM::clamp<float,6>(VECTOR<float,6> const&,VECTOR<float,6> const&,VECTOR<float,6> const&);
template VECTOR<float,6> PhysBAM::clamp_max<float,6>(VECTOR<float,6> const&,VECTOR<float,6> const&);
template VECTOR<float,6> PhysBAM::clamp_min<float,6>(VECTOR<float,6> const&,VECTOR<float,6> const&);
template VECTOR<int,4> PhysBAM::clamp<int,4>(VECTOR<int,4> const&,VECTOR<int,4> const&,VECTOR<int,4> const&);
/*template std::ostream& PhysBAM::operator<< <VECTOR<double,2>,4>(std::ostream&,VECTOR<VECTOR<double,2>,4> const&);
template std::ostream& PhysBAM::operator<< <VECTOR<double,3>,4>(std::ostream&,VECTOR<VECTOR<double,3>,4> const&);
template std::ostream& PhysBAM::operator<< <VECTOR<double,4>,4>(std::ostream&,VECTOR<VECTOR<double,4>,4> const&);
template std::ostream& PhysBAM::operator<< <VECTOR<float,2>,4>(std::ostream&,VECTOR<VECTOR<float,2>,4> const&);
template std::ostream& PhysBAM::operator<< <VECTOR<float,3>,4>(std::ostream&,VECTOR<VECTOR<float,3>,4> const&);
template std::ostream& PhysBAM::operator<< <VECTOR<float,4>,4>(std::ostream&,VECTOR<VECTOR<float,4>,4> const&);
template std::ostream& PhysBAM::operator<< <bool,4>(std::ostream&,VECTOR<bool,4> const&);
template std::ostream& PhysBAM::operator<< <bool,6>(std::ostream&,VECTOR<bool,6> const&);
template std::ostream& PhysBAM::operator<< <double,4>(std::ostream&,VECTOR<double,4> const&);
template std::ostream& PhysBAM::operator<< <double,5>(std::ostream&,VECTOR<double,5> const&);
template std::ostream& PhysBAM::operator<< <double,6>(std::ostream&,VECTOR<double,6> const&);
template std::ostream& PhysBAM::operator<< <float,4>(std::ostream&,VECTOR<float,4> const&);
template std::ostream& PhysBAM::operator<< <float,5>(std::ostream&,VECTOR<float,5> const&);
template std::ostream& PhysBAM::operator<< <float,6>(std::ostream&,VECTOR<float,6> const&);
template std::ostream& PhysBAM::operator<< <int,4>(std::ostream&,VECTOR<int,4> const&);
*/
