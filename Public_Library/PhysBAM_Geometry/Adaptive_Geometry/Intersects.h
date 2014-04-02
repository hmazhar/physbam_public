//#####################################################################
// Copyright 2008, Jeffrey Hellrung, Andrew Selle, Eftychios Sifakis.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#ifndef __INTERSECTS__
#define __INTERSECTS__
#include <PhysBAM_Tools/Adaptive_Arithmetic/ADAPTIVE_SIGNED_VOLUME.h>
#include <PhysBAM_Geometry/Adaptive_Geometry/Intersection_Coordinates.h>
#include <PhysBAM_Geometry/Adaptive_Geometry/Intersects_And_Intersection_Coordinates_Forward.h>
#include <PhysBAM_Geometry/Adaptive_Geometry/Intersects_Forward.h>

namespace PhysBAM{
namespace GEOMETRIC_PREDICATES_DETAIL{
    //#####################################################################
    // Has_Separating_Hyperplane
    //#####################################################################
    template<class T_EXACT,class T,int m,int n> bool
    Has_Separating_Hyperplane(const VECTOR<VECTOR<T,1>,m>& simplex1,const VECTOR<VECTOR<T,1>,n>& simplex2,bool* is_degenerate_p);
    template<class T_EXACT,class T,int m,int n> bool
    Has_Separating_Hyperplane(const VECTOR<VECTOR<T,2>,m>& simplex1,const VECTOR<VECTOR<T,2>,n>& simplex2,bool* is_degenerate_p);
    template<class T_EXACT,class T,int m,int n> bool
    Has_Separating_Hyperplane(const VECTOR<VECTOR<T,3>,m>& simplex1,const VECTOR<VECTOR<T,3>,n>& simplex2,bool* is_degenerate_p);
    template<class T_EXACT,class T,int d,int m,int n> bool
    Is_Cohyperplanar(const VECTOR<VECTOR<T,d>,m>& simplex1,const VECTOR<VECTOR<T,d>,n>& simplex2);

    //#####################################################################
    // Intersects_Or_Is_Degenerate_After_Projection
    //#####################################################################
    template<class T_EXACT,class T> bool
    Intersects_Or_Is_Degenerate_After_Projection(const VECTOR<VECTOR<T,1>,1>& simplex1,const VECTOR<VECTOR<T,1>,1>& simplex2)
    {assert(Is_Cohyperplanar<T_EXACT>(simplex1,simplex2));return true;}
    
    template<class T_EXACT,class T,int d,int m,int n> bool
    Intersects_Or_Is_Degenerate_After_Projection(const VECTOR<VECTOR<T,d>,m>& simplex1,const VECTOR<VECTOR<T,d>,n>& simplex2);
    template<class T_EXACT,class T,int d,int n> bool
    Is_Any_Inside_Helper(const VECTOR<VECTOR<T,d>,d+1>& simplex1,const VECTOR<VECTOR<T,d>,n>& simplex2);
    
    template<class T_EXACT,class T,int d,int n> typename ENABLE_IF<(n<d+1),bool>::TYPE
    Is_Any_Inside(const VECTOR<VECTOR<T,d>,d+1>& simplex1,const VECTOR<VECTOR<T,d>,n>& simplex2)
    {return Is_Any_Inside_Helper<T_EXACT>(simplex1,simplex2);}
    
    template<class T_EXACT,class T,int d> bool
    Is_Any_Inside(const VECTOR<VECTOR<T,d>,d+1>& simplex1,const VECTOR<VECTOR<T,d>,d+1>& simplex2)
    {return Is_Any_Inside_Helper<T_EXACT>(simplex1,simplex2)||Is_Any_Inside_Helper<T_EXACT>(simplex2,simplex1);}

//#################################################################
// Intersects
//#################################################################
template<class T_EXACT,class T,int d> bool
Intersects(const VECTOR<T,d>& point,const VECTOR<VECTOR<T,d>,d+1>& simplex,bool* is_degenerate_p)
{
    return Intersects<T_EXACT>(simplex,point,is_degenerate_p);
}

template<class T_EXACT,class T,int d> bool
Intersects(const VECTOR<VECTOR<T,d>,d+1>& simplex,const VECTOR<T,d>& point,bool* is_degenerate_p);

template<class T_EXACT,class T,int d,int m,int n> typename ENABLE_IF<(m<n),bool>::TYPE
Intersects(const VECTOR<VECTOR<T,d>,m>& simplex1,const VECTOR<VECTOR<T,d>,n>& simplex2,bool* is_degenerate_p)
{
    return Intersects<T_EXACT>(simplex2,simplex1,is_degenerate_p);
}

template<class T_EXACT,class T,int d> bool
Intersects(const VECTOR<VECTOR<T,d>,d+1>& simplex1,const VECTOR<VECTOR<T,d>,1>& point,bool* is_degenerate_p)
{
    return Intersects<T_EXACT>(simplex1,point[1],is_degenerate_p);
}

template<class T_EXACT,class T,int d,int m,int n> typename ENABLE_IF<(n<=m&&m<=d&&m+n<=d+2),bool>::TYPE
Intersects(const VECTOR<VECTOR<T,d>,m>& simplex1,const VECTOR<VECTOR<T,d>,n>& simplex2,bool* is_degenerate_p);
template<class T_EXACT,class T,int d,int m,int n> typename ENABLE_IF<(n<=m&&m<=d&&m+n>=d+3),bool>::TYPE
Intersects(const VECTOR<VECTOR<T,d>,m>& simplex1,const VECTOR<VECTOR<T,d>,n>& simplex2,bool* is_degenerate_p);
template<class T_EXACT,class T,int d,int n> bool
Intersects(const VECTOR<VECTOR<T,d>,d+1>& simplex1,const VECTOR<VECTOR<T,d>,n>& simplex2,bool* is_degenerate_p);
template<class T_EXACT,class T>
bool Intersects(const VECTOR<VECTOR<T,3>,3>& triangle1,const VECTOR<VECTOR<T,3>,3>& triangle2,const VECTOR<VECTOR<T,3>,3>& triangle3,bool* is_degenerate_p);
//#################################################################
}
using GEOMETRIC_PREDICATES_DETAIL::Intersects;
}
#endif
