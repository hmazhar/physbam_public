//#####################################################################
// Copyright 2007-2009, Jon Gretarsson, Geoffrey Irving, Nipun Kwatra, Andrew Selle, Jonathan Su.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class POINT_SIMPLEX_1D  
//##################################################################### 
#ifndef __POINT_SIMPLEX_1D__
#define __POINT_SIMPLEX_1D__

#include <PhysBAM_Tools/Arrays/ARRAYS_FORWARD.h>
#include <PhysBAM_Tools/Log/DEBUG_UTILITIES.h>
#include <PhysBAM_Tools/Math_Tools/ONE.h>
#include <PhysBAM_Tools/Utilities/TYPE_UTILITIES.h>
#include <PhysBAM_Tools/Vectors/VECTOR_1D.h>
namespace PhysBAM{

template<class TV> class RANGE;

template<class T>
class POINT_SIMPLEX_1D
{
    typedef VECTOR<T,1> TV;
public:
    typedef TV VECTOR_T;

    TV x1;
    bool direction; // false for -1 and true for 1 direction
 
    POINT_SIMPLEX_1D()
        :x1(0),direction(true)
    {}

    POINT_SIMPLEX_1D(const TV& x1,const bool direction)
        :x1(x1),direction(direction)
    {}

    template<class T_ARRAY>
    explicit POINT_SIMPLEX_1D(const T_ARRAY& X_input)
        :x1(X_input(1)),direction(true)
    {
        STATIC_ASSERT(T_ARRAY::m==1); // TODO: initialize directions properly
    }

    const TV& Center() const
    {return x1;}

    T Size() const
    {return (T)1;}

    static TV Normal(const TV& x1,const int direction)
    {return TV((T)(direction?1:-1));}

    TV Normal() const
    {return POINT_SIMPLEX_1D<T>::Normal(x1,direction);}

    TV Normal(const TV& location) const
    {return Normal();}

    static ONE Clamped_Barycentric_Coordinates(const TV& location,const TV& x1)
    {return ONE();}

    template<class T_ARRAY>
    static typename ENABLE_IF<AND<IS_SAME<typename T_ARRAY::ELEMENT,TV>::value,T_ARRAY::m==1>::value,ONE>::TYPE
    Clamped_Barycentric_Coordinates(const TV& location,const T_ARRAY& X)
    {return ONE();}

    TV Sum_Barycentric_Coordinates(const POINT_SIMPLEX_1D<T>& embedded_point_simplex) const
    {return TV::All_Ones_Vector();} // TODO

    ONE Barycentric_Coordinates(const TV& location) const
    {return ONE();}

    static TV Point_From_Barycentric_Coordinates(const ONE,const TV& x1)
    {return x1;}

    template<class T_ARRAY>
    static typename ENABLE_IF<AND<IS_SAME<typename T_ARRAY::ELEMENT,TV>::value,T_ARRAY::m==1>::value,TV>::TYPE
    Point_From_Barycentric_Coordinates(const ONE,const T_ARRAY& X)
    {return X(1);}

    TV Point_From_Barycentric_Coordinates(const ONE) const
    {return x1;}

    RANGE<TV> Bounding_Box() const
    {return RANGE<TV>(x1);}

    T Signed_Distance(const TV& location) const
    {return (direction?(T)1:(T)-1)*(location.x-x1.x);}

    const TV& X(const int i) const
    {assert(1==i);return x1;}

    TV& X(const int i)
    {assert(1==i);return x1;}

    void Clip_To_Box(const RANGE<TV>& box,ARRAY<POINT_SIMPLEX_1D<T> >& clipped_simplices) const
    {clipped_simplices.Remove_All();
    if(box.Lazy_Inside(x1)) clipped_simplices.Append(*this);}

    bool Inside(const TV& point,const T thickness_over_two=0) const
    {if(x1.x-thickness_over_two <= point.x && x1.x+thickness_over_two >= point.x) return true;
    else return false;}

    VECTOR<T,0> Principal_Curvatures(const TV& X) const PHYSBAM_OVERRIDE
    {return VECTOR<T,0>();}

    static std::string Name()
    {return "POINT_SIMPLEX_1D<T>";}
};
}
#endif
