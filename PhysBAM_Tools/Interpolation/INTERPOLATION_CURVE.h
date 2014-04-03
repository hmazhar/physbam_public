//#####################################################################
// Copyright 2004-2008, Ron Fedkiw, Eran Guendelman, Geoffrey Irving, Craig Schroeder.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class INTERPOLATION_CURVE
//#####################################################################
//
// A control point's value controls the interval to the right of the control point.
//
//#####################################################################
#ifndef __INTERPOLATION_CURVE__
#define __INTERPOLATION_CURVE__

#include <PhysBAM_Tools/Arrays/ARRAY.h>
#include <PhysBAM_Tools/Data_Structures/TRIPLE.h>
#include <PhysBAM_Tools/Matrices/MATRIX_POLICY.h>
#include <PhysBAM_Tools/Vectors/VECTOR_FORWARD.h>
namespace PhysBAM{

template<class T,class S> struct INTERPOLATION_CURVE_HELPER
{
    typedef S DERIVATIVE;

    static S Value(const T t_left,const T t_right,const S& u_left,const S& u_right,const T t)
    {return u_left+(t-t_left)*(u_right-u_left)/(t_right-t_left);}

    static S Difference(const S& u_left,const S& u_right)
    {return u_right-u_left;}
};

template<class T> struct INTERPOLATION_CURVE_HELPER<T,bool> // enforce constant interpolation for bools
{
    struct UNUSABLE{};
    typedef UNUSABLE DERIVATIVE;

    static bool Value(const T t_left,const T t_right,const bool u_left,const bool u_right,const T t)
    {return u_left;}
};

template<class T,class TV> struct INTERPOLATION_CURVE_HELPER<T,ROTATION<TV> >
{
    typedef typename TV::SPIN DERIVATIVE;

    static ROTATION<TV> Value(const T t_left,const T t_right,const ROTATION<TV>& u_left,const ROTATION<TV>& u_right,const T t)
    {return ROTATION<TV>::Spherical_Linear_Interpolation(u_left,u_right,(t-t_left)/(t_right-t_left));}

    static DERIVATIVE Difference(const ROTATION<TV>& u_left,const ROTATION<TV>& u_right)
    {return (u_right*u_left.Inverse()).Rotation_Vector();}
};

template<class T,class TV> struct INTERPOLATION_CURVE_HELPER<T,FRAME<TV> >
{
    typedef TWIST<TV> DERIVATIVE;

    static FRAME<TV> Value(const T t_left,const T t_right,const FRAME<TV>& u_left,const FRAME<TV>& u_right,const T t)
    {return FRAME<TV>::Interpolation(u_left,u_right,(t-t_left)/(t_right-t_left));}

    static TWIST<TV> Difference(const FRAME<TV>& u_left,const FRAME<TV>& u_right)
    {return TWIST<TV>(u_right.t-u_left.t,(u_right.r*u_left.r.Inverse()).Rotation_Vector());}
};

template<class T,class T2>
class INTERPOLATION_CURVE
{
public:
    enum INTERPOLATION_TYPE {CONSTANT,LINEAR};
    struct CONTROL_POINT
    {
        T t;T2 value;INTERPOLATION_TYPE type;

        CONTROL_POINT()
            :t(),value(),type(LINEAR)
        {}

        CONTROL_POINT(const T t,const T2& value,const INTERPOLATION_TYPE type=LINEAR)
            :t(t),value(value),type(type)
        {}
    };

    ARRAY<CONTROL_POINT> control_points;

    INTERPOLATION_CURVE()
    {}

    template<class S> struct HELPER:public INTERPOLATION_CURVE_HELPER<T,S>{};
    typedef typename HELPER<T2>::DERIVATIVE DERIVATIVE;

    int Locate_Interval(const T t) const
    {int lo=0,hi=control_points.m+1;
    while(hi-lo>1){int m=(lo+hi)/2;if(t<control_points(m).t) hi=m;else lo=m;}
    return lo;}

    T2 Value(T t) const
    {int left_index=Locate_Interval(t);
    if(left_index==0) return control_points.m>0?control_points(1).value:T2();
    else if(left_index==control_points.m) return control_points(left_index).value;
    else if(control_points(left_index).type==CONSTANT) return control_points(left_index).value;
    else return HELPER<T2>::Value(control_points(left_index).t,control_points(left_index+1).t,control_points(left_index).value,control_points(left_index+1).value,t);}

    DERIVATIVE Derivative(T t) const
    {int left_index=Locate_Interval(t);
    if(left_index==0 || left_index==control_points.m || control_points(left_index).type==CONSTANT) return DERIVATIVE();
    return 1/(control_points(left_index+1).t-control_points(left_index).t)*HELPER<T2>::Difference(control_points(left_index).value,control_points(left_index+1).value);}

    void Add_Control_Point(T t,const T2& value,INTERPOLATION_TYPE type=LINEAR)
    {int index=Locate_Interval(t);
    control_points.Insert(CONTROL_POINT(t,value,type),index+1);}

//#####################################################################
};
}
#endif
