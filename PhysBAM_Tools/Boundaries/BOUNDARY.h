//#####################################################################
// Copyright 2002-2005, Ronald Fedkiw, Geoffrey Irving, Frank Losasso, Andrew Selle.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class BOUNDARY
//#####################################################################
#ifndef __BOUNDARY__
#define __BOUNDARY__

#include <PhysBAM_Tools/Utilities/NONCOPYABLE.h>
#include <PhysBAM_Tools/Utilities/PHYSBAM_OVERRIDE.h>
#include <PhysBAM_Tools/Vectors/VECTOR_3D.h>
#include <cassert>
namespace PhysBAM{

template<class TV,class T2>
class BOUNDARY:public NONCOPYABLE
{
    typedef typename TV::SCALAR T;typedef VECTOR<bool,2> TV_BOOL2;typedef VECTOR<TV_BOOL2,TV::dimension> TV_SIDES;
public:
    bool use_fixed_boundary,clamp_below,clamp_above;
    T2 fixed_boundary_value,lower_threshold,upper_threshold;
    VECTOR<VECTOR<bool,2>,TV::dimension> constant_extrapolation;

    BOUNDARY()
    {
        Set_Constant_Extrapolation();
        Set_Fixed_Boundary(false);
        Limit_Minimum_Boundary_Value(false);
        Limit_Maximum_Boundary_Value(false);
    }

    virtual ~BOUNDARY()
    {}

    virtual void Set_Constant_Extrapolation(const TV_SIDES& constant_extrapolation_input=TV_SIDES::Constant_Vector(TV_BOOL2::Constant_Vector(true)))
    {constant_extrapolation=constant_extrapolation_input;}

    void Turn_Off_Constant_Extrapolation()
    {Set_Constant_Extrapolation(TV_SIDES());}

    virtual bool Constant_Extrapolation(const int side) const
    {assert(1<=side && side <=6);int axis=(side+1)/2;return constant_extrapolation(axis)(side&1?1:2);}

    virtual void Set_Fixed_Boundary(const bool use_fixed_boundary_input=true,const T2 fixed_boundary_value_input=T2())
    {use_fixed_boundary=use_fixed_boundary_input;fixed_boundary_value=fixed_boundary_value_input;
    if(use_fixed_boundary) clamp_above=clamp_below=false;}

    void Limit_Minimum_Boundary_Value(const bool clamp_below_input=true,const T2 lower_threshold_input=T2())
    {clamp_below=clamp_below_input;lower_threshold=lower_threshold_input;}

    void Limit_Maximum_Boundary_Value(const bool clamp_above_input=true,const T2 upper_threshold_input=T2())
    {clamp_above=clamp_above_input;upper_threshold=upper_threshold_input;}

//#####################################################################
};
}
#endif
