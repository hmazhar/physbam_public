//#####################################################################
// Copyright 2009, Jon Gretarsson.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class LEVELSET_UTILITIES
//#####################################################################
#ifndef __LEVELSET_UTILITIES__
#define __LEVELSET_UTILITIES__

#include <PhysBAM_Tools/Math_Tools/constants.h>
#include <PhysBAM_Tools/Utilities/NONCOPYABLE.h>
#include <PhysBAM_Tools/Utilities/PHYSBAM_OVERRIDE.h>
#include <PhysBAM_Tools/Vectors/SCALAR_POLICY.h>
#include <cassert>
#include <cfloat>
namespace PhysBAM{

template<class T>
class LEVELSET_UTILITIES:public NONCOPYABLE
{
    STATIC_ASSERT((IS_SCALAR<T>::value));
protected:
    LEVELSET_UTILITIES() {};
    virtual ~LEVELSET_UTILITIES() {};

public:
    static T Sign(const T phi)
    {if(phi <= 0) return -1;else return 1;}

    static int Interface(const T phi_1,const T phi_2)
    {if((phi_1 > 0 && phi_2 <=0) || (phi_1 <= 0 && phi_2 > 0)) return 1;else return 0;}

    static int Thin_Shells_Interface(const T phi_1,const T phi_2)
    {if(phi_1==FLT_MAX||phi_2==FLT_MAX) return 0;else return Interface(phi_1,phi_2);}

    static T Heaviside(const T phi,const T half_width=0)
    {if(phi <= -half_width) return 0;if(phi >= half_width) return 1;
    T phi_over_half_width=phi/half_width;return (T).5*(1+phi_over_half_width+sin((T)pi*phi_over_half_width)/(T)pi);}

    static T Heaviside(const T phi,const T value_minus,const T value_plus,const T half_width=0)
    {return value_minus+(value_plus-value_minus)*Heaviside(phi,half_width);}

    static T Delta(const T phi,const T half_width)
    {if(phi <= -half_width || phi >= half_width) return 0;else return (T)(1+cos(pi*phi/half_width))/(2*half_width);}

    // finds the interface as an average of the neighbors
    static T Average(const T phi_left,const T value_left,const T phi_right,const T value_right)
    {T left=std::abs(phi_left),right=std::abs(phi_right),sum=left+right;
    if(right<sum) return value_right+right/sum*(value_left-value_right);
    return value_left;}

    // finds the interface value, while sorting out the sign of phi
    static T Convex_Average(const T phi_1,const T phi_2,const T value_minus,const T value_plus)
    {if(phi_1 <= 0) return Average(phi_1,value_minus,phi_2,value_plus);else return Average(phi_2,value_minus,phi_1,value_plus);}

    static T Theta(const T phi_left,const T phi_right)
    {assert(phi_left!=phi_right);return phi_left/(phi_left-phi_right);}

//#####################################################################
    static T Theta_Quadratic(const T phi_left_left,const T phi_left,const T phi_right,const T dx);
    static T Theta_Cubic(const T phi_left_left,const T phi_left,const T phi_right,const T phi_right_right,const T dx);
    static T Negative_Cell_Fraction(const T phi_lower_left,const T phi_lower_right,const T phi_upper_left,const T phi_upper_right,const T aspect_ratio=1);
//#####################################################################
};
}
#endif
