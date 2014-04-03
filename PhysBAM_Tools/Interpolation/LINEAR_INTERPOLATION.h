//#####################################################################
// Copyright 2003-2005, Ronald Fedkiw, Eran Guendelman, Geoffrey Irving, Frank Losasso, Andrew Selle.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class LINEAR_INTERPOLATION 
//#####################################################################
#ifndef __LINEAR_INTERPOLATION__
#define __LINEAR_INTERPOLATION__

#include <PhysBAM_Tools/Vectors/VECTOR.h>
namespace PhysBAM{

template<class T,class T2>
class LINEAR_INTERPOLATION
{
public:
    static T2 Linear(const T2& u_left,const T2& u_right,const T x)
    {return (1-x)*u_left+x*u_right;}

    static T2 Linear(const T2& u_left,const T2& u_right,const VECTOR<T,1>& X)
    {return Linear(u_left,u_right,X.x);}

    static T2 Linear(const T x_left,const T x_right,const T2& u_left,const T2& u_right,const T x)
    {return u_left+(x-x_left)*(u_right-u_left)/(x_right-x_left);}

    static T2 Linear(const T x_left,const T2& u_left,const T2& u_slope,const T x)
    {return u_left+(x-x_left)*u_slope;}

    static T2 Linear(const T x_left,const T2& u_left,const T2& u_slope,const VECTOR<T,1> X)
    {return u_left+(X.x-x_left)*u_slope;}

    static T2 Linear_Predivided(const T x_left,const T one_over_x_right_minus_x_left,const T2& u_left,const T2& u_right,const T x)
    {return u_left+(x-x_left)*one_over_x_right_minus_x_left*(u_right-u_left);}

    static T2 Linear_Normalized(const T2& u_left,const T2& u_slope,const T x)
    {return u_left+x*u_slope;}

    static T2 Linear(const T2 nodes[2],const VECTOR<T,1>& X)
    {return Linear(nodes[0],nodes[1],X.x);}

    static T2 Linear(const T2 nodes[4],const VECTOR<T,2>& X)
    {return Bilinear(nodes[0],nodes[1],nodes[2],nodes[3],X);}

    static T2 Linear(const T2 nodes[8],const VECTOR<T,3>& X)
    {return Trilinear(nodes[0],nodes[1],nodes[2],nodes[3],nodes[4],nodes[5],nodes[6],nodes[7],X);}

    static T2 Bilinear(const T2& u1,const T2& u2,const T2& u3,const T2& u4,const VECTOR<T,2>& minimum_corner,const VECTOR<T,2>& maximum_corner,const VECTOR<T,2>& X);
    // X in [0,1]x[0,1]
    static T2 Bilinear(const T2& u1,const T2& u2,const T2& u3,const T2& u4,const VECTOR<T,2>& X);
    static T2 Bilinear(const T2& u1,const T2& u3,T one_over_y_top_minus_y_bottom,const T x_left,const T y_bottom,const T2& slope12,const T2& slope34,const VECTOR<T,2>& X);
    static T2 Trilinear(const T2& u1,const T2& u2,const T2& u3,const T2& u4,const T2& u5,const T2& u6,const T2& u7,const T2& u8,
        const VECTOR<T,3>& minimum_corner,const VECTOR<T,3>& maximum_corner,const VECTOR<T,3>& X);
    // X in [0,1]x[0,1]x[0,1]
    static T2 Trilinear(const T2& u1,const T2& u2,const T2& u3,const T2& u4,const T2& u5,const T2& u6,const T2& u7,const T2& u8,const VECTOR<T,3>& X);
    static T2 Trilinear(const T2& u1,const T2& u3,const T2& u5,const T2& u7,T one_over_y_top_minus_y_bottom,T one_over_z_back_minus_z_front,const T x_left,const T y_bottom,const T z_front,
        const T2& slope12,const T2& slope34,const T2& slope56,const T2& slope78,const VECTOR<T,3>& X);

//#####################################################################
};
}
#endif
