//#####################################################################
// Copyright 2005-2007, Geoffrey Irving, Nipun Kwatra, Andrew Selle.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class IMAGE
//#####################################################################
#ifndef COMPILE_WITHOUT_READ_WRITE_SUPPORT
#ifndef __IMAGE__
#define __IMAGE__

#include <PhysBAM_Tools/Grids_Uniform_Arrays/ARRAYS_ND.h>
#include <PhysBAM_Tools/Math_Tools/exchange.h>
#include <PhysBAM_Tools/Vectors/VECTOR_3D.h>
namespace PhysBAM{

template<class T>
class IMAGE
{
public:
    IMAGE()
    {}

    template<int d> static VECTOR<T,d> Byte_Color_To_Scalar_Color(const VECTOR<unsigned char,d> color_in)
    {return (VECTOR<T,d>(color_in)+VECTOR<T,d>::All_Ones_Vector())/(T)512;}

    static T Byte_Color_To_Scalar_Color(const unsigned char color_in)
    {return ((T)color_in+(T).5)/256;}
    
    template<int d> static VECTOR<unsigned char,d> Scalar_Color_To_Byte_Color(const VECTOR<T,d> color_in)
    {return VECTOR<unsigned char,d>(clamp(VECTOR<int,d>(VECTOR<T,d>((T)256*color_in)),VECTOR<int,d>(),VECTOR<int,d>(255*VECTOR<int,d>::All_Ones_Vector())));}

    static unsigned char Scalar_Color_To_Byte_Color(const T color_in)
    {return clamp((int)((T)256*color_in),0,255);}

    template<int d> static void Flip_X(ARRAY<VECTOR<T,d> ,VECTOR<int,2> >& image)
    {for(int i=1;i<=image.m/2;i++)for(int j=1;j<=image.n;j++)exchange(image(i,j),image(image.m-i+1,j));}

    template<int d> static void Flip_Y(ARRAY<VECTOR<T,d> ,VECTOR<int,2> >& image)
    {for(int j=1;j<=image.n/2;j++)for(int i=1;i<=image.m;i++)exchange(image(i,j),image(i,image.n-j+1));}

    template<int d> static void Invert(ARRAY<VECTOR<T,d> ,VECTOR<int,2> >& image)
    {for(int i=image.domain.min_corner.x;i<=image.domain.max_corner.x;i++) for(int j=image.domain.min_corner.y;j<=image.domain.max_corner.y;j++) image(i,j)=VECTOR<T,d>::All_Ones_Vector()-image(i,j);}

    static void Threshold(ARRAY<VECTOR<T,3> ,VECTOR<int,2> >& image,const T threshold,const VECTOR<T,3>& low_color,const VECTOR<T,3>& high_color)
    {for(int i=image.domain.min_corner.x;i<=image.domain.max_corner.x;i++) for(int j=image.domain.min_corner.y;j<=image.domain.max_corner.y;j++) if(image(i,j).Magnitude()<threshold) image(i,j)=low_color;else image(i,j)=high_color;}

//#####################################################################
    template<int d> static void Read(const std::string& filename,ARRAY<VECTOR<T,d> ,VECTOR<int,2> >& image);
    template<int d> static void Write(const std::string& filename,const ARRAY<VECTOR<T,d> ,VECTOR<int,2> >& image,const T gamma=1,const T dither_amplitude=0);
    static bool Is_Supported(const std::string& filename);
//#####################################################################
};
}
#endif
#endif
