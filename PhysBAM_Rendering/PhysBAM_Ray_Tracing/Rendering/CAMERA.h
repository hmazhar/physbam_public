//#####################################################################
// Copyright 2002-2007, Robert Bridson, Ronald Fedkiw, Geoffrey Irving, Sergey Koltakov, Andrew Selle.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class CAMERA 
//#####################################################################
//
// Left handed coordinate system: x from left to right, y from bottom to top, z from front to back.
// The camera is positioned at position pointing at the focal point with the film oriented according to the up vector.
// The image plane is located at the focal point.
// The image plane, heigth and width (along with the position, focal_point and vertical_vector) define the viewing frustrum.
// The front and back clipping planes reduce the viewing frustrum to a finite size. 
//
//#####################################################################
#ifndef __CAMERA__
#define __CAMERA__

#include <PhysBAM_Tools/Matrices/MATRIX_4X4.h>
#include <PhysBAM_Tools/Vectors/VECTOR_3D.h>
#include <PhysBAM_Rendering/PhysBAM_Ray_Tracing/Rendering/FILM.h>
namespace PhysBAM{

template<class T>
class CAMERA:public NONCOPYABLE
{
    typedef VECTOR<T,3> TV;

public:
    VECTOR<T,3> position; // camera position 
    VECTOR<T,3> focal_point; // where the image plane is located 
    VECTOR<T,3> look_vector; // points from the position to the focal point - normalized
    VECTOR<T,3> vertical_vector; // point up in the image plane - normalized
    VECTOR<T,3> horizontal_vector; // points to the right on the omage plane - normalized
    FILM<T> film;
    RANDOM_NUMBERS<T> random;

    CAMERA()
        :position(0,0,-1),focal_point(0,0,0),vertical_vector(0,1,0)
    {}

    void Position_And_Aim_Camera(const VECTOR<T,3>& position_input,const VECTOR<T,3>& look_at_point,const VECTOR<T,3>& pseudo_up_vector)
    {position=position_input;look_vector=(look_at_point-position).Normalized();
    horizontal_vector=TV::Cross_Product(look_vector,pseudo_up_vector).Normalized();
    vertical_vector=TV::Cross_Product(horizontal_vector,look_vector).Normalized();}

    void Focus_Camera(const T focal_distance,const T aspect_ratio,const T field_of_view)
    {focal_point=position+focal_distance*look_vector;
    film.width=(T)2*focal_distance*tan((T).5*field_of_view);film.height=film.width/aspect_ratio;}

//#####################################################################
};
}
#endif
