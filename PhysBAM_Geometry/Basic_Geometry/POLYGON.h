//#####################################################################
// Copyright 2002-2005, Robert Bridson, Ronald Fedkiw, Geoffrey Irving.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class POLYGON  
//##################################################################### 
#ifndef __POLYGON__
#define __POLYGON__

#include <PhysBAM_Tools/Arrays/ARRAY.h>
#include <PhysBAM_Tools/Math_Tools/RANGE.h>
#include <PhysBAM_Tools/Vectors/VECTOR_3D.h>
#include <PhysBAM_Geometry/Basic_Geometry/SEGMENT_2D.h>
#include <cfloat>
namespace PhysBAM{

template<class TV>
class POLYGON 
{
    typedef typename TV::SCALAR T;
public:
    ARRAY<TV> X; // vertex coordinates

    POLYGON();
    POLYGON(const int number_of_vertices);
    POLYGON(const RANGE<TV>& range);

//#####################################################################
    T Area() const;
    TV Find_Closest_Point_On_Polygon(const TV& X_point,int& side) const;
    T Distance_From_Polygon_To_Point(const TV& X_point) const;
    bool Inside_Polygon(const TV& X_point) const;
    T Signed_Distance(const TV& location) const;
//#####################################################################
};   
//#####################################################################
}
#endif
