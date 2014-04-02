//#####################################################################
// Copyright 2002-2007, Robert Bridson, Ronald Fedkiw, Geoffrey Irving, Sergey Koltakov, Nipun Kwatra, Andrew Selle, Jonathan Su, Joseph Teran.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class SEGMENT_3D  
//##################################################################### 
#ifndef __SEGMENT_3D__
#define __SEGMENT_3D__

#include <PhysBAM_Tools/Arrays/ARRAYS_FORWARD.h>
#include <PhysBAM_Tools/Vectors/VECTOR_3D.h>
namespace PhysBAM{

template<class T>
class SEGMENT_3D
{
    typedef VECTOR<T,3> TV;
public:
    VECTOR<T,3> x1,x2;

    SEGMENT_3D()
        :x1(0,0,0),x2(1,0,0)
    {}

    SEGMENT_3D(const VECTOR<T,3>& x1_input,const VECTOR<T,3>& x2_input)
        :x1(x1_input),x2(x2_input)
    {}

    template<class T_ARRAY>
    SEGMENT_3D(const T_ARRAY& X_input)
        :x1(X_input(1)),x2(X_input(2))
    {
        STATIC_ASSERT(T_ARRAY::m==2);
    }

    T Length() const
    {return (x2-x1).Magnitude();}

    T Size() const
    {return Length();}

    static T Size(const TV& x1,const TV& x2)
    {return (x2-x1).Magnitude();}

    template<class T_ARRAY>
    static T Signed_Size(const T_ARRAY& X)
    {STATIC_ASSERT(T_ARRAY::m==2);return Size(X(1),X(2));}

    template<class T_ARRAY>
    static T Size(const T_ARRAY& X)
    {return Signed_Size(X);}

    static VECTOR<T,3> Normal(const VECTOR<T,3>& x1,const VECTOR<T,3>& x2,const VECTOR<T,3>& x3) 
    {VECTOR<T,3> v=x2-x1,face_normal=VECTOR<T,3>::Cross_Product(v,x3-x1);
    return VECTOR<T,3>::Cross_Product(face_normal,v).Normalized();} // rotate by 90 degrees clockwise

    static VECTOR<T,3> Normal_Direction(const VECTOR<T,3>& x1,const VECTOR<T,3>& x2,const VECTOR<T,3>& x3) 
    {VECTOR<T,3> v=x2-x1,face_normal=VECTOR<T,3>::Cross_Product(v,x3-x1);
    return VECTOR<T,3>::Cross_Product(face_normal,v);} // can have any magnitude

    template<class T_ARRAY>
    void Edge_Edge_Interaction_Data(const SEGMENT_3D<T>& segment,const INDIRECT_ARRAY<T_ARRAY,VECTOR<int,4>&> V_edges,const T& distance,VECTOR<T,3>& normal,const VECTOR<T,2>& weights,
        const T small_number=0,const bool verbose=true) const
    {Edge_Edge_Interaction_Data(segment,V_edges(1),V_edges(2),V_edges(3),V_edges(4),distance,normal,weights,small_number,verbose);}

    template<class T_ARRAY>
    bool Edge_Edge_Interaction(const SEGMENT_3D<T>& segment,const INDIRECT_ARRAY<T_ARRAY,VECTOR<int,4>&> V_edges,const T interaction_distance,T& distance,VECTOR<T,3>& normal,
        VECTOR<T,2>& weights,T& relative_speed,bool allow_negative_weights,const T small_number=0,const bool exit_early=false,const bool verbose=true) const
    {return Edge_Edge_Interaction(segment,V_edges(1),V_edges(2),V_edges(3),V_edges(4),interaction_distance,distance,normal,weights,relative_speed,small_number,exit_early,verbose);}

    template<class T_ARRAY>
    bool Edge_Edge_Collision(const SEGMENT_3D<T>& segment,const INDIRECT_ARRAY<T_ARRAY,VECTOR<int,4>&> V_edges,const T dt,const T collision_thickness,T& collision_time,
        VECTOR<T,3>& normal,VECTOR<T,2>& weights,T& relative_speed,const T small_number=0,const bool exit_early=false) const
    {return Edge_Edge_Collision(segment,V_edges(1),V_edges(2),V_edges(3),V_edges(4),dt,collision_thickness,collision_time,normal,weights,relative_speed,small_number,exit_early);}

//#####################################################################
    VECTOR<T,3> Closest_Point_On_Segment(const VECTOR<T,3>& point) const;
    T Distance_From_Point_To_Segment(const VECTOR<T,3>& point) const;
    VECTOR<T,3> Closest_Point_On_Line(const VECTOR<T,3>& point) const;
    T Distance_From_Point_To_Line(const VECTOR<T,3>& point) const;
    VECTOR<T,3> Shortest_Vector_Between_Lines(const SEGMENT_3D<T>& segment,VECTOR<T,2>& weights) const;
    VECTOR<T,3> Shortest_Vector_Between_Segments(const SEGMENT_3D<T>& segment,VECTOR<T,2>& weights) const;
    bool Edge_Edge_Interaction(const SEGMENT_3D<T>& segment,const T interaction_distance,T& distance,VECTOR<T,3>& normal,VECTOR<T,2>& weights,bool allow_negative_weights) const;
    void Edge_Edge_Interaction_Data(const SEGMENT_3D<T>& segment,const TV& v1,const TV& v2,const TV& v3,const TV& v4,
        const T& distance,TV& normal,const VECTOR<T,2>& weights,const T small_number=0,const bool verbose=true) const;
    bool Edge_Edge_Interaction(const SEGMENT_3D<T>& segment,const TV& v1,const TV& v2,const TV& v3,const TV& v4,
        const T interaction_distance,T& distance,VECTOR<T,3>& normal,VECTOR<T,2>& weights,T& relative_speed,bool allow_negative_weights,const T small_number=0,
        const bool exit_early=false,const bool verbose=true) const;
    bool Edge_Edge_Collision(const SEGMENT_3D<T>& segment,const TV& v1,const TV& v2,const TV& v3,const TV& v4,const T dt,
        const T collision_thickness,T& collision_time,VECTOR<T,3>& normal,VECTOR<T,2>& weights,T& relative_speed,const T small_number=0,const bool exit_early=false) const;
    T Interpolation_Fraction(const VECTOR<T,3>& location) const;
    VECTOR<T,2> Barycentric_Coordinates(const VECTOR<T,3>& location) const;
    VECTOR<T,2> Clamped_Barycentric_Coordinates(const VECTOR<T,3>& location,const T tolerance=1e-7) const;
    static T Interpolation_Fraction(const VECTOR<T,3>& location,const VECTOR<T,3>& x1,const VECTOR<T,3>& x2);
    static VECTOR<T,2> Barycentric_Coordinates(const VECTOR<T,3>& location,const VECTOR<T,3>& x1,const VECTOR<T,3>& x2);
    static VECTOR<T,2> Clamped_Barycentric_Coordinates(const VECTOR<T,3>& location,const VECTOR<T,3>& x1,const VECTOR<T,3>& x2,const T tolerance=1e-7);
//#####################################################################
};   
}
#endif
