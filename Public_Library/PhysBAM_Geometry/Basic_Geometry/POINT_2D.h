//#####################################################################
// Copyright 2007, Nipun Kwatra, Jonathan Su
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class POINT_2D  
//##################################################################### 
#ifndef __POINT_2D__
#define __POINT_2D__

#include <PhysBAM_Tools/Arrays/INDIRECT_ARRAY.h>
#include <PhysBAM_Tools/Log/DEBUG_UTILITIES.h>
#include <PhysBAM_Tools/Log/LOG.h>
#include <PhysBAM_Tools/Polynomials/QUADRATIC.h>
#include <PhysBAM_Tools/Utilities/TYPE_UTILITIES.h>
#include <PhysBAM_Tools/Vectors/VECTOR_2D.h>
#include <PhysBAM_Geometry/Continuous_Collision_Detection/EDGE_EDGE_COLLISION.h>
namespace PhysBAM{

template<class T>
class POINT_2D:public VECTOR<T,2>
{
    typedef VECTOR<T,2> TV;
public:
    POINT_2D()
    {}

    POINT_2D(const TV& x_input)
        :TV(x_input)
    {}

    bool Edge_Edge_Interaction(const POINT_2D<T>& point,const T interaction_distance,T& distance,VECTOR<T,2>& normal) const
    {normal=*this-point;distance=normal.Magnitude();return distance<=interaction_distance;}

    void Edge_Edge_Interaction_Data(const POINT_2D<T>& point,const TV& v1,const TV& v2,const T& distance,TV& normal,const T small_number) const
    {if(distance > small_number) normal/=distance;
    else normal=(v1-v2).Normalized();}

    template<class T_ARRAY>
    bool Edge_Edge_Interaction(const POINT_2D<T>& point,const INDIRECT_ARRAY<T_ARRAY,VECTOR<int,2>&> V_edges,const T interaction_distance,
        T& distance,VECTOR<T,2>& normal,VECTOR<T,2>& weights,T& relative_speed,bool allow_negaqtive_weights,const T small_number,const bool exit_early=false) const
    {if(!Edge_Edge_Interaction(point,interaction_distance,distance,normal)) return false;
    if(!exit_early){
        Edge_Edge_Interaction_Data(point,V_edges(1),V_edges(2),distance,normal,small_number);
        relative_speed=VECTOR<T,2>::Dot_Product(V_edges(1)-V_edges(2),normal);} // relative speed is in the normal direction
    return true;}

    template<class T_ARRAY>
    bool Edge_Edge_Collision(const POINT_2D<T>& point,const INDIRECT_ARRAY<T_ARRAY,VECTOR<int,2>&> V_edges,const T dt,const T collision_thickness,T& collision_time,
        VECTOR<T,2>& normal,VECTOR<T,2>& weights,T& relative_speed,bool allow_negative_weights,const T small_number=0,const bool exit_early=false) const
    {return CONTINUOUS_COLLISION_DETECTION_COMPUTATIONS::Edge_Edge_Collision(*this,point,V_edges,dt,collision_thickness,collision_time,normal,weights,relative_speed,
            allow_negative_weights,small_number,exit_early);}
};
}
#endif
