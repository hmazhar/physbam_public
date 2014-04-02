//#####################################################################
// Copyright 2003-2006, Ronald Fedkiw, Eran Guendelman, Geoffrey Irving.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class DEEP_WATER_EVOLUTION
//#####################################################################
#include <PhysBAM_Tools/Grids_Uniform/UNIFORM_GRID_ITERATOR_NODE.h>
#include <PhysBAM_Geometry/Fourier_Transforms_Calculations/DEEP_WATER_EVOLUTION_HEIGHTS.h>
using namespace PhysBAM;
//#####################################################################
// Function Intersect_With_Bodies
//#####################################################################
template<class TV> void DEEP_WATER_EVOLUTION_HEIGHTS<TV>::
Intersect_With_Geometry(const RIGID_GEOMETRY<TV_FULL>& rigid_geometry,const T collision_thickness)
{
    bool intersected=false;
    for(NODE_ITERATOR iterator(grid);iterator.Valid();iterator.Next()){TV_INT node=iterator.Node_Index();
        TV pos=!lambda?grid.X(node):Xh(node);
        RAY<TV_FULL> ray(pos.Insert(-10000,2),TV_FULL::Axis_Vector(2));
        ray.semi_infinite=true;

        if(rigid_geometry.Simplex_Intersection(ray,collision_thickness)){
            T object_bottom=ray.Point(ray.t_max).y;
            if(object_bottom<h(node)){h(node)=object_bottom;intersected=true;}}}

    if(intersected) Set_H_Hats_From_Height();
}
//#####################################################################
template class DEEP_WATER_EVOLUTION_HEIGHTS<VECTOR<float,1> >;
template class DEEP_WATER_EVOLUTION_HEIGHTS<VECTOR<float,2> >;
#ifndef COMPILE_WITHOUT_DOUBLE_SUPPORT
template class DEEP_WATER_EVOLUTION_HEIGHTS<VECTOR<double,1> >;
template class DEEP_WATER_EVOLUTION_HEIGHTS<VECTOR<double,2> >;
#endif
