//#####################################################################
// Copyright 2004-2005, Andrew Selle.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class UNIFORM_BOX_PARTITION
//#####################################################################
#ifndef __UNIFORM_BOX_PARTITION__
#define __UNIFORM_BOX_PARTITION__

#include <PhysBAM_Tools/Arrays/ARRAY.h>
#include <PhysBAM_Tools/Data_Structures/PAIR.h>
#include <PhysBAM_Tools/Grids_Uniform/GRID.h>
#include <PhysBAM_Tools/Utilities/NONCOPYABLE.h>
#include <PhysBAM_Geometry/Basic_Geometry/RAY.h>
#include <PhysBAM_Geometry/Basic_Geometry_Intersections/RAY_BOX_INTERSECTION.h>
namespace PhysBAM{

template<class T,class DATA_T>
class UNIFORM_BOX_PARTITION:public NONCOPYABLE
{
    typedef VECTOR<T,3> TV;typedef VECTOR<int,3> TV_INT;
public:
    GRID<TV> grid;
    ARRAY<ARRAY<DATA_T>*,VECTOR<int,3> > cells;
    RANGE<TV> bounding_box;
    bool initialized;
    
    UNIFORM_BOX_PARTITION()
        :initialized(false)
    {}

    ~UNIFORM_BOX_PARTITION()
    {for(int i=1;i<=cells.counts.x;i++) for(int j=1;j<=cells.counts.y;j++) for(int ij=1;ij<=cells.counts.z;ij++) delete cells(i,j,ij);}

    void Initialize(ARRAY<PAIR<RANGE<TV>,DATA_T> >& boxes_input,const T thickness_over_two=1e-6)
    {if(boxes_input.m==0){grid.Initialize(2,2,2,RANGE<TV>(0,1,0,1,0,1));cells.Resize(grid.Domain_Indices());return;}
    bounding_box=boxes_input(1).x;for(int k=2;k<=boxes_input.m;k++) bounding_box.Enlarge_To_Include_Box(boxes_input(k).x.Thickened(thickness_over_two));
    bounding_box=bounding_box.Thickened(thickness_over_two);
    VECTOR<T,3> lengths=bounding_box.Edge_Lengths();
    VECTOR<int,3> dimensions;
    int max_axis=lengths.Dominant_Axis();T max_dimension=((T)3*pow((T)boxes_input.m,(T)one_third));
    dimensions[max_axis]=(int)max_dimension; 
    int other_axis_1=max_axis%3+1;int other_axis_2=other_axis_1%3+1; // rotate to get other axis indices
    dimensions[other_axis_1]=(int)(max_dimension*T(lengths[other_axis_1])/T(lengths[max_axis]));
    dimensions[other_axis_2]=(int)(max_dimension*T(lengths[other_axis_2])/T(lengths[max_axis]));
    dimensions=clamp_min(dimensions,VECTOR<int,3>(1,1,1));
    grid.Initialize(dimensions[1]+1,dimensions[2]+1,dimensions[3]+1,bounding_box);
    if(initialized) for(int i=1;i<=cells.counts.x;i++) for(int j=1;j<=cells.counts.y;j++) for(int ij=1;ij<=cells.counts.z;ij++) delete cells(i,j,ij);
    cells.Resize(grid.Get_MAC_Grid().Domain_Indices(),false,false);cells.Fill(0);initialized=true;
    for(int k=1;k<=boxes_input.m;k++){
        RANGE<TV_INT> domain;
        grid.Cell(boxes_input(k).x.Minimum_Corner(),domain.min_corner,0);
        grid.Cell(boxes_input(k).x.Maximum_Corner(),domain.max_corner,0);
        domain.max_corner=TV_INT::Componentwise_Min(domain.max_corner,grid.numbers_of_cells);
        for(int i=domain.min_corner.x;i<=domain.max_corner.x;i++)for(int j=domain.min_corner.y;j<=domain.max_corner.y;j++)for(int ij=domain.min_corner.z;ij<=domain.max_corner.z;ij++){
            if(!cells(i,j,ij))cells(i,j,ij)=new ARRAY<DATA_T>;
            cells(i,j,ij)->Append(boxes_input(k).y);}}}

//#####################################################################
// Function Map_Intersection
//#####################################################################
template<class HELPER_T>
bool Map_Intersection(RAY<VECTOR<T,3> >& ray,HELPER_T pointer) const
{
    T t_start=0;
    if(bounding_box.Lazy_Outside(ray.endpoint)){RAY<VECTOR<T,3> > ray_temp=ray;if(INTERSECTION::Intersects(ray_temp,bounding_box)) t_start=ray_temp.t_max;else return false;}
    VECTOR<T,3> point=ray.Point(t_start);
    VECTOR<int,3> index=grid.Clamped_Index_End_Minus_One(point);
    VECTOR<T,3> cross,dt;VECTOR<int,3> step,end;T parallel_tolerance=(T)1e-6*grid.dX.Min();
    if(abs(ray.direction.x) < parallel_tolerance){cross.x=FLT_MAX;dt.x=0;end.x=0;}
    else{
        T one_over_direction_x=1/ray.direction.x;
        if(ray.direction.x > 0){cross.x=t_start+(grid.Axis_X(index.x+1,1)-point.x)*one_over_direction_x;dt.x=grid.dX.x*one_over_direction_x;step.x=1;end.x=grid.counts.x;}
        else{cross.x=t_start+(grid.Axis_X(index.x,1)-point.x)*one_over_direction_x;dt.x=-grid.dX.x*one_over_direction_x;step.x=-1;end.x=0;}}
    if(abs(ray.direction.y) < parallel_tolerance){cross.y=FLT_MAX;dt.y=0;end.y=0;}
    else{
        T one_over_direction_y=1/ray.direction.y;
        if(ray.direction.y > 0){cross.y=t_start+(grid.Axis_X(index.y+1,2)-point.y)*one_over_direction_y;dt.y=grid.dX.y*one_over_direction_y;step.y=1;end.y=grid.counts.y;}
        else{cross.y=t_start+(grid.Axis_X(index.y,2)-point.y)*one_over_direction_y;dt.y=-grid.dX.y*one_over_direction_y;step.y=-1;end.y=0;}}
    if(abs(ray.direction.z) < parallel_tolerance){cross.z=FLT_MAX;dt.z=0;end.z=0;}
    else{
        T one_over_direction_z=1/ray.direction.z;
        if(ray.direction.z > 0){cross.z=t_start+(grid.Axis_X(index.z+1,3)-point.z)*one_over_direction_z;dt.z=grid.dX.z*one_over_direction_z;step.z=1;end.z=grid.counts.z;}
        else{cross.z=t_start+(grid.Axis_X(index.z,3)-point.z)*one_over_direction_z;dt.z=-grid.dX.z*one_over_direction_z;step.z=-1;end.z=0;}}
    bool intersection=false;
    for(;;){
        int compare_bits=((cross[1]<cross[2])<<2)+((cross[1]<cross[3])<<1)+((cross[2]<cross[3]));
        const int compare_bits_to_axis[8]={3,2,3,2,3,3,1,1};int axis=compare_bits_to_axis[compare_bits];
        if(cells(index)&&pointer->Callback(ray,*cells(index),cross[axis]))return true;
        if(!ray.semi_infinite&&ray.t_max<cross[axis])break;
        index[axis]+=step[axis];
        if(index[axis]==end[axis])break;
        cross[axis]+=dt[axis];}
    return intersection;
}
//#####################################################################
// Function Map_Inside
//#####################################################################
template<class HELPER_T>
void Map_Inside(const VECTOR<T,3>& location,HELPER_T pointer) const
{
    int i,j,ij;grid.Cell(location,i,j,ij,0);if(cells(i,j,ij)) pointer->Callback(location,*cells(i,j,ij));
}
//#####################################################################
// Function Map_Inside
//#####################################################################
template<class HELPER_T>
void Map_Inside(const RANGE<VECTOR<T,3> >& box,HELPER_T pointer) const
{
    int i_min,j_min,ij_min,i_max,j_max,ij_max;
    grid.Cell(box.Minimum_Corner(),i_min,j_min,ij_min,0);grid.Cell(box.Maximum_Corner(),i_max,j_max,ij_max,0);
    for(int i=i_min;i<=i_max;i++)for(int j=j_min;j<=j_max;j++)for(int ij=ij_min;ij<=ij_max;ij++) if(cells(i,j,ij)) pointer->Callback(box,*cells(i,j,ij));
}
//#####################################################################
};
}
#endif
