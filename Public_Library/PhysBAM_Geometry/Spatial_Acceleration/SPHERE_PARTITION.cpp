//#####################################################################
// Copyright 2002-2005, Ronald Fedkiw, Geoffrey Irving.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#include <PhysBAM_Geometry/Spatial_Acceleration/SPHERE_PARTITION.h>
using namespace PhysBAM;
//#####################################################################
// Function Normal
//#####################################################################
template<class T> VECTOR<T,3> SPHERE_PARTITION<T>::
Normal(const VECTOR<T,3>& location,const int aggregate) const
{
    return spheres(aggregate).Normal(location); 
}
//#####################################################################
// Function Inside
//#####################################################################
template<class T> bool SPHERE_PARTITION<T>::
Inside(const VECTOR<T,3>& location,const T thickness_over_two) const 
{
    if(box.Outside(location,thickness_over_two)) return false; 
    
    // find the voxel that contains the location - left borders
    int i,j,ij;Find_Voxel(location,i,j,ij);
    // check all the spheres in that voxel
    if(voxel_sphere_list(i,j,ij)) for(int k=1;k<=voxel_sphere_list(i,j,ij)->m;k++) if(spheres((*voxel_sphere_list(i,j,ij))(k)).Inside(location,thickness_over_two)) return true;

    return false; // otherwise
}
//#####################################################################
// Function Outside
//#####################################################################
template<class T> bool SPHERE_PARTITION<T>::
Outside(const VECTOR<T,3>& location,const T thickness_over_two) const  
{
    if(box.Outside(location,thickness_over_two)) return true; 
    
    // find the voxel that contains the location - left borders
    int i,j,ij;Find_Voxel(location,i,j,ij);
    // check all the spheres in that voxel
    if(voxel_sphere_list(i,j,ij)) for(int k=1;k<=voxel_sphere_list(i,j,ij)->m;k++) if(!spheres((*voxel_sphere_list(i,j,ij))(k)).Outside(location,thickness_over_two)) return false;

    return true; // otherwise
}
//#####################################################################
// Function Boundary
//#####################################################################
template<class T> bool SPHERE_PARTITION<T>::
Boundary(const VECTOR<T,3>& location,const T thickness_over_two) const
{
    return !Inside(location,thickness_over_two) && !Outside(location,thickness_over_two);
}
//#####################################################################
// Function Set_Up_Grid
//#####################################################################
template<class T> void SPHERE_PARTITION<T>::
Set_Up_Grid(const int m,const int n,const int mn)
{
    box=RANGE<TV>::Empty_Box();
    for(int k=1;k<=spheres.m;k++)box.Enlarge_To_Include_Box(spheres(k).Bounding_Box());
    grid.Initialize(m,n,mn,box);
    
    for(int i=1;i<=voxel_sphere_list.counts.x;i++) for(int j=1;j<=voxel_sphere_list.counts.y;j++) for(int ij=1;ij<=voxel_sphere_list.counts.z;ij++){
        delete voxel_sphere_list(i,j,ij);voxel_sphere_list(i,j,ij)=0;}
    voxel_sphere_list.Resize(1,grid.counts.x-1,1,grid.counts.y-1,1,grid.counts.z-1);
    for(int k=1;k<=spheres.m;k++){
        // find the grid cell that contains the sphere - left borders
        int i_left,j_bottom,ij_front;Find_Voxel(spheres(k).center,i_left,j_bottom,ij_front);
        // find the starting and stopping values for the sphere extent
        T radius=spheres(k).radius;
        int i_range=(int)(radius/grid.dX.x+1),j_range=(int)(radius/grid.dX.y+1),ij_range=(int)(radius/grid.dX.z+1);
        int i_start=max(1,i_left-i_range),i_end=min(grid.counts.x-1,i_left+i_range),j_start=max(1,j_bottom-j_range),
             j_end=min(grid.counts.y-1,j_bottom+j_range),ij_start=max(1,ij_front-ij_range),ij_end=min(grid.counts.z-1,ij_front+ij_range);
        // map out the cells that get the sphere
        for(int i=i_start;i<=i_end;i++) for(int j=j_start;j<=j_end;j++) for(int ij=ij_start;ij<=ij_end;ij++){
            if(!voxel_sphere_list(i,j,ij)) voxel_sphere_list(i,j,ij)=new ARRAY<int>();
            voxel_sphere_list(i,j,ij)->Append(k);}}
}
//#####################################################################
// Function Find_Voxel
//#####################################################################
template<class T> void SPHERE_PARTITION<T>::
Find_Voxel(const VECTOR<T,3>& location,int& i_left,int& j_bottom,int& ij_front) const   
{
    typedef VECTOR<int,3> TV_INT;
    TV_INT index=clamp(TV_INT((location-grid.X(TV_INT::All_Ones_Vector()))/grid.dX)+1,TV_INT::All_Ones_Vector(),grid.counts-1);
    i_left=index.x;j_bottom=index.y;ij_front=index.z;
}
//#####################################################################
template class SPHERE_PARTITION<float>;
#ifndef COMPILE_WITHOUT_DOUBLE_SUPPORT
template class SPHERE_PARTITION<double>;
#endif
