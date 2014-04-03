//#####################################################################
// Copyright 2009, Jon Gretarsson.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Namespace SIGNED_DISTANCE
//##################################################################### 
#include <PhysBAM_Tools/Grids_Uniform/UNIFORM_GRID_ITERATOR_NODE.h>
#include <PhysBAM_Tools/Log/LOG.h>
#include <PhysBAM_Geometry/Basic_Geometry/BOX.h>
#include <PhysBAM_Geometry/Grids_Uniform_Computations/TRIANGULATED_SURFACE_SIGNED_DISTANCE_UNIFORM.h>
#include <PhysBAM_Geometry/Grids_Uniform_Level_Sets/LEVELSET_3D.h>
#include <PhysBAM_Geometry/Spatial_Acceleration/TRIANGLE_HIERARCHY.h>
#include <PhysBAM_Geometry/Topology_Based_Geometry/TRIANGULATED_SURFACE.h>
namespace PhysBAM{

namespace SIGNED_DISTANCE{
//#####################################################################
// Function Calculate
//#####################################################################
template<class T> void Calculate(TRIANGULATED_SURFACE<T>& surface,const GRID<VECTOR<T,3> >& grid,ARRAY<T,VECTOR<int,3> >& phi,bool print_progress)
{
    // initialize acceleration structures
    bool triangle_list_defined=surface.triangle_list!=0;if(!triangle_list_defined) surface.Update_Triangle_List();
    bool hierarchy_defined=surface.hierarchy!=0;if(!hierarchy_defined) surface.Initialize_Hierarchy();
    bool bounding_box_defined=surface.bounding_box!=0;if(!bounding_box_defined) surface.Update_Bounding_Box();
    bool adjacent_elements_defined=surface.mesh.adjacent_elements!=0;if(!adjacent_elements_defined) surface.mesh.Initialize_Adjacent_Elements();
    bool incident_elements_defined=surface.mesh.incident_elements!=0;if(!incident_elements_defined) surface.mesh.Initialize_Incident_Elements();

    // ensure the phi array size matches the grid
    phi.Resize(grid.Domain_Indices(),false,false);

    T epsilon=(T)1e-8*grid.min_dX;
    int total_cells=grid.counts.x*grid.counts.y*grid.counts.z,cells_done=0,progress=-1;
    for(int i=1;i<=grid.counts.x;i++) for(int j=1;j<=grid.counts.y;j++) for(int k=1;k<=grid.counts.z;k++){
        phi(i,j,k)=surface.Calculate_Signed_Distance(grid.X(i,j,k),epsilon);
        if(print_progress){
            cells_done++;int new_progress=(int)((T)100*cells_done/total_cells);
            if(new_progress > progress){
#ifndef COMPILE_WITHOUT_READ_WRITE_SUPPORT
                LOG::cout<<new_progress<<"% "<<std::flush;
#endif
                progress=new_progress;}}}
#ifndef COMPILE_WITHOUT_READ_WRITE_SUPPORT
    if(print_progress) LOG::cout<<std::endl;
#endif

    // delete acceleration structures if defined in this function
    if(!incident_elements_defined){delete surface.mesh.incident_elements;surface.mesh.incident_elements=0;}
    if(!adjacent_elements_defined){delete surface.mesh.adjacent_elements;surface.mesh.adjacent_elements=0;}
    if(!bounding_box_defined){delete surface.bounding_box;surface.bounding_box=0;}
    if(!hierarchy_defined){delete surface.hierarchy;surface.hierarchy=0;}
    if(!triangle_list_defined){delete surface.triangle_list;surface.triangle_list=0;}

}
//##################################################################### 
template void Calculate(TRIANGULATED_SURFACE<float>&,const GRID<VECTOR<float,3> >&,ARRAY<float,VECTOR<int,3> >&,bool);
#ifndef COMPILE_WITHOUT_DOUBLE_SUPPORT
template void Calculate(TRIANGULATED_SURFACE<double>&,const GRID<VECTOR<double,3> >&,ARRAY<double,VECTOR<int,3> >&,bool);
#endif
};
};
