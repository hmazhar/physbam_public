//#####################################################################
// Copyright 2009, Jon Gretarsson.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Namespace SIGNED_DISTANCE
//##################################################################### 
#ifndef COMPILE_WITHOUT_DYADIC_SUPPORT
#include <PhysBAM_Tools/Grids_Dyadic/DYADIC_GRID_ITERATOR_NODE.h>
#include <PhysBAM_Tools/Grids_Dyadic/OCTREE_GRID.h>
#include <PhysBAM_Tools/Grids_Uniform/UNIFORM_GRID_ITERATOR_NODE.h>
#include <PhysBAM_Tools/Log/LOG.h>
#include <PhysBAM_Geometry/Grids_Dyadic_Computations/TRIANGULATED_SURFACE_SIGNED_DISTANCE_DYADIC.h>
#include <PhysBAM_Geometry/Grids_Dyadic_Level_Sets/LEVELSET_OCTREE.h>
#include <PhysBAM_Geometry/Grids_Uniform_Level_Sets/LEVELSET_3D.h>
#include <PhysBAM_Geometry/Spatial_Acceleration/TRIANGLE_HIERARCHY.h>
#include <PhysBAM_Geometry/Topology_Based_Geometry/TRIANGULATED_SURFACE.h>
namespace PhysBAM{

namespace SIGNED_DISTANCE{
//#####################################################################
// Function Calculate
//#####################################################################
template<class T> void Calculate(TRIANGULATED_SURFACE<T>& surface,const GRID<VECTOR<T,3> >& uniform_grid,OCTREE_GRID<T>& grid,ARRAY<T>& phi,const int maximum_depth,const T half_bandwidth,bool verbose)
{
    // initialize octree grid
    grid.Initialize(uniform_grid,maximum_depth,3,true,false);
    LEVELSET_OCTREE<T> octree_levelset(grid,phi);
    T epsilon=(T)1e-8*grid.Minimum_Edge_Length();
    
    // initialize acceleration structures
    bool triangle_list_defined=surface.triangle_list!=0;if(!triangle_list_defined) surface.Update_Triangle_List();
    bool hierarchy_defined=surface.hierarchy!=0;if(!hierarchy_defined) surface.Initialize_Hierarchy();
    bool adjacent_elements_defined=surface.mesh.adjacent_elements!=0;if(!adjacent_elements_defined) surface.mesh.Initialize_Adjacent_Elements();
    bool incident_elements_defined=surface.mesh.incident_elements!=0;if(!incident_elements_defined) surface.mesh.Initialize_Incident_Elements();

    // put all cells in a list
    ARRAY<OCTREE_CELL<T>*> live_cells(uniform_grid.numbers_of_cells.Product(),false),new_cells;
    int index=0;for(int i=1;i<uniform_grid.counts.x;i++)for(int j=1;j<uniform_grid.counts.y;j++)for(int ij=1;ij<uniform_grid.counts.z;ij++) live_cells(++index)=grid.cells(i,j,ij);
    int phi_index=0;

    // generate levelset
    for(int depth=1;;depth++){
#ifndef COMPILE_WITHOUT_READ_WRITE_SUPPORT
        if(verbose) LOG::cout<<"octree depth = "<<depth<<", cells =  "<<grid.number_of_cells<<" cells"<<std::endl;
#endif
        // compute phi values
        phi.Resize(grid.number_of_nodes);
        while(phi_index<grid.number_of_nodes){phi_index++;phi(phi_index)=surface.Calculate_Signed_Distance(grid.Node_Location(phi_index),epsilon);}
        // refine cells
        if(depth==maximum_depth) break;
        new_cells.Remove_All();new_cells.Preallocate(live_cells.m);
        for(int c=1;c<=live_cells.m;c++)if(octree_levelset.Whitney_Criteria(live_cells(c),half_bandwidth))
            live_cells(c)->Create_Children(grid.number_of_cells,&new_cells,grid.number_of_nodes,0,grid.number_of_faces,0,&grid);
        grid.Tree_Topology_Changed();
        ARRAY<OCTREE_CELL<T>*>::Exchange_Arrays(live_cells,new_cells);}
#ifndef COMPILE_WITHOUT_READ_WRITE_SUPPORT
    if(verbose) LOG::cout<<"total cells = "<<grid.number_of_cells<<" out of "<<(1<<(maximum_depth-1))<<"^3 = "<<(1<<3*(maximum_depth-1))<<" possible"<<std::endl;
#endif

    // delete acceleration structures if defined in this function
    if(!incident_elements_defined){delete surface.mesh.incident_elements;surface.mesh.incident_elements=0;}
    if(!adjacent_elements_defined){delete surface.mesh.adjacent_elements;surface.mesh.adjacent_elements=0;}
    if(!hierarchy_defined){delete surface.hierarchy;surface.hierarchy=0;}
    if(!triangle_list_defined){delete surface.triangle_list;surface.triangle_list=0;}
}
//##################################################################### 
template void Calculate(TRIANGULATED_SURFACE<float>&,const GRID<VECTOR<float,3> >&,OCTREE_GRID<float>&,ARRAY<float>&,const int,const float,bool);
#ifndef COMPILE_WITHOUT_DOUBLE_SUPPORT
template void Calculate(TRIANGULATED_SURFACE<double>&,const GRID<VECTOR<double,3> >&,OCTREE_GRID<double>&,ARRAY<double>&,const int,const double,bool);
#endif
};
};
#endif
