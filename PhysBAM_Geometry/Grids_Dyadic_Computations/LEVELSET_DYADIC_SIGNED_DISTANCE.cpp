//#####################################################################
// Copyright 2009, Jon Gretarsson.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Namespace SIGNED_DISTANCE
//##################################################################### 
#ifndef COMPILE_WITHOUT_DYADIC_SUPPORT
#include <PhysBAM_Tools/Data_Structures/PAIR.h>
#include <PhysBAM_Tools/Grids_Dyadic/MAP_OCTREE_MESH.h>
#include <PhysBAM_Tools/Grids_Dyadic/OCTREE_GRID.h>
#include <PhysBAM_Tools/Grids_Dyadic/QUADTREE_GRID.h>
#include <PhysBAM_Tools/Grids_Dyadic_Arrays/GRID_ARRAYS_POLICY_DYADIC.h>
#include <PhysBAM_Tools/Grids_Uniform/UNIFORM_GRID_ITERATOR_NODE.h>
#include <PhysBAM_Tools/Grids_Uniform_Arrays/ARRAYS_ND.h>
#include <PhysBAM_Tools/Log/LOG.h>
#include <PhysBAM_Geometry/Grids_Dyadic_Computations/LEVELSET_DYADIC_SIGNED_DISTANCE.h>
#include <PhysBAM_Geometry/Grids_Dyadic_Level_Sets/LEVELSET_DYADIC.h>
#include <PhysBAM_Geometry/Grids_Dyadic_Level_Sets/LEVELSET_OCTREE.h>
#include <PhysBAM_Geometry/Grids_Uniform_Level_Sets/LEVELSET_3D.h>

namespace PhysBAM{
namespace SIGNED_DISTANCE{
//#####################################################################
// Function Calculate
//#####################################################################
template<class T_GRID,class T_UNIFORM,class T_ARRAYS> void Calculate(const LEVELSET_DYADIC<T_GRID>& levelset,const T_UNIFORM& grid,T_ARRAYS& phi,bool verbose)
{
    typedef typename T_UNIFORM::NODE_ITERATOR NODE_ITERATOR;
    for(NODE_ITERATOR iterator(grid,0);iterator.Valid();iterator.Next()) phi(iterator.Node_Index())=levelset.Phi(iterator.Location());
}
//#####################################################################
// Function Calculate
//#####################################################################
template<class T> struct FIND_CELLS_TO_REFINE_HELPER{const ARRAY<T>* phi;ARRAY<OCTREE_CELL<T>*>* cells_to_refine;const T* refinement_distance;int* maximum_tree_depth;};
template<class T> static void Calculate_Helper(void* data,const OCTREE_CELL<T>* cell)
{
    FIND_CELLS_TO_REFINE_HELPER<T>* helper=(FIND_CELLS_TO_REFINE_HELPER<T>*)data;
    if(cell->Depth_Of_This_Cell()>=*helper->maximum_tree_depth) return;
    T dx=(T).5*(T)root_three*cell->DX().x+*helper->refinement_distance;
    for(int i=0;i<8;i++) if(abs((*helper->phi)(cell->Node(i))) < dx){helper->cells_to_refine->Append((OCTREE_CELL<T>*)cell);return;}
}
template<class T_GRID,class T> void Calculate(const LEVELSET_3D<T_GRID>& levelset,OCTREE_GRID<T>& grid,ARRAY<T>& phi,bool verbose)
{
    // gets the dimensions, unifirm grid resolution and number of levelset directly from grid
    T refinement_distance=3*grid.Minimum_Edge_Length();
    ARRAY<OCTREE_CELL<T>*> cells_to_refine;
    FIND_CELLS_TO_REFINE_HELPER<T> helper;helper.phi=&phi;helper.cells_to_refine=&cells_to_refine;
    helper.refinement_distance=&refinement_distance;helper.maximum_tree_depth=&grid.maximum_depth;
    while(true){
        LEVELSET_OCTREE<T> levelset(grid,phi);
        phi.Resize(grid.number_of_nodes);
        for(int i=1;i<=grid.number_of_nodes;i++)
            phi(i)=levelset.Phi(grid.Node_Location(i));
        levelset.Fast_Marching_Method();
        cells_to_refine.Resize(0);
        MAP_OCTREE_MESH<T>::Map_Cells(grid.uniform_grid,grid.cells,0,&helper,Calculate_Helper);
        if(cells_to_refine.m==0) break;
        for(int i=1;i<=cells_to_refine.m;i++){
            assert(!cells_to_refine(i)->Has_Children());
            cells_to_refine(i)->Create_Children(grid.number_of_cells,0,grid.number_of_nodes,0,grid.number_of_faces,0,&grid);}
        grid.Tree_Topology_Changed();
#ifndef COMPILE_WITHOUT_READ_WRITE_SUPPORT
        if(verbose) LOG::cout << "total cells = " << grid.number_of_cells << std::endl;
#endif
    }
}
//#####################################################################
// Function Calculate_Based_On_Curvature
//#####################################################################
namespace{
template<class T,class T_GRID> bool Compare_Curvature_To_DX(const LEVELSET_3D<T_GRID>& levelset,OCTREE_GRID<T>& octree_grid,const OCTREE_CELL<T>& cell_input){
    PHYSBAM_NOT_IMPLEMENTED(); // needs to be updated for new octrees
    assert(levelset.curvature);
    VECTOR<T,3> sample=octree_grid.Node_Location(1,cell_input),DX=octree_grid.DX(cell_input);
    T dx_max_element=DX.Max();
    int x_samples=(int)ceil(DX.x/levelset.grid.min_dX),y_samples=(int)ceil(DX.y/levelset.grid.min_dX),z_samples=(int)ceil(DX.z/levelset.grid.min_dX);
    for(int i=1;i<=x_samples;i++){
        for(int j=1;j<=y_samples;j++){
            for(int k=1;k<=z_samples;k++){
                if(1/abs(levelset.Curvature(sample)) < dx_max_element) return true;
                sample.z+=levelset.grid.min_dX;}
            sample.y+=levelset.grid.min_dX;}
        sample.x+=levelset.grid.min_dX;}
    return false;
}
}
template<class T,class T_GRID> void Calculate_Based_On_Curvature(const LEVELSET_3D<T_GRID>& levelset,OCTREE_GRID<T>& octree_grid,ARRAY<T>& phi,const int levels,bool output_statistics)
{
    PHYSBAM_NOT_IMPLEMENTED(); // needs to be updated for new octrees
#if 0
    assert(&octree_grid && &phi);

    // initialize octree variables
    octree_grid.Initialize(levelset.grid.Domain(),levels);
    octree_grid.root_cell=*(OCTREE_CELL<T>::Create_Root(true,false,octree_grid.number_of_cells,octree_grid.number_of_nodes,octree_grid.number_of_faces,
                                            (VECTOR<float,3>)octree_grid.domain_center,(VECTOR<float,3>)octree_grid.domain_length));

    // need the curvature of the level set
    bool curvature_defined=levelset.curvature!=0;if(!curvature_defined) levelset.Compute_Curvature();

    // update phi for the original 8 nodes
    phi.Resize(1,8);
    for(int k=0;k<8;k++) phi(octree_grid.root_cell.Node(k))=levelset.Phi(octree_grid.Node_Location(k,&(octree_grid.root_cell)));

    // refine one level at a time
    LEVELSET_OCTREE<T> octree_levelset(octree_grid,phi);
    ARRAY<OCTREE_CELL<T>*> cells_to_refine(1,1),new_children;cells_to_refine(1)=&(octree_grid.root_cell);
    ARRAY<PAIR<OCTREE_CELL<T>*,int> > compute_phi;
    int level=1;
    while(cells_to_refine.m && level < levels){level++;
#ifndef COMPILE_WITHOUT_READ_WRITE_SUPPORT
        if(output_statistics) LOG::cout << "octree level = " << level << " - refining " << cells_to_refine.m << " out of " << octree_grid.number_of_cells << " cells" << std::endl;
#endif
        if(level < levels)
            octree_grid.root_cell.Add_Children(cells_to_refine,octree_grid.number_of_cells,&new_children,octree_grid.number_of_nodes,&compute_phi,octree_grid.number_of_faces,0);
        else octree_grid.root_cell.Add_Children(cells_to_refine,octree_grid.number_of_cells,0,octree_grid.number_of_nodes,&compute_phi,octree_grid.number_of_faces,0); // don't need children for the last time through
        cells_to_refine.Remove_All();phi.Resize(1,octree_grid.number_of_nodes);
        // compute phi for the new nodes
        for(int k=1;k<=compute_phi.m;k++){
            OCTREE_CELL<T>* cell=compute_phi(k).x;int node_number=compute_phi(k).y,phi_node=cell->Node(node_number);
            phi(phi_node)=levelset.Phi(octree_grid.Node_Location(node_number,cell));}
        if(level < levels){
            int index=0;for(int c=1;c<=new_children.m;c++)
                if(octree_levelset.Whitney_Criteria(new_children(c)) && Compare_Curvature_To_DX(octree_grid,new_children(c))){index++;cells_to_refine.Append(new_children(c));}
            cells_to_refine.Resize(1,index);}
        new_children.Remove_All();compute_phi.Remove_All();}
#ifndef COMPILE_WITHOUT_READ_WRITE_SUPPORT
    if(output_statistics) LOG::cout << "total cells = " << octree_grid.number_of_cells << " out of " << pow(2.,level-1) << "^3 = " << pow(pow(2.,level-1),3.) << " possible" << std::endl;
#endif

    // delete curvature if defined in this function
    if(!curvature_defined){delete levelset.curvature;levelset.curvature=0;}
#endif
}
//#####################################################################
#define LEVELSET_DYADIC_SIGNED_DISTANCE_HELPER(T) \
    template void Calculate(const LEVELSET_3D<GRID<VECTOR<T,3> > >& levelset,OCTREE_GRID<T>&,ARRAY<T>&,bool); \
    template void Calculate(const LEVELSET_DYADIC<OCTREE_GRID<T> >&,const OCTREE_GRID<T>::UNIFORM_GRID&,GRID_ARRAYS_POLICY<OCTREE_GRID<T>::UNIFORM_GRID>::ARRAYS_SCALAR&,bool); \
    template void Calculate(const LEVELSET_DYADIC<QUADTREE_GRID<T> >&,const QUADTREE_GRID<T>::UNIFORM_GRID&,GRID_ARRAYS_POLICY<QUADTREE_GRID<T>::UNIFORM_GRID>::ARRAYS_SCALAR&,bool); \
    template void Calculate_Based_On_Curvature(const LEVELSET_3D<GRID<VECTOR<T,3> > >& levelset,OCTREE_GRID<T>&,ARRAY<T>&,const int,bool);

LEVELSET_DYADIC_SIGNED_DISTANCE_HELPER(float);
#ifndef COMPILE_WITHOUT_DOUBLE_SUPPORT
LEVELSET_DYADIC_SIGNED_DISTANCE_HELPER(double);
#endif
};
};
#endif
