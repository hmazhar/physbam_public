//#####################################################################
// Copyright 2004-2005, Zhaosheng Bao, Frank Losasso, Sergey Koltakov, Andrew Selle.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class DUALCONTOUR_OCTREE
//##################################################################### 
#ifndef COMPILE_WITHOUT_DYADIC_SUPPORT
#include <PhysBAM_Tools/Grids_Dyadic/DYADIC_GRID_ITERATOR_CELL.h>
#include <PhysBAM_Tools/Grids_Dyadic/MAP_OCTREE_MESH.h>
#include <PhysBAM_Tools/Grids_Dyadic/OCTREE_GRID.h>
#include <PhysBAM_Tools/Grids_Dyadic_Interpolation/INTERPOLATION_DYADIC.h>
#include <PhysBAM_Tools/Grids_Dyadic_Interpolation/LINEAR_INTERPOLATION_DYADIC_HELPER.h>
#include <PhysBAM_Tools/Grids_Uniform/UNIFORM_GRID_ITERATOR_CELL.h>
#include <PhysBAM_Geometry/Grids_Dyadic_Computations/DUALCONTOUR_OCTREE.h>
#include <PhysBAM_Geometry/Grids_Dyadic_Level_Sets/LEVELSET_OCTREE.h>
#include <PhysBAM_Geometry/Topology_Based_Geometry/TRIANGULATED_SURFACE.h>
using namespace PhysBAM;
//#####################################################################
// Function Edge_Process
//#####################################################################
template <class T> static void
Edge_Process(void* data,const OCTREE_CELL<T>* cell1,const OCTREE_CELL<T>* cell2,const OCTREE_CELL<T>* cell3,const OCTREE_CELL<T>* cell4,int orient)
{
    DUALCONTOUR_OCTREE<T>* this_ptr=(DUALCONTOUR_OCTREE<T>*)data;
    if(cell1 == 0 || cell2 == 0 || cell3 == 0 || cell4 == 0)return;
    const OCTREE_CELL<T>* cells[4]={cell1,cell2,cell3,cell4};
    int deepest_depth;int deepest_tree;MAP_OCTREE_MESH<T>::Get_Deepest_Cell(4,cells,deepest_tree, deepest_depth);
    int node_table[3][4][2]={{{6,7},{2,3},{0,1},{4,5}},{{5,7},{4,6},{0,2},{1,3}},{{3,7},{1,5},{0,4},{2,6}}};

    int node_index1=cells[deepest_tree]->Node(node_table[orient][deepest_tree][0]),node_index2=cells[deepest_tree]->Node(node_table[orient][deepest_tree][1]);
    T phi1=this_ptr->phi_nodes(node_index1),phi2=this_ptr->phi_nodes(node_index2);
    if(this_ptr->topology.m <= this_ptr->num_topology+2)this_ptr->topology.Resize(this_ptr->topology.m+2500);
    if(phi1<=0&&phi2>0)this_ptr->topology(++this_ptr->num_topology).Set(cell4->Cell(),cell3->Cell(),cell2->Cell(),cell1->Cell());
    else if(phi1>0&&phi2<=0)this_ptr->topology(++this_ptr->num_topology).Set(cell1->Cell(),cell2->Cell(),cell3->Cell(),cell4->Cell());
}
//#####################################################################
// Function Dualcontour_Octree
//#####################################################################
template <class T> void DUALCONTOUR_OCTREE<T>::
Dualcontour_Octree()
{
    geometry.Resize(levelset->grid.number_of_cells);
    normals.Resize(levelset->grid.number_of_cells);levelset->Compute_Normals();
    phi_nodes.Resize(levelset->grid.number_of_nodes,false);
    LINEAR_INTERPOLATION_DYADIC_HELPER<OCTREE_GRID<T> >::Interpolate_From_Cells_To_Nodes(levelset->grid,levelset->phi,phi_nodes);

    for(DYADIC_GRID_ITERATOR_CELL<OCTREE_GRID<T> > iterator(levelset->grid,0);iterator.Valid();iterator.Next()){
        int cell_number=iterator.Cell_Index();OCTREE_CELL<T>* cell=iterator.Cell_Pointer();
        TV position=iterator.Location();
        TV normal=levelset->Normal_From_Leaf_Cell(cell,position,&phi_nodes);
        T phi=levelset->Phi_From_Close_Cell(cell,position,&phi_nodes);
        int iterations=0;
        while(abs(phi)>1e-5 && (iterations++)<10){
            position-=phi*normal;
            phi=levelset->Phi_From_Close_Cell(cell,position,&phi_nodes);
            normal=levelset->Extended_Normal(position,true,&phi_nodes);}
        // clamp vertices (analogous to Ensure_Vertices_In_Correct_Cells in DUALCONTOUR_3D)
        position=cell->Bounding_Box().Clamp(position);
        geometry(cell_number)=position;normals(cell_number)=normal;}

    MAP_OCTREE_MESH<T>::Map_Edges(levelset->grid.uniform_grid,levelset->grid.cells,0,this,Edge_Process);
    topology.Resize(num_topology);
}
//#####################################################################
// Function Get_Triangulated_Surface
//#####################################################################
template <class T> TRIANGULATED_SURFACE<T>* DUALCONTOUR_OCTREE<T>::
Get_Triangulated_Surface()
{
    TRIANGULATED_SURFACE<T>* surface=TRIANGULATED_SURFACE<T>::Create();

    GEOMETRY_PARTICLES<TV>& particles=surface->particles;
    TRIANGLE_MESH& mesh=surface->mesh;
    mesh.elements.Resize(num_topology*2);

    particles.array_collection->Add_Elements(geometry.m);
    particles.X=geometry;
    ARRAY<TV>* vertex_normals=new ARRAY<TV>(normals);
    surface->Update_Number_Nodes();
    
    int current_triangle=0;
    for(int i=1;i<=num_topology;i++){
        int p1,p2,p3,p4;topology(i).Get(p1,p2,p3,p4);
        TV x1=particles.X(p1),x2=particles.X(p2),x3=particles.X(p3),x4=particles.X(p4);
        if((TV::Cross_Product(x1-x2,x1-x4)).Magnitude_Squared()>0) mesh.elements(++current_triangle).Set(p1,p2,p4);
        if((TV::Cross_Product(x2-x3,x2-x4)).Magnitude_Squared()>0) mesh.elements(++current_triangle).Set(p2,p3,p4);}
    mesh.elements.Resize(current_triangle);
    surface->Update_Triangle_List();surface->Use_Vertex_Normals();surface->vertex_normals=vertex_normals;
    return surface;
}
//#####################################################################
template class DUALCONTOUR_OCTREE<float>;
#ifndef COMPILE_WITHOUT_DOUBLE_SUPPORT
template class DUALCONTOUR_OCTREE<double>;
#endif
#endif
