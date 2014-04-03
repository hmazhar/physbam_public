//#####################################################################
// Copyright 2009.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#ifndef __TRIANGLE_SUBDIVISION_REFINE__
#define __TRIANGLE_SUBDIVISION_REFINE__

#include <PhysBAM_Tools/Math_Tools/cyclic_shift.h>
#include <PhysBAM_Tools/Vectors/VECTOR_3D.h>
#include <PhysBAM_Geometry/Topology/TRIANGLE_MESH.h>
#include <PhysBAM_Geometry/Topology_Based_Geometry/TRIANGLE_SUBDIVISION.h>
namespace PhysBAM{

namespace TOPOLOGY_BASED_GEOMETRY_COMPUTATIONS
{
//#####################################################################
// Function Refine_Mesh
//#####################################################################
void
Refine_Mesh(TRIANGLE_SUBDIVISION& ts,TRIANGLE_MESH& refined_triangle_mesh,const int start_index_for_new_nodes_input)
{
    ts.start_index_for_new_nodes=start_index_for_new_nodes_input?start_index_for_new_nodes_input:ts.triangle_mesh.number_nodes+1;
    if(!ts.triangle_mesh.segment_mesh) ts.delete_segment_mesh=true; // make it below, but delete it when this class is deleted
    bool triangle_edges_defined=ts.triangle_mesh.element_edges!=0;
    if(!triangle_edges_defined)ts.triangle_mesh.Initialize_Element_Edges(); // this makes a segment mesh as well

    refined_triangle_mesh.elements.Exact_Resize(4*ts.triangle_mesh.elements.m);
    
    int new_t=0;
    for(int t=1;t<=ts.triangle_mesh.elements.m;t++){
        int i,j,k;ts.triangle_mesh.elements(t).Get(i,j,k);
        int edge_ij,edge_jk,edge_ki;(*ts.triangle_mesh.element_edges)(t).Get(edge_ij,edge_jk,edge_ki);
        int ij=ts.start_index_for_new_nodes-1+edge_ij,jk=ts.start_index_for_new_nodes-1+edge_jk,ki=ts.start_index_for_new_nodes-1+edge_ki;
        refined_triangle_mesh.elements(++new_t).Set(i,ij,ki);
        refined_triangle_mesh.elements(++new_t).Set(j,jk,ij);
        refined_triangle_mesh.elements(++new_t).Set(k,ki,jk);
        refined_triangle_mesh.elements(++new_t).Set(jk,ki,ij);
        refined_triangle_mesh.number_nodes=max(refined_triangle_mesh.number_nodes,i,j,k,ij,jk,ki);} // update the number of nodes
    
    if(!triangle_edges_defined){delete ts.triangle_mesh.element_edges;ts.triangle_mesh.element_edges=0;} 
    // we keep segment_mesh for later use, i.e. it's not deleted!
}
//#####################################################################
// Function Refine_Mesh_Dual
//#####################################################################
void
Refine_Mesh_Dual(TRIANGLE_SUBDIVISION& ts,TRIANGLE_MESH& refined_triangle_mesh,const int start_index_for_new_nodes_input)
{
    ts.start_index_for_new_nodes=start_index_for_new_nodes_input?start_index_for_new_nodes_input:ts.triangle_mesh.number_nodes+1;
    bool adjacent_elements_defined=ts.triangle_mesh.adjacent_elements!=0;if(!adjacent_elements_defined)ts.triangle_mesh.Initialize_Adjacent_Elements();

    refined_triangle_mesh.elements.Exact_Resize(3*ts.triangle_mesh.elements.m);
    
    int new_t=0;
    for(int t=1;t<=ts.triangle_mesh.elements.m;t++){
        assert((*ts.triangle_mesh.adjacent_elements)(t).m==3);  // boundary subdivision unimplemented
        int tv=ts.start_index_for_new_nodes-1+t;
        int ti,tj,tk;ts.triangle_mesh.elements(t).Get(ti,tj,tk);
        for(int a=1;a<=3;a++){
            int s=(*ts.triangle_mesh.adjacent_elements)(t)(a);if(s<t)continue;
            int si,sj,sk;ts.triangle_mesh.elements(s).Get(si,sj,sk);
            if(ti==si||ti==sj||ti==sk){cyclic_shift(ti,tj,tk);if(ti==si||ti==sj||ti==sk)cyclic_shift(ti,tj,tk);}
            int sv=ts.start_index_for_new_nodes-1+s;
            refined_triangle_mesh.elements(++new_t).Set(tv,tj,sv);
            refined_triangle_mesh.elements(++new_t).Set(tv,sv,tk);}}
    refined_triangle_mesh.number_nodes=ts.triangle_mesh.number_nodes+ts.triangle_mesh.elements.m;
    
    if(!adjacent_elements_defined){delete ts.triangle_mesh.adjacent_elements;ts.triangle_mesh.adjacent_elements=0;} 
}

}
}
#endif
