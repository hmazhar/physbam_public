//#####################################################################
// Copyright 2002-2006, Robert Bridson, Ronald Fedkiw, Eilene Hao, Geoffrey Irving, Igor Neverov, Andrew Selle.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#include <PhysBAM_Tools/Arrays/ARRAY.h>
#include <PhysBAM_Tools/Math_Tools/cube.h>
#include <PhysBAM_Tools/Math_Tools/cyclic_shift.h>
#include <PhysBAM_Tools/Random_Numbers/RANDOM_NUMBERS.h>
#include <PhysBAM_Tools/Vectors/VECTOR_3D.h>
#include <PhysBAM_Geometry/Basic_Geometry/TRIANGLE_3D.h>
#include <PhysBAM_Geometry/Topology/SEGMENT_MESH.h>
#include <PhysBAM_Geometry/Topology/TRIANGLE_MESH.h>
#include <PhysBAM_Geometry/Topology_Based_Geometry/TRIANGLE_SUBDIVISION.h>
#include <PhysBAM_Geometry/Topology_Based_Geometry_Computations/TRIANGLE_SUBDIVISION_FRACTAL.h>
#include <PhysBAM_Geometry/Topology_Based_Geometry_Computations/TRIANGLE_SUBDIVISION_LINEAR.h>
#include <PhysBAM_Geometry/Topology_Based_Geometry_Computations/TRIANGLE_SUBDIVISION_LOOP.h>
#include <PhysBAM_Geometry/Topology_Based_Geometry_Computations/TRIANGLE_SUBDIVISION_REFINE.h>
#include <PhysBAM_Geometry/Topology_Based_Geometry_Computations/TRIANGLE_SUBDIVISION_ROOT_THREE.h>
namespace PhysBAM{
//#####################################################################
// Function Clean_Memory
//#####################################################################
void TRIANGLE_SUBDIVISION::
Clean_Memory()
{
    if(delete_segment_mesh){delete triangle_mesh.segment_mesh;triangle_mesh.segment_mesh=0;}
    if(delete_neighbor_nodes){delete triangle_mesh.neighbor_nodes;triangle_mesh.neighbor_nodes=0;}
    if(delete_topologically_sorted_neighbor_nodes){delete triangle_mesh.topologically_sorted_neighbor_nodes;triangle_mesh.topologically_sorted_neighbor_nodes=0;}
    if(delete_boundary_mesh){delete triangle_mesh.boundary_mesh;triangle_mesh.boundary_mesh=0;}
    if(delete_incident_elements){delete triangle_mesh.incident_elements;triangle_mesh.incident_elements=0;}
    delete_segment_mesh=delete_neighbor_nodes=delete_topologically_sorted_neighbor_nodes=delete_boundary_mesh=delete_incident_elements=false;
}
//#####################################################################
// Function Refine_Mesh
//#####################################################################
void TRIANGLE_SUBDIVISION::
Refine_Mesh(TRIANGLE_MESH& refined_triangle_mesh,const int start_index_for_new_nodes_input)
{
    TOPOLOGY_BASED_GEOMETRY_COMPUTATIONS::Refine_Mesh(*this,refined_triangle_mesh,start_index_for_new_nodes_input);
}
//#####################################################################
// Function Refine_Mesh_Dual
//#####################################################################
void TRIANGLE_SUBDIVISION::
Refine_Mesh_Dual(TRIANGLE_MESH& refined_triangle_mesh,const int start_index_for_new_nodes_input)
{
    TOPOLOGY_BASED_GEOMETRY_COMPUTATIONS::Refine_Mesh_Dual(*this,refined_triangle_mesh,start_index_for_new_nodes_input);
}
//#####################################################################
// Function Apply_Linear_Subdivision
//#####################################################################
template<class TV> void TRIANGLE_SUBDIVISION::
Apply_Linear_Subdivision(ARRAY_VIEW<const TV> base_values,ARRAY_VIEW<TV> subdivided_values)
{
    TOPOLOGY_BASED_GEOMETRY_COMPUTATIONS::Apply_Linear_Subdivision(*this,base_values,subdivided_values);
}
//#####################################################################
// Function Apply_Fractal_Subdivision
//#####################################################################
template<class TV> void TRIANGLE_SUBDIVISION::
Apply_Fractal_Subdivision(ARRAY_VIEW<const TV> base_values,ARRAY_VIEW<TV> subdivided_values,const float power)
{
    TOPOLOGY_BASED_GEOMETRY_COMPUTATIONS::Apply_Fractal_Subdivision(*this,base_values,subdivided_values,power);
}
//#####################################################################
// Function Apply_Loop_Subdivision
//#####################################################################
template<class TV> void TRIANGLE_SUBDIVISION::
Apply_Loop_Subdivision(ARRAY_VIEW<const TV> base_values,ARRAY_VIEW<TV> subdivided_values)
{
    TOPOLOGY_BASED_GEOMETRY_COMPUTATIONS::Apply_Loop_Subdivision(*this,base_values,subdivided_values);
}
//#####################################################################
// Function Apply_Root_Three_Subdivision
//#####################################################################
template<class TV> void TRIANGLE_SUBDIVISION::
Apply_Root_Three_Subdivision(ARRAY_VIEW<const TV> base_values,ARRAY_VIEW<TV> subdivided_values)
{
    TOPOLOGY_BASED_GEOMETRY_COMPUTATIONS::Apply_Root_Three_Subdivision(*this,base_values,subdivided_values);
}
//#####################################################################
template void TRIANGLE_SUBDIVISION::Apply_Linear_Subdivision(ARRAY_VIEW<const VECTOR<float,3> > base_values,ARRAY_VIEW<VECTOR<float,3> > subdivided_values);
template void TRIANGLE_SUBDIVISION::Apply_Fractal_Subdivision(ARRAY_VIEW<const VECTOR<float,3> > base_values,ARRAY_VIEW<VECTOR<float,3> > subdivided_values,const float power);
template void TRIANGLE_SUBDIVISION::Apply_Loop_Subdivision(ARRAY_VIEW<const VECTOR<float,3> > base_values,ARRAY_VIEW<VECTOR<float,3> > subdivided_values);
template void TRIANGLE_SUBDIVISION::Apply_Root_Three_Subdivision(ARRAY_VIEW<const VECTOR<float,3> > base_values,ARRAY_VIEW<VECTOR<float,3> > subdivided_values);
#ifndef COMPILE_WITHOUT_DOUBLE_SUPPORT
template void TRIANGLE_SUBDIVISION::Apply_Linear_Subdivision(ARRAY_VIEW<const VECTOR<double,3> > base_values,ARRAY_VIEW<VECTOR<double,3> > subdivided_values);
template void TRIANGLE_SUBDIVISION::Apply_Fractal_Subdivision(ARRAY_VIEW<const VECTOR<double,3> > base_values,ARRAY_VIEW<VECTOR<double,3> > subdivided_values,const float power);
template void TRIANGLE_SUBDIVISION::Apply_Loop_Subdivision(ARRAY_VIEW<const VECTOR<double,3> > base_values,ARRAY_VIEW<VECTOR<double,3> > subdivided_values);
template void TRIANGLE_SUBDIVISION::Apply_Root_Three_Subdivision(ARRAY_VIEW<const VECTOR<double,3> > base_values,ARRAY_VIEW<VECTOR<double,3> > subdivided_values);
#endif
}
