//#####################################################################
// Copyright 2009.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#ifndef __TRIANGULATED_SURFACE_SUBDIVISION__
#define __TRIANGULATED_SURFACE_SUBDIVISION__

#include <PhysBAM_Tools/Arrays/ARRAY.h>
#include <PhysBAM_Tools/Vectors/VECTOR_3D.h>
#include <PhysBAM_Geometry/Topology/TRIANGLE_MESH.h>
#include <PhysBAM_Geometry/Topology_Based_Geometry/TRIANGLE_SUBDIVISION.h>
#include <PhysBAM_Geometry/Topology_Based_Geometry/TRIANGULATED_SURFACE.h>

namespace PhysBAM{

namespace TOPOLOGY_BASED_GEOMETRY_COMPUTATIONS{

template<class T>
void Linearly_Subdivide(TRIANGULATED_SURFACE<T>& ts)
{
    TRIANGLE_SUBDIVISION subdivision(ts.mesh);
    TRIANGLE_MESH refined_mesh;
    subdivision.Refine_Mesh(refined_mesh);
    ARRAY<VECTOR<T,3> > X_save(ts.particles.X),V_save(ts.particles.V);
    ts.particles.array_collection->Add_Elements(refined_mesh.number_nodes-ts.particles.array_collection->Size());
    if(X_save.Size()) subdivision.Apply_Linear_Subdivision(X_save,ts.particles.X);
    if(V_save.Size()) subdivision.Apply_Linear_Subdivision(V_save,ts.particles.V);
    ts.mesh.Initialize_Mesh(refined_mesh);
}
template<class T>
void Loop_Subdivide(TRIANGULATED_SURFACE<T>& ts)
{
    TRIANGLE_SUBDIVISION subdivision(ts.mesh);
    TRIANGLE_MESH refined_mesh;
    subdivision.Refine_Mesh(refined_mesh);
    ARRAY<VECTOR<T,3> > X_save(ts.particles.X),V_save(ts.particles.V);
    ts.particles.array_collection->Add_Elements(refined_mesh.number_nodes-ts.particles.array_collection->Size());
    if(X_save.Size()) subdivision.Apply_Loop_Subdivision(X_save,ts.particles.X);
    if(V_save.Size()) subdivision.Apply_Loop_Subdivision(V_save,ts.particles.V);
    ts.mesh.Initialize_Mesh(refined_mesh);
}
template<class T>
void Root_Three_Subdivide(TRIANGULATED_SURFACE<T>& ts)
{
    TRIANGLE_SUBDIVISION subdivision(ts.mesh);
    TRIANGLE_MESH refined_mesh;
    subdivision.Refine_Mesh_Dual(refined_mesh);
    ARRAY<VECTOR<T,3> > X_save(ts.particles.X),V_save(ts.particles.V);
    ts.particles.array_collection->Add_Elements(refined_mesh.number_nodes-ts.particles.array_collection->Size());
    if(X_save.Size()) subdivision.Apply_Root_Three_Subdivision(X_save,ts.particles.X);
    if(V_save.Size()) subdivision.Apply_Root_Three_Subdivision(V_save,ts.particles.V);
    ts.mesh.Initialize_Mesh(refined_mesh);
}
}
}
#endif
