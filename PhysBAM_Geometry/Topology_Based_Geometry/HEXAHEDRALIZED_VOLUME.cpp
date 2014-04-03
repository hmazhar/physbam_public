//#####################################################################
// Copyright 2002-2006, Robert Bridson, Ronald Fedkiw, Geoffrey Irving, Avi Robinson-Mosher, Eftychios Sifakis, Joseph Teran.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class HEXAHEDRALIZED_VOLUME
//##################################################################### 
#include <PhysBAM_Tools/Arrays/INDIRECT_ARRAY.h>
#include <PhysBAM_Tools/Data_Structures/HASHTABLE.h>
#include <PhysBAM_Tools/Grids_Uniform/GRID.h>
#include <PhysBAM_Tools/Log/LOG.h>
#include <PhysBAM_Tools/Math_Tools/RANGE.h>
#include <PhysBAM_Geometry/Basic_Geometry/HEXAHEDRON.h>
#include <PhysBAM_Geometry/Spatial_Acceleration/BOX_HIERARCHY.h>
#include <PhysBAM_Geometry/Topology_Based_Geometry/HEXAHEDRALIZED_VOLUME.h>
#include <PhysBAM_Geometry/Topology_Based_Geometry/TETRAHEDRALIZED_VOLUME.h>
#include <PhysBAM_Geometry/Topology_Based_Geometry/TRIANGULATED_SURFACE.h>
#include <PhysBAM_Geometry/Topology_Based_Geometry_Computations/HEXAHEDRALIZED_VOLUME_REFRESH.h>
using namespace PhysBAM;
//#####################################################################
// Constructor
//#####################################################################
template<class T> HEXAHEDRALIZED_VOLUME<T>::
HEXAHEDRALIZED_VOLUME(HEXAHEDRON_MESH& mesh_input,GEOMETRY_PARTICLES<TV>& particles_input)
    :MESH_OBJECT<TV,HEXAHEDRON_MESH>(mesh_input,particles_input),hexahedron_list(0),tetrahedralized_volume(0),triangulated_surface(0),hierarchy(0)
{}
//#####################################################################
// Destructor
//#####################################################################
template<class T> HEXAHEDRALIZED_VOLUME<T>::
~HEXAHEDRALIZED_VOLUME()
{
    Clean_Memory();
}
//#####################################################################
// Function Clean_Memory
//#####################################################################
template<class T> void HEXAHEDRALIZED_VOLUME<T>::
Clean_Memory()
{
    MESH_OBJECT<TV,HEXAHEDRON_MESH>::Clean_Memory();
    delete hexahedron_list;hexahedron_list=0;
    delete tetrahedralized_volume;tetrahedralized_volume=0;
    delete triangulated_surface;triangulated_surface=0;
    delete hierarchy;hierarchy=0;
}
//#####################################################################
// Function Update_Hexahedron_List
//#####################################################################
template<class T> void HEXAHEDRALIZED_VOLUME<T>::
Update_Hexahedron_List()
{
    TOPOLOGY_BASED_GEOMETRY_COMPUTATIONS::Update_Hexahedron_List(*this);
}
//#####################################################################
// Funcion Initialize_Tetrahedralized_Volume
//#####################################################################
template<class T> void HEXAHEDRALIZED_VOLUME<T>::
Initialize_Tetrahedralized_Volume()
{
    TOPOLOGY_BASED_GEOMETRY_COMPUTATIONS::Initialize_Tetrahedralized_Volume(*this);
}
//#####################################################################
// Funcion Initialize_Cube_Mesh_And_Particles
//#####################################################################
template<class T> void HEXAHEDRALIZED_VOLUME<T>::
Initialize_Cube_Mesh_And_Particles(const GRID<TV>& grid)
{
    TOPOLOGY_BASED_GEOMETRY_COMPUTATIONS::Initialize_Cube_Mesh_And_Particles(*this,grid);
}
//#####################################################################
// Funcion Total_Volume
//#####################################################################
template<class T> T HEXAHEDRALIZED_VOLUME<T>::
Total_Volume() const
{
    T volume=0;
    for(int h=1;h<=mesh.elements.m;h++){int p1,p2,p3,p4,p5,p6,p7,p8;mesh.elements(h).Get(p1,p2,p3,p4,p5,p6,p7,p8);
        volume+=HEXAHEDRON<T>::Signed_Volume(particles.X(p1),particles.X(p2),particles.X(p3),particles.X(p4),particles.X(p5),particles.X(p6),particles.X(p7),particles.X(p8));}
    return volume;
}
//#####################################################################
// Funcion Initialize_Triangulated_Surface
//#####################################################################
template<class T> void HEXAHEDRALIZED_VOLUME<T>::
Initialize_Triangulated_Surface()
{
    mesh.Initialize_Boundary_Mesh();
    triangulated_surface=new TRIANGULATED_SURFACE<T>(*mesh.boundary_mesh,particles);
}
//#####################################################################
template class HEXAHEDRALIZED_VOLUME<float>;
#ifndef COMPILE_WITHOUT_DOUBLE_SUPPORT
template class HEXAHEDRALIZED_VOLUME<double>;
#endif
