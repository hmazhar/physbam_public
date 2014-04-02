//#####################################################################
// Copyright 2009, Jon Gretarsson.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class RING_TESSELLATION
//##################################################################### 
#include <PhysBAM_Tools/Matrices/MATRIX_3X2.h>
#include <PhysBAM_Tools/Vectors/COMPLEX.h>
#include <PhysBAM_Geometry/Basic_Geometry/RING.h>
#include <PhysBAM_Geometry/Tessellation/RING_TESSELLATION.h>
#include <PhysBAM_Geometry/Topology_Based_Geometry/TRIANGULATED_SURFACE.h>

namespace PhysBAM{
namespace TESSELLATION{
//#####################################################################
// Function Generate_Triangles
//#####################################################################
template<class T> TRIANGULATED_SURFACE<T>* Generate_Triangles(const RING<T>& ring,const int n)
{
    typedef VECTOR<T,3> TV;
    TRIANGULATED_SURFACE<T>* surface=TRIANGULATED_SURFACE<T>::Create();
    surface->mesh.Initialize_Torus_Mesh(4,n);
    GEOMETRY_PARTICLES<TV>& particles=surface->particles;particles.array_collection->Add_Elements(4*n);
    MATRIX<T,3,2> radial_basis;
    radial_basis.Column(1)=ring.plane1.normal.Orthogonal_Vector();
    radial_basis.Column(2)=TV::Cross_Product(ring.plane1.normal,radial_basis.Column(1));
    for(int i=1,p=0;i<=n;++i){
        TV radial=radial_basis*COMPLEX<T>::Unit_Polar(T(2*pi/n)*i).Vector();
        particles.X(++p)=ring.plane1.x1+ring.inner_radius*radial;
        particles.X(++p)=ring.plane1.x1+ring.outer_radius*radial;
        particles.X(++p)=ring.plane2.x1+ring.outer_radius*radial;
        particles.X(++p)=ring.plane2.x1+ring.inner_radius*radial;}
    return surface;
}
//#####################################################################
template TRIANGULATED_SURFACE<float>* Generate_Triangles(const RING<float>&,const int);
#ifndef COMPILE_WITHOUT_DOUBLE_SUPPORT
template TRIANGULATED_SURFACE<double>* Generate_Triangles(const RING<double>&,const int);
#endif
}
}
