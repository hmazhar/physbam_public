//#####################################################################
// Copyright 2006, Geoffrey Irving, Tamar Shinar, Eftychios Sifakis.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class EMBEDDING_POLICY
//#####################################################################
#ifndef __EMBEDDING_POLICY__
#define __EMBEDDING_POLICY__

namespace PhysBAM{

template<class TV,int d> class EMBEDDED_OBJECT;
template<class TV,int d> class EMBEDDED_MATERIAL_SURFACE;
template<class TV> class EMBEDDED_TRIANGULATED_OBJECT;
template<class TV> class TRIANGLES_OF_MATERIAL;
template<class T> class EMBEDDED_TETRAHEDRALIZED_VOLUME;
template<class T> class EMBEDDED_TETRAHEDRALIZED_VOLUME_BOUNDARY_SURFACE;
template<class T,int d> class VECTOR;

template<class TV,int d> struct EMBEDDING_POLICY;

//#####################################################################
// 2D
//#####################################################################
template<class T>
struct EMBEDDING_POLICY<VECTOR<T,2>,2>
{
    typedef EMBEDDED_TRIANGULATED_OBJECT<VECTOR<T,2> > EMBEDDED_OBJECT;
    typedef TRIANGLES_OF_MATERIAL<VECTOR<T,2> > EMBEDDING;
    typedef PhysBAM::EMBEDDED_MATERIAL_SURFACE<VECTOR<T,2>,2> EMBEDDED_MATERIAL_SURFACE;
};
//#####################################################################
// 3D
//#####################################################################
template<class T>
struct EMBEDDING_POLICY<VECTOR<T,3>,2>
{
    typedef EMBEDDED_TRIANGULATED_OBJECT<VECTOR<T,3> > EMBEDDED_OBJECT;
    typedef TRIANGLES_OF_MATERIAL<VECTOR<T,3> > EMBEDDING;
    typedef PhysBAM::EMBEDDED_MATERIAL_SURFACE<VECTOR<T,3>,2> EMBEDDED_MATERIAL_SURFACE;
};
template<class T>
struct EMBEDDING_POLICY<VECTOR<T,3>,3>
{
    typedef EMBEDDED_TETRAHEDRALIZED_VOLUME<T> EMBEDDED_OBJECT;
    typedef EMBEDDED_TETRAHEDRALIZED_VOLUME_BOUNDARY_SURFACE<T> EMBEDDING;
    typedef PhysBAM::EMBEDDED_MATERIAL_SURFACE<VECTOR<T,3>,3> EMBEDDED_MATERIAL_SURFACE;
};
//#####################################################################
}
#endif
