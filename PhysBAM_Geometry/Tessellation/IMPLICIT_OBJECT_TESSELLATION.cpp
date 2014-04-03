//#####################################################################
// Copyright 2009, Jon Gretarsson.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// namespace TESSELLATION
//##################################################################### 
#include <PhysBAM_Tools/Log/LOG.h>
#include <PhysBAM_Tools/Math_Tools/RANGE.h>
#include <PhysBAM_Tools/Matrices/FRAME.h>
#include <PhysBAM_Geometry/Basic_Geometry/BOUNDED_HORIZONTAL_PLANE.h>
#include <PhysBAM_Geometry/Basic_Geometry/CYLINDER.h>
#include <PhysBAM_Geometry/Basic_Geometry/ORIENTED_BOX.h>
#include <PhysBAM_Geometry/Basic_Geometry/PLANE.h>
#include <PhysBAM_Geometry/Basic_Geometry/RING.h>
#include <PhysBAM_Geometry/Basic_Geometry/SPHERE.h>
#include <PhysBAM_Geometry/Basic_Geometry/TORUS.h>
#include <PhysBAM_Geometry/Grids_Uniform_Computations/DUALCONTOUR_3D.h>
#include <PhysBAM_Geometry/Grids_Uniform_Level_Sets/LEVELSET_3D.h>
#include <PhysBAM_Geometry/Implicit_Objects/ANALYTIC_IMPLICIT_OBJECT.h>
#include <PhysBAM_Geometry/Implicit_Objects/IMPLICIT_OBJECT.h>
#include <PhysBAM_Geometry/Implicit_Objects/IMPLICIT_OBJECT_TRANSFORMED.h>
#include <PhysBAM_Geometry/Implicit_Objects_Uniform/IMPLICIT_OBJECT_COMBINED.h>
#include <PhysBAM_Geometry/Implicit_Objects_Uniform/IMPLICIT_OBJECT_COMBINED_EULERIAN.h>
#include <PhysBAM_Geometry/Implicit_Objects_Uniform/LEVELSET_IMPLICIT_OBJECT.h>
#include <PhysBAM_Geometry/Implicit_Objects_Uniform/MULTIBODY_LEVELSET_IMPLICIT_OBJECT.h>
#include <PhysBAM_Geometry/Solids_Geometry/RIGID_GEOMETRY.h>
#include <PhysBAM_Geometry/Tessellation/BOUNDED_HORIZONTAL_PLANE_TESSELLATION.h>
#include <PhysBAM_Geometry/Tessellation/CYLINDER_TESSELLATION.h>
#include <PhysBAM_Geometry/Tessellation/IMPLICIT_OBJECT_TESSELLATION.h>
#include <PhysBAM_Geometry/Tessellation/ORIENTED_BOX_TESSELLATION.h>
#include <PhysBAM_Geometry/Tessellation/PLANE_TESSELLATION.h>
#include <PhysBAM_Geometry/Tessellation/RANGE_TESSELLATION.h>
#include <PhysBAM_Geometry/Tessellation/RING_TESSELLATION.h>
#include <PhysBAM_Geometry/Tessellation/SPHERE_TESSELLATION.h>
#include <PhysBAM_Geometry/Tessellation/TORUS_TESSELLATION.h>
#include <PhysBAM_Geometry/Topology_Based_Geometry/POINT_SIMPLICES_1D.h>
#include <PhysBAM_Geometry/Topology_Based_Geometry/SEGMENTED_CURVE_2D.h>
#include <PhysBAM_Geometry/Topology_Based_Geometry/TRIANGULATED_SURFACE.h>
#ifndef COMPILE_WITHOUT_DYADIC_SUPPORT
#include <PhysBAM_Geometry/Grids_Dyadic_Computations/DUALCONTOUR_OCTREE.h>
#include <PhysBAM_Geometry/Grids_Dyadic_Level_Sets/LEVELSET_OCTREE.h>
#include <PhysBAM_Geometry/Implicit_Objects_Dyadic/DYADIC_IMPLICIT_OBJECT.h>
#endif
#ifndef COMPILE_WITHOUT_RLE_SUPPORT
#include <PhysBAM_Geometry/Grids_RLE_Computations/DUALCONTOUR_RLE_3D.h>
#include <PhysBAM_Geometry/Grids_RLE_Level_Sets/LEVELSET_RLE.h>
#include <PhysBAM_Geometry/Implicit_Objects_RLE/RLE_IMPLICIT_OBJECT.h>
#endif

namespace PhysBAM{
namespace TESSELLATION{
//#####################################################################
// Function Generate_Triangles_Helper
//#####################################################################
namespace {
template<class T> TRIANGULATED_SURFACE<T>* Generate_Triangles_Helper(const LEVELSET_IMPLICIT_OBJECT<VECTOR<T,3> >* implicit)
{ // TODO(jontg): There should be a better implementation for this...
    TRIANGULATED_SURFACE<T>* surface=TRIANGULATED_SURFACE<T>::Create();
    implicit->levelset.Calculate_Triangulated_Surface_From_Marching_Tetrahedra(*surface);
    return surface;
}

template<class T> TRIANGULATED_SURFACE<T>* Generate_Triangles_Helper(const MULTIBODY_LEVELSET_IMPLICIT_OBJECT<VECTOR<T,3> >* implicit)
{
    ARRAY<TRIANGULATED_SURFACE<T>*> triangulated_surfaces(implicit->levelsets->m);
    for(int i=1;i<=implicit->levelsets->m;++i) triangulated_surfaces(i)=Generate_Triangles(*(*implicit->levelsets)(i));
    ARRAY<FRAME<VECTOR<T,3> > > identity_frames(implicit->levelsets->m);
    TRIANGULATED_SURFACE<T>* union_object=TRIANGULATED_SURFACE<T>::Union_Mesh_Objects_Relatively(triangulated_surfaces,identity_frames);
    triangulated_surfaces.Delete_Pointers_And_Clean_Memory();
    return union_object;
}

#ifndef COMPILE_WITHOUT_DYADIC_SUPPORT
template<class T> TRIANGULATED_SURFACE<T>* Generate_Triangles_Helper(const DYADIC_IMPLICIT_OBJECT<VECTOR<T,3> >* implicit)
{DUALCONTOUR_OCTREE<T> dual_contour(const_cast<LEVELSET_OCTREE<T>*>(&implicit->levelset));return dual_contour.Get_Triangulated_Surface();}
#endif

#ifndef COMPILE_WITHOUT_RLE_SUPPORT
template<class T> TRIANGULATED_SURFACE<T>* Generate_Triangles_Helper(const RLE_IMPLICIT_OBJECT<VECTOR<T,3> >* implicit)
{
    typedef typename RLE_GRID_POLICY<VECTOR<T,3> >::RLE_GRID T_GRID;
    return DUALCONTOUR_RLE_3D<T>::Create_Triangulated_Surface_From_Levelset(
        const_cast<LEVELSET_RLE<T_GRID>&>(implicit->levelset),implicit->levelset.grid.number_of_ghost_cells);
}
#endif

template<class T> TRIANGULATED_SURFACE<T>* Generate_Triangles_Helper(const IMPLICIT_OBJECT_COMBINED<VECTOR<T,3> >* implicit)
{
    TRIANGULATED_SURFACE<T>* triangles_1=Generate_Triangles(*implicit->implicit_object1);
    TRIANGULATED_SURFACE<T>* triangles_2=Generate_Triangles(*implicit->implicit_object2);
    assert(triangles_1->particles.array_collection->Size()==triangles_2->particles.array_collection->Size());
    for(int i=1;i<=triangles_1->particles.array_collection->Size();++i)
        triangles_1->particles.X(i)=(1-implicit->alpha)*triangles_1->particles.X(i)+implicit->alpha*triangles_2->particles.X(i);
    return triangles_1;
}

template<class T> TRIANGULATED_SURFACE<T>* Generate_Triangles_Helper(const IMPLICIT_OBJECT_COMBINED_EULERIAN<VECTOR<T,3> >* implicit)
{PHYSBAM_NOT_IMPLEMENTED();}

template<class T> TRIANGULATED_SURFACE<T>* Generate_Triangles_Helper(const IMPLICIT_OBJECT_TRANSFORMED<VECTOR<T,3>,RIGID_GEOMETRY<VECTOR<T,3> > >* implicit)
{ // TODO(jontg): Better way to templatize this?
    TRIANGULATED_SURFACE<T>* surface=Generate_Triangles(*implicit->object_space_implicit_object);
    for(int i=1;i<=surface->particles.array_collection->Size();++i) surface->particles.X(i)=implicit->World_Space_Point(surface->particles.X(i));
    return surface;
}

template<class T> TRIANGULATED_SURFACE<T>* Generate_Triangles_Helper(const IMPLICIT_OBJECT_TRANSFORMED<VECTOR<T,3>,FRAME<VECTOR<T,3> > >* implicit)
{ // TODO(jontg): Better way to templatize this?
    TRIANGULATED_SURFACE<T>* surface=Generate_Triangles(*implicit->object_space_implicit_object);
    for(int i=1;i<=surface->particles.array_collection->Size();++i) surface->particles.X(i)=implicit->World_Space_Point(surface->particles.X(i));
    return surface;
}
}
//#####################################################################
// Function Generate_Triangles
//#####################################################################
template<class T> POINT_SIMPLICES_1D<T>* Generate_Triangles(IMPLICIT_OBJECT<VECTOR<T,1> > const& implicit) {PHYSBAM_FATAL_ERROR();}
template<class T> SEGMENTED_CURVE_2D<T>* Generate_Triangles(IMPLICIT_OBJECT<VECTOR<T,2> > const& implicit) {PHYSBAM_FATAL_ERROR();}
template<class T> TRIANGULATED_SURFACE<T>* Generate_Triangles(IMPLICIT_OBJECT<VECTOR<T,3> > const& implicit_input)
{
    typedef VECTOR<T,3> TV;
    if(const LEVELSET_IMPLICIT_OBJECT<TV>* implicit=dynamic_cast<const LEVELSET_IMPLICIT_OBJECT<TV>*>(&implicit_input)) return Generate_Triangles_Helper(implicit);
    else if(const MULTIBODY_LEVELSET_IMPLICIT_OBJECT<TV>* implicit=dynamic_cast<const MULTIBODY_LEVELSET_IMPLICIT_OBJECT<TV>*>(&implicit_input)) return Generate_Triangles_Helper(implicit);
    else if(const IMPLICIT_OBJECT_COMBINED<TV>* implicit=dynamic_cast<const IMPLICIT_OBJECT_COMBINED<TV>*>(&implicit_input)) return Generate_Triangles_Helper(implicit); 
    else if(const IMPLICIT_OBJECT_COMBINED_EULERIAN<TV>* implicit=dynamic_cast<const IMPLICIT_OBJECT_COMBINED_EULERIAN<TV>*>(&implicit_input)) return Generate_Triangles_Helper(implicit); 
    else if(const IMPLICIT_OBJECT_TRANSFORMED<TV,FRAME<TV> >* implicit=dynamic_cast<const IMPLICIT_OBJECT_TRANSFORMED<TV,FRAME<TV> >*>(&implicit_input)) return Generate_Triangles_Helper(implicit); 
    else if(const IMPLICIT_OBJECT_TRANSFORMED<TV,RIGID_GEOMETRY<TV> >* implicit=dynamic_cast<const IMPLICIT_OBJECT_TRANSFORMED<TV,RIGID_GEOMETRY<TV> >*>(&implicit_input)) return Generate_Triangles_Helper(implicit); 
    else if(const ANALYTIC_IMPLICIT_OBJECT<PLANE<T> >* implicit=dynamic_cast<const ANALYTIC_IMPLICIT_OBJECT<PLANE<T> >*>(&implicit_input)) return Generate_Triangles(implicit->analytic);
    else if(const ANALYTIC_IMPLICIT_OBJECT<RING<T> >* implicit=dynamic_cast<const ANALYTIC_IMPLICIT_OBJECT<RING<T> >*>(&implicit_input)) return Generate_Triangles(implicit->analytic);
    else if(const ANALYTIC_IMPLICIT_OBJECT<SPHERE<TV> >* implicit=dynamic_cast<const ANALYTIC_IMPLICIT_OBJECT<SPHERE<TV> >*>(&implicit_input)) return Generate_Triangles(implicit->analytic);
    else if(const ANALYTIC_IMPLICIT_OBJECT<BOUNDED_HORIZONTAL_PLANE<TV> >* implicit=dynamic_cast<const ANALYTIC_IMPLICIT_OBJECT<BOUNDED_HORIZONTAL_PLANE<TV> >*>(&implicit_input)) return Generate_Triangles(implicit->analytic);
    else if(const ANALYTIC_IMPLICIT_OBJECT<ORIENTED_BOX<TV> >* implicit=dynamic_cast<const ANALYTIC_IMPLICIT_OBJECT<ORIENTED_BOX<TV> >*>(&implicit_input)) return Generate_Triangles(implicit->analytic);
    else if(const ANALYTIC_IMPLICIT_OBJECT<BOX<TV> >* implicit=dynamic_cast<const ANALYTIC_IMPLICIT_OBJECT<BOX<TV> >*>(&implicit_input)) return Generate_Triangles(implicit->analytic);
    else if(const ANALYTIC_IMPLICIT_OBJECT<TORUS<T> >* implicit=dynamic_cast<const ANALYTIC_IMPLICIT_OBJECT<TORUS<T> >*>(&implicit_input)) return Generate_Triangles(implicit->analytic);
    else if(const ANALYTIC_IMPLICIT_OBJECT<CYLINDER<T> >* implicit=dynamic_cast<const ANALYTIC_IMPLICIT_OBJECT<CYLINDER<T> >*>(&implicit_input)) return Generate_Triangles(implicit->analytic);
#ifndef COMPILE_WITHOUT_DYADIC_SUPPORT
    else if(const DYADIC_IMPLICIT_OBJECT<TV>* implicit=dynamic_cast<const DYADIC_IMPLICIT_OBJECT<TV>*>(&implicit_input)) return Generate_Triangles_Helper(implicit); 
#endif
#ifndef COMPILE_WITHOUT_RLE_SUPPORT
    else if(const RLE_IMPLICIT_OBJECT<TV>* implicit=dynamic_cast<const RLE_IMPLICIT_OBJECT<TV>*>(&implicit_input)) return Generate_Triangles_Helper(implicit); 
#endif
    else{
#ifndef COMPILE_WITHOUT_READ_WRITE_SUPPORT
        LOG::cout<<"Trying to generate triangles on an object of type:\n\t"<<typeid(implicit_input).name()<<std::endl;
#endif
        PHYSBAM_FATAL_ERROR();}
}
//#####################################################################
#define INSTANTIATION_HELPER(T) \
    template POINT_SIMPLICES_1D<T>* Generate_Triangles(const IMPLICIT_OBJECT<VECTOR<T,1> >&); \
    template SEGMENTED_CURVE_2D<T>* Generate_Triangles(const IMPLICIT_OBJECT<VECTOR<T,2> >&); \
    template TRIANGULATED_SURFACE<T>* Generate_Triangles(const IMPLICIT_OBJECT<VECTOR<T,3> >&); 

INSTANTIATION_HELPER(float);
#ifndef COMPILE_WITHOUT_DOUBLE_SUPPORT
INSTANTIATION_HELPER(double);
#endif
}
}
