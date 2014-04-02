//#####################################################################
// Copyright 2007, Craig Schroeder, Tamar Shinar.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class OPENGL_FREE_PARTICLES
//#####################################################################
#ifndef __OPENGL_FREE_PARTICLES__
#define __OPENGL_FREE_PARTICLES__

#include <PhysBAM_Geometry/Solids_Geometry/DEFORMABLE_GEOMETRY_COLLECTION.h>
#include <PhysBAM_Rendering/PhysBAM_OpenGL/OpenGL/OPENGL_POINTS_2D.h>
#include <PhysBAM_Rendering/PhysBAM_OpenGL/OpenGL/OPENGL_POINTS_3D.h>
namespace PhysBAM{

template<class TV> struct OPENGL_POLICY;
template<class T> struct OPENGL_POLICY<VECTOR<T,2> >{typedef OPENGL_POINTS_2D<T,INDIRECT_ARRAY<ARRAY_VIEW<VECTOR<T,2> > > > OPENGL_POINTS;};
template<class T> struct OPENGL_POLICY<VECTOR<T,3> >{typedef OPENGL_POINTS_3D<T,INDIRECT_ARRAY<ARRAY_VIEW<VECTOR<T,3> > > > OPENGL_POINTS;};

template<class TV_input>
class OPENGL_FREE_PARTICLES:public OPENGL_POLICY<TV_input>::OPENGL_POINTS
{
    typedef TV_input TV;typedef typename TV::SCALAR T;typedef typename OPENGL_POLICY<TV>::OPENGL_POINTS BASE;
    using BASE::points;
public:
    DEFORMABLE_GEOMETRY_COLLECTION<TV>& deformable_geometry_collection;

    OPENGL_FREE_PARTICLES(DEFORMABLE_GEOMETRY_COLLECTION<TV>& deformable_geometry_collection,INDIRECT_ARRAY<ARRAY_VIEW<TV> >& points,const OPENGL_COLOR& color=OPENGL_COLOR::White(),const T point_size=5);
    ~OPENGL_FREE_PARTICLES(){}

//#####################################################################
    void Print_Selection_Info(std::ostream &output_stream,OPENGL_SELECTION *selection) const PHYSBAM_OVERRIDE;
    int Particle_Index(const int index) const PHYSBAM_OVERRIDE {return points.indices(index);}
//#####################################################################
};
}
#endif
