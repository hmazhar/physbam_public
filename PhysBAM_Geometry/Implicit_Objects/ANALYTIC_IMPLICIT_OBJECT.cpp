//#####################################################################
// Copyright 2007, Geoffrey Irving, Andrew Selle.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class ANALYTIC_IMPLICIT_OBJECT
//#####################################################################
#include <PhysBAM_Geometry/Basic_Geometry/BOUNDED_HORIZONTAL_PLANE.h>
#include <PhysBAM_Geometry/Basic_Geometry/BOX.h>
#include <PhysBAM_Geometry/Basic_Geometry/CYLINDER.h>
#include <PhysBAM_Geometry/Basic_Geometry/LINE_2D.h>
#include <PhysBAM_Geometry/Basic_Geometry/ORIENTED_BOX.h>
#include <PhysBAM_Geometry/Basic_Geometry/PLANE.h>
#include <PhysBAM_Geometry/Basic_Geometry/POINT_SIMPLEX_1D.h>
#include <PhysBAM_Geometry/Basic_Geometry/RING.h>
#include <PhysBAM_Geometry/Basic_Geometry/SPHERE.h>
#include <PhysBAM_Geometry/Basic_Geometry/TORUS.h>
#include <PhysBAM_Geometry/Implicit_Objects/ANALYTIC_IMPLICIT_OBJECT.h>
#include <PhysBAM_Geometry/Registry/STRUCTURE_REGISTRY.h>
using namespace PhysBAM;
//#####################################################################
// Register this class as read-write
//#####################################################################
namespace PhysBAM{
bool Register_Analytic_Implicit_Object(){
    static bool done=false;if(done) return true;done=true;
    STRUCTURE_REGISTRY<VECTOR<float,1> >::Register<ANALYTIC_IMPLICIT_OBJECT<BOX<VECTOR<float,1> > > >();
    STRUCTURE_REGISTRY<VECTOR<float,1> >::Register<ANALYTIC_IMPLICIT_OBJECT<POINT_SIMPLEX_1D<float> > >();
    STRUCTURE_REGISTRY<VECTOR<float,2> >::Register<ANALYTIC_IMPLICIT_OBJECT<BOUNDED_HORIZONTAL_PLANE<VECTOR<float,2> > > >();
    STRUCTURE_REGISTRY<VECTOR<float,2> >::Register<ANALYTIC_IMPLICIT_OBJECT<BOX<VECTOR<float,2> > > >();
    STRUCTURE_REGISTRY<VECTOR<float,2> >::Register<ANALYTIC_IMPLICIT_OBJECT<LINE_2D<float> > >();
    STRUCTURE_REGISTRY<VECTOR<float,2> >::Register<ANALYTIC_IMPLICIT_OBJECT<ORIENTED_BOX<VECTOR<float,2> > > >();
    STRUCTURE_REGISTRY<VECTOR<float,2> >::Register<ANALYTIC_IMPLICIT_OBJECT<SPHERE<VECTOR<float,2> > > >();
    STRUCTURE_REGISTRY<VECTOR<float,3> >::Register<ANALYTIC_IMPLICIT_OBJECT<BOX<VECTOR<float,3> > > >();
    STRUCTURE_REGISTRY<VECTOR<float,3> >::Register<ANALYTIC_IMPLICIT_OBJECT<ORIENTED_BOX<VECTOR<float,3> > > >();
    STRUCTURE_REGISTRY<VECTOR<float,3> >::Register<ANALYTIC_IMPLICIT_OBJECT<PLANE<float> > >();
    STRUCTURE_REGISTRY<VECTOR<float,3> >::Register<ANALYTIC_IMPLICIT_OBJECT<BOUNDED_HORIZONTAL_PLANE<VECTOR<float,3> > > >();
    STRUCTURE_REGISTRY<VECTOR<float,3> >::Register<ANALYTIC_IMPLICIT_OBJECT<CYLINDER<float> > >();
    STRUCTURE_REGISTRY<VECTOR<float,3> >::Register<ANALYTIC_IMPLICIT_OBJECT<RING<float> > >();
    STRUCTURE_REGISTRY<VECTOR<float,3> >::Register<ANALYTIC_IMPLICIT_OBJECT<SPHERE<VECTOR<float,3> > > >();
    STRUCTURE_REGISTRY<VECTOR<float,3> >::Register<ANALYTIC_IMPLICIT_OBJECT<TORUS<float> > >();
#ifndef COMPILE_WITHOUT_DOUBLE_SUPPORT
    STRUCTURE_REGISTRY<VECTOR<double,1> >::Register<ANALYTIC_IMPLICIT_OBJECT<BOX<VECTOR<double,1> > > >();
    STRUCTURE_REGISTRY<VECTOR<double,1> >::Register<ANALYTIC_IMPLICIT_OBJECT<POINT_SIMPLEX_1D<double> > >();
    STRUCTURE_REGISTRY<VECTOR<double,2> >::Register<ANALYTIC_IMPLICIT_OBJECT<BOUNDED_HORIZONTAL_PLANE<VECTOR<double,2> > > >();
    STRUCTURE_REGISTRY<VECTOR<double,2> >::Register<ANALYTIC_IMPLICIT_OBJECT<BOX<VECTOR<double,2> > > >();
    STRUCTURE_REGISTRY<VECTOR<double,2> >::Register<ANALYTIC_IMPLICIT_OBJECT<LINE_2D<double> > >();
    STRUCTURE_REGISTRY<VECTOR<double,2> >::Register<ANALYTIC_IMPLICIT_OBJECT<ORIENTED_BOX<VECTOR<double,2> > > >();
    STRUCTURE_REGISTRY<VECTOR<double,2> >::Register<ANALYTIC_IMPLICIT_OBJECT<SPHERE<VECTOR<double,2> > > >();
    STRUCTURE_REGISTRY<VECTOR<double,3> >::Register<ANALYTIC_IMPLICIT_OBJECT<BOX<VECTOR<double,3> > > >();
    STRUCTURE_REGISTRY<VECTOR<double,3> >::Register<ANALYTIC_IMPLICIT_OBJECT<ORIENTED_BOX<VECTOR<double,3> > > >();
    STRUCTURE_REGISTRY<VECTOR<double,3> >::Register<ANALYTIC_IMPLICIT_OBJECT<PLANE<double> > >();
    STRUCTURE_REGISTRY<VECTOR<double,3> >::Register<ANALYTIC_IMPLICIT_OBJECT<BOUNDED_HORIZONTAL_PLANE<VECTOR<double,3> > > >();
    STRUCTURE_REGISTRY<VECTOR<double,3> >::Register<ANALYTIC_IMPLICIT_OBJECT<CYLINDER<double> > >();
    STRUCTURE_REGISTRY<VECTOR<double,3> >::Register<ANALYTIC_IMPLICIT_OBJECT<RING<double> > >();
    STRUCTURE_REGISTRY<VECTOR<double,3> >::Register<ANALYTIC_IMPLICIT_OBJECT<SPHERE<VECTOR<double,3> > > >();
    STRUCTURE_REGISTRY<VECTOR<double,3> >::Register<ANALYTIC_IMPLICIT_OBJECT<TORUS<double> > >();
#endif
    return true;
}
static bool registered=Register_Analytic_Implicit_Object();
}
