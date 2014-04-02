//#####################################################################
// Copyright 2009, Avi Robinson-Mosher.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class READ_WRITE_ANALYTIC_IMPLICIT_OBJECT
//#####################################################################
#ifndef COMPILE_WITHOUT_READ_WRITE_SUPPORT
#include <PhysBAM_Geometry/Read_Write/Implicit_Objects/READ_WRITE_ANALYTIC_IMPLICIT_OBJECT.h>
namespace PhysBAM{

void Register_Read_Write_Analytic_Implicit_Object()
{
#define DIMENSION_READ_WRITE_HELPER(T,RW) \
    Read_Write<STRUCTURE<VECTOR<T,1> >,RW>::Register_Read_Write<ANALYTIC_IMPLICIT_OBJECT<POINT_SIMPLEX_1D<T> > >(); \
    Read_Write<STRUCTURE<VECTOR<T,1> >,RW>::Register_Read_Write<ANALYTIC_IMPLICIT_OBJECT<BOX<VECTOR<T,1> > > >(); \
    Read_Write<STRUCTURE<VECTOR<T,2> >,RW>::Register_Read_Write<ANALYTIC_IMPLICIT_OBJECT<BOUNDED_HORIZONTAL_PLANE<VECTOR<T,2> > > >(); \
    Read_Write<STRUCTURE<VECTOR<T,2> >,RW>::Register_Read_Write<ANALYTIC_IMPLICIT_OBJECT<BOX<VECTOR<T,2> > > >(); \
    Read_Write<STRUCTURE<VECTOR<T,2> >,RW>::Register_Read_Write<ANALYTIC_IMPLICIT_OBJECT<LINE_2D<T> > >(); \
    Read_Write<STRUCTURE<VECTOR<T,2> >,RW>::Register_Read_Write<ANALYTIC_IMPLICIT_OBJECT<ORIENTED_BOX<VECTOR<T,2> > > >(); \
    Read_Write<STRUCTURE<VECTOR<T,2> >,RW>::Register_Read_Write<ANALYTIC_IMPLICIT_OBJECT<SPHERE<VECTOR<T,2> > > >(); \
    Read_Write<STRUCTURE<VECTOR<T,3> >,RW>::Register_Read_Write<ANALYTIC_IMPLICIT_OBJECT<BOX<VECTOR<T,3> > > >(); \
    Read_Write<STRUCTURE<VECTOR<T,3> >,RW>::Register_Read_Write<ANALYTIC_IMPLICIT_OBJECT<ORIENTED_BOX<VECTOR<T,3> > > >(); \
    Read_Write<STRUCTURE<VECTOR<T,3> >,RW>::Register_Read_Write<ANALYTIC_IMPLICIT_OBJECT<PLANE<T> > >(); \
    Read_Write<STRUCTURE<VECTOR<T,3> >,RW>::Register_Read_Write<ANALYTIC_IMPLICIT_OBJECT<BOUNDED_HORIZONTAL_PLANE<VECTOR<T,3> > > >(); \
    Read_Write<STRUCTURE<VECTOR<T,3> >,RW>::Register_Read_Write<ANALYTIC_IMPLICIT_OBJECT<CYLINDER<T> > >(); \
    Read_Write<STRUCTURE<VECTOR<T,3> >,RW>::Register_Read_Write<ANALYTIC_IMPLICIT_OBJECT<RING<T> > >(); \
    Read_Write<STRUCTURE<VECTOR<T,3> >,RW>::Register_Read_Write<ANALYTIC_IMPLICIT_OBJECT<SPHERE<VECTOR<T,3> > > >(); \
    Read_Write<STRUCTURE<VECTOR<T,3> >,RW>::Register_Read_Write<ANALYTIC_IMPLICIT_OBJECT<TORUS<T> > >();

#ifdef COMPILE_WITHOUT_DOUBLE_SUPPORT
#define RW_HELPER(RW) \
    DIMENSION_READ_WRITE_HELPER(float,RW)
#else
#define RW_HELPER(RW) \
    DIMENSION_READ_WRITE_HELPER(float,RW) \
    DIMENSION_READ_WRITE_HELPER(double,RW)
#endif

    RW_HELPER(float);
#ifndef COMPILE_WITHOUT_DOUBLE_SUPPORT
    RW_HELPER(double);
#endif
}
}
#endif
