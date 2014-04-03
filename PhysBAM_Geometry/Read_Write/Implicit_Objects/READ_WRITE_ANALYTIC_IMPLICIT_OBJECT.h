//#####################################################################
// Copyright 2009, Avi Robinson-Mosher.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class READ_WRITE_ANALYTIC_IMPLICIT_OBJECT
//#####################################################################
#ifndef COMPILE_WITHOUT_READ_WRITE_SUPPORT
#ifndef __READ_WRITE_ANALYTIC_IMPLICIT_OBJECT__
#define __READ_WRITE_ANALYTIC_IMPLICIT_OBJECT__

#include <PhysBAM_Tools/Read_Write/Vectors/READ_WRITE_VECTOR.h>
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
#include <PhysBAM_Geometry/Read_Write/Geometry/READ_WRITE_BOUNDED_HORIZONTAL_PLANE.h>
#include <PhysBAM_Geometry/Read_Write/Geometry/READ_WRITE_BOX.h>
#include <PhysBAM_Geometry/Read_Write/Geometry/READ_WRITE_CYLINDER.h>
#include <PhysBAM_Geometry/Read_Write/Geometry/READ_WRITE_LINE_2D.h>
#include <PhysBAM_Geometry/Read_Write/Geometry/READ_WRITE_ORIENTED_BOX.h>
#include <PhysBAM_Geometry/Read_Write/Geometry/READ_WRITE_PLANE.h>
#include <PhysBAM_Geometry/Read_Write/Geometry/READ_WRITE_POINT_SIMPLEX_1D.h>
#include <PhysBAM_Geometry/Read_Write/Geometry/READ_WRITE_RING.h>
#include <PhysBAM_Geometry/Read_Write/Geometry/READ_WRITE_SPHERE.h>
#include <PhysBAM_Geometry/Read_Write/Geometry/READ_WRITE_TORUS.h>
#include <PhysBAM_Geometry/Read_Write/Implicit_Objects/READ_WRITE_IMPLICIT_OBJECT.h>
namespace PhysBAM{

void Register_Read_Write_Analytic_Implicit_Object();

template<class RW,class T_ANALYTIC>
class Read_Write<ANALYTIC_IMPLICIT_OBJECT<T_ANALYTIC>,RW>:public Read_Write<IMPLICIT_OBJECT<typename T_ANALYTIC::VECTOR_T>,RW>
{
    typedef typename T_ANALYTIC::VECTOR_T TV;
public:
    static void Read_Helper(std::istream& input,STRUCTURE<TV>& structure_object)
    {ANALYTIC_IMPLICIT_OBJECT<T_ANALYTIC>& object=dynamic_cast<ANALYTIC_IMPLICIT_OBJECT<T_ANALYTIC>&>(structure_object);
    Read_Binary<RW>(input,object.analytic);object.Update_Box();}

    static void Read_Structure_Helper(std::istream& input,STRUCTURE<TV>& structure_object)
    {ANALYTIC_IMPLICIT_OBJECT<T_ANALYTIC>& object=dynamic_cast<ANALYTIC_IMPLICIT_OBJECT<T_ANALYTIC>&>(structure_object);
    Read_Binary<RW>(input,object.analytic);object.Update_Box();}

    static void Write_Helper(std::ostream& output,const STRUCTURE<TV>& structure_object)
    {const ANALYTIC_IMPLICIT_OBJECT<T_ANALYTIC>& object=dynamic_cast<const ANALYTIC_IMPLICIT_OBJECT<T_ANALYTIC>&>(structure_object);
    Write_Binary<RW>(output,object.analytic);}

    static void Write_Structure_Helper(std::ostream& output,const STRUCTURE<TV>& structure_object)
    {const ANALYTIC_IMPLICIT_OBJECT<T_ANALYTIC>& object=dynamic_cast<const ANALYTIC_IMPLICIT_OBJECT<T_ANALYTIC>&>(structure_object);
    Write_Binary<RW>(output,object.analytic);}
};
}
#endif
#endif
