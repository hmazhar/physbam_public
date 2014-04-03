//#####################################################################
// Copyright 2009, Geoffrey Irving, Michael Lentine.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class REGISTER_GEOMETRY_READ_WRITE
//#####################################################################
#include <PhysBAM_Tools/Arrays/ARRAY_COLLECTION.h>
#include <PhysBAM_Tools/Matrices/FRAME.h>
#include <PhysBAM_Tools/Point_Clouds/POINT_CLOUD_FORWARD.h>
#include <PhysBAM_Tools/Read_Write/Arrays/READ_WRITE_ARRAY_COLLECTION.h>
#include <PhysBAM_Tools/Vectors/TWIST.h>
#include <PhysBAM_Tools/Vectors/VECTOR.h>
#include <PhysBAM_Geometry/Geometry_Particles/GEOMETRY_PARTICLES_FORWARD.h>
#include <PhysBAM_Geometry/Solids_Geometry/RIGID_GEOMETRY.h>
namespace PhysBAM{
bool Register_Implicit_Object_Transformed(); // TODO(jontg): Move these out of here
bool Register_Analytic_Implicit_Object();
void Register_Attribute_Name(const ATTRIBUTE_ID id,const char* name);
void Register_Free_Particles();

void Initialize_Geometry_Particle()
{
    static bool done=false;if(done) return;done=true;
    Register_Attribute_Name(ATTRIBUTE_ID_ROTATION,"rotation");
    Register_Attribute_Name(ATTRIBUTE_ID_TWIST,"twist");
    Register_Attribute_Name(ATTRIBUTE_ID_RIGID_GEOMETRY,"rigid_geometry");
    Register_Attribute_Name(ATTRIBUTE_ID_X,"X");
    Register_Attribute_Name(ATTRIBUTE_ID_V,"V");
    Register_Attribute_Name(ATTRIBUTE_ID_STRUCTURE_IDS,"structure_ids");
    Register_Attribute_Name(ATTRIBUTE_ID_ID,"id");
    Register_Attribute_Name(ATTRIBUTE_ID_COLOR,"color");
    Register_Attribute_Name(ATTRIBUTE_ID_RADIUS,"radius");
    Register_Attribute_Name(ATTRIBUTE_ID_DISPLAY_SIZE,"display_size");

    #define READ_WRITE_VECTOR_HELPER(T,RW,d) \
        Read_Write<ARRAY_COLLECTION,RW>::Register_Read_Write<VECTOR<T,d> >(); \
        Read_Write<ARRAY_COLLECTION,RW>::Register_Read_Write<ROTATION<VECTOR<T,d> > >(); \
        Read_Write<ARRAY_COLLECTION,RW>::Register_Read_Write<TWIST<VECTOR<T,d> > >(); \
        Read_Write<ARRAY_COLLECTION,RW>::Register_Read_Write<RIGID_GEOMETRY<VECTOR<T,d> >*>();

    #define READ_WRITE_SCALAR_HELPER(T,RW) \
        Read_Write<ARRAY_COLLECTION,RW>::Register_Read_Write<T>(); \
        Read_Write<ARRAY_COLLECTION,RW>::Register_Read_Write<VECTOR<T,0> >(); \
        READ_WRITE_VECTOR_HELPER(T,RW,1);READ_WRITE_VECTOR_HELPER(T,RW,2);READ_WRITE_VECTOR_HELPER(T,RW,3);

    #define READ_WRITE_HELPER(RW) \
        Read_Write<ARRAY_COLLECTION,RW>::Register_Read_Write<VECTOR<int,1> >(); \
        Read_Write<ARRAY_COLLECTION,RW>::Register_Read_Write<VECTOR<int,2> >(); \
        Read_Write<ARRAY_COLLECTION,RW>::Register_Read_Write<VECTOR<int,3> >(); \
        Read_Write<ARRAY_COLLECTION,RW>::Register_Read_Write<int>();

    READ_WRITE_SCALAR_HELPER(float,float);
    READ_WRITE_SCALAR_HELPER(float,double);
    READ_WRITE_HELPER(float);
    #ifndef COMPILE_WITHOUT_DOUBLE_SUPPORT
    READ_WRITE_SCALAR_HELPER(double,float);
    READ_WRITE_SCALAR_HELPER(double,double);
    READ_WRITE_HELPER(double);
    #endif

    Register_Implicit_Object_Transformed();
    Register_Analytic_Implicit_Object();
    Register_Free_Particles();
}
}
