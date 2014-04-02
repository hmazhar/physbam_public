//#####################################################################
// Copyright 2009, Avi Robinson-Mosher.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class READ_WRITE_STRUCTURE
//#####################################################################
#ifndef COMPILE_WITHOUT_READ_WRITE_SUPPORT
#include <PhysBAM_Tools/Vectors/VECTOR.h>
#include <PhysBAM_Geometry/Read_Write/Geometry/READ_WRITE_STRUCTURE.h>
namespace PhysBAM{
//#####################################################################
// Function Read_Registry
//#####################################################################
template<class RW,class TV> HASHTABLE<std::string,void (*)(std::istream&,STRUCTURE<TV>&)>& Read_Write<STRUCTURE<TV>,RW>::
Read_Registry()
{
    static HASHTABLE<std::string,void (*)(std::istream&,STRUCTURE<TV>&)> read_registry;
    return read_registry;
}
//#####################################################################
// Function Read_Structure_Registry
//#####################################################################
template<class RW,class TV> HASHTABLE<std::string,void (*)(std::istream&,STRUCTURE<TV>&)>& Read_Write<STRUCTURE<TV>,RW>::
Read_Structure_Registry()
{
    static HASHTABLE<std::string,void (*)(std::istream&,STRUCTURE<TV>&)> read_structure_registry;
    return read_structure_registry;
}
//#####################################################################
// Function Write_Registry
//#####################################################################
template<class RW,class TV> HASHTABLE<std::string,void (*)(std::ostream&,const STRUCTURE<TV>&)>& Read_Write<STRUCTURE<TV>,RW>::
Write_Registry()
{
    static HASHTABLE<std::string,void (*)(std::ostream&,const STRUCTURE<TV>&)> write_registry;
    return write_registry;
}
//#####################################################################
// Function Write_Structure_Registry
//#####################################################################
template<class RW,class TV> HASHTABLE<std::string,void (*)(std::ostream&,const STRUCTURE<TV>&)>& Read_Write<STRUCTURE<TV>,RW>::
Write_Structure_Registry()
{
    static HASHTABLE<std::string,void (*)(std::ostream&,const STRUCTURE<TV>&)> write_structure_registry;
    return write_structure_registry;
}
//#####################################################################

void Register_Read_Write_Structure();
void Register_Read_Write_Embedded_Object();
void Register_Read_Write_Embedded_Tetrahedralized_Volume();
void Register_Read_Write_Embedded_Triangulated_Object();
void Register_Read_Write_Mesh_Object();
void Register_Read_Write_Point_Simplices_1d();
void Register_Read_Write_Segmented_Curve_2d();
void Register_Read_Write_Segmented_Curve();
void Register_Read_Write_Hexahedralized_Volume();
void Register_Read_Write_Tetrahedralized_Volume();
void Register_Read_Write_Triangulated_Area();
void Register_Read_Write_Triangulated_Surface();
void Register_Read_Write_Analytic_Implicit_Object();
#ifndef COMPILE_WITHOUT_DYADIC_SUPPORT
void Register_Read_Write_Dyadic_Implicit_Object();
#endif
void Register_Read_Write_Implicit_Object();
void Register_Read_Write_Implicit_Object_Transformed();
void Register_Read_Write_Levelset_Implicit_Object();
void Register_Read_Write_Multibody_Levelset_Implicit_Object();
#ifndef COMPILE_WITHOUT_RLE_SUPPORT
void Register_Read_Write_RLE_Implicit_Object();
#endif
void Register_Read_Write_Free_Particles();

void Initialize_Read_Write_Structures()
{
    static bool done=false;if(done) return;done=true;
    Register_Read_Write_Embedded_Object();
    Register_Read_Write_Embedded_Tetrahedralized_Volume();
    Register_Read_Write_Embedded_Triangulated_Object();
    Register_Read_Write_Mesh_Object();
    Register_Read_Write_Point_Simplices_1d();
    Register_Read_Write_Segmented_Curve_2d();
    Register_Read_Write_Segmented_Curve();
    //Register_Read_Write_Structure();
    Register_Read_Write_Hexahedralized_Volume();
    Register_Read_Write_Tetrahedralized_Volume();
    Register_Read_Write_Triangulated_Area();
    Register_Read_Write_Triangulated_Surface();
    Register_Read_Write_Analytic_Implicit_Object();
#ifndef COMPILE_WITHOUT_DYADIC_SUPPORT
    Register_Read_Write_Dyadic_Implicit_Object();
#endif
    //Register_Read_Write_Implicit_Object();
    Register_Read_Write_Implicit_Object_Transformed();
    Register_Read_Write_Levelset_Implicit_Object();
    Register_Read_Write_Multibody_Levelset_Implicit_Object();
#ifndef COMPILE_WITHOUT_RLE_SUPPORT
    Register_Read_Write_RLE_Implicit_Object();
#endif
    Register_Read_Write_Free_Particles();
}

#define INSTANTIATION_HELPER_DIMENSION(T,RW) \
    template class Read_Write<STRUCTURE<VECTOR<T,1> >,RW>; \
    template class Read_Write<STRUCTURE<VECTOR<T,2> >,RW>; \
    template class Read_Write<STRUCTURE<VECTOR<T,3> >,RW>;

INSTANTIATION_HELPER_DIMENSION(float,float)
INSTANTIATION_HELPER_DIMENSION(float,double)
#ifndef COMPILE_WITHOUT_DOUBLE_SUPPORT
INSTANTIATION_HELPER_DIMENSION(double,float)
INSTANTIATION_HELPER_DIMENSION(double,double)
#endif

}
#endif
