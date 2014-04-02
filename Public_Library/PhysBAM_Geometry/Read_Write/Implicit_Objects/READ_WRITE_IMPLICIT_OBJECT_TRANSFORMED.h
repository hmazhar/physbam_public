//#####################################################################
// Copyright 2009, Avi Robinson-Mosher.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class READ_WRITE_IMPLICIT_OBJECT_TRANSFORMED
//#####################################################################
#ifndef COMPILE_WITHOUT_READ_WRITE_SUPPORT
#ifndef __READ_WRITE_IMPLICIT_OBJECT_TRANSFORMED__
#define __READ_WRITE_IMPLICIT_OBJECT_TRANSFORMED__

#include <PhysBAM_Tools/Read_Write/Matrices_And_Vectors/READ_WRITE_FRAME.h>
#include <PhysBAM_Geometry/Implicit_Objects/IMPLICIT_OBJECT_TRANSFORMED.h>
#include <PhysBAM_Geometry/Read_Write/Implicit_Objects/READ_WRITE_IMPLICIT_OBJECT.h>
namespace PhysBAM{

template<class RW,class TV,class TRANSFORM>
class Read_Write<IMPLICIT_OBJECT_TRANSFORMED_HELPER<TV,TRANSFORM>,RW>
{
public:
    template<class READ_TRANSFORM> static void Read_Transform_Helper(std::istream& input,IMPLICIT_OBJECT_TRANSFORMED_HELPER<TV,TRANSFORM>& object,const READ_TRANSFORM& transform_input)
    {PHYSBAM_FATAL_ERROR();}

    static void Read_Transform_Helper(std::istream& input,IMPLICIT_OBJECT_TRANSFORMED_HELPER<TV,TRANSFORM>& object,const FRAME<TV>& transform_input)
    {object.transform=new FRAME<TV>();Read_Binary<RW>(input,(FRAME<TV>&)*object.transform);object.owns_transform=true;}

    static void Read_Transform_Helper(std::istream& input,IMPLICIT_OBJECT_TRANSFORMED_HELPER<TV,TRANSFORM>& object)
    {Read_Transform_Helper(input,object,*object.transform);}

    static void Write_Transform_Helper(std::ostream& output,const IMPLICIT_OBJECT_TRANSFORMED_HELPER<TV,TRANSFORM>& object)
    {Write_Binary<RW>(output,*object.transform);}
};

template<class RW,class TV>
class Read_Write<IMPLICIT_OBJECT_TRANSFORMED_HELPER<TV,typename TV::SCALAR>,RW>
{
    typedef typename TV::SCALAR TRANSFORM;
public:
    static void Read_Transform_Helper(std::istream& input,IMPLICIT_OBJECT_TRANSFORMED_HELPER<TV,TRANSFORM>& object)
    {Read_Binary<RW>(input,object.center,object.scale,object.one_over_scale,object.center_adjustment);}

    static void Write_Transform_Helper(std::ostream& output,const IMPLICIT_OBJECT_TRANSFORMED_HELPER<TV,TRANSFORM>& object)
    {Write_Binary<RW>(output,object.center,object.scale,object.one_over_scale,object.center_adjustment);}
};

template<class RW,class TV,class TRANSFORM>
class Read_Write<IMPLICIT_OBJECT_TRANSFORMED<TV,TRANSFORM>,RW>:public Read_Write<IMPLICIT_OBJECT<TV>,RW>
{
public:
    static void Read_Helper(std::istream& input,STRUCTURE<TV>& structure_object)
    {IMPLICIT_OBJECT_TRANSFORMED<TV,TRANSFORM>& object=dynamic_cast<IMPLICIT_OBJECT_TRANSFORMED<TV,TRANSFORM>&>(structure_object);
    Read_Write<typename IMPLICIT_OBJECT_TRANSFORMED<TV,TRANSFORM>::BASE_HELPER,RW>::Read_Transform_Helper(input,object);std::string name;Read_Binary<RW>(input,name);
    object.object_space_implicit_object=dynamic_cast<IMPLICIT_OBJECT<TV>*>(STRUCTURE<TV>::Create_From_Name(name));
    Read_Binary<RW>(input,*object.object_space_implicit_object);}

    static void Write_Helper(std::ostream& output,const STRUCTURE<TV>& structure_object)
    {const IMPLICIT_OBJECT_TRANSFORMED<TV,TRANSFORM>& object=dynamic_cast<const IMPLICIT_OBJECT_TRANSFORMED<TV,TRANSFORM>&>(structure_object);
    Read_Write<typename IMPLICIT_OBJECT_TRANSFORMED<TV,TRANSFORM>::BASE_HELPER,RW>::Write_Transform_Helper(output,object);
    Write_Binary<RW>(output,object.object_space_implicit_object->Name());Write_Binary<RW>(output,*object.object_space_implicit_object);}
};
}
#endif
#endif
