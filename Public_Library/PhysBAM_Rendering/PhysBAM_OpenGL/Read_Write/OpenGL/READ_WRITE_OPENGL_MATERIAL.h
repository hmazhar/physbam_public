//#####################################################################
// Copyright 2009, Avi Robinson-Mosher.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class READ_WRITE_OPENGL_MATERIAL
//#####################################################################
#ifndef __READ_WRITE_OPENGL_MATERIAL__
#define __READ_WRITE_OPENGL_MATERIAL__

#include <PhysBAM_Tools/Read_Write/Utilities/READ_WRITE.h>
#include <PhysBAM_Rendering/PhysBAM_OpenGL/OpenGL/OPENGL_MATERIAL.h>
namespace PhysBAM{

template<class RW>
class Read_Write<OPENGL_MATERIAL,RW>
{
public:
    static void Read(std::istream& input,OPENGL_MATERIAL& object)
    {char version;Read_Binary<RW>(input,version);
    if(version!=1) throw READ_ERROR(STRING_UTILITIES::string_sprintf("Unrecognized OPENGL_MATERIAL version %d",version));
    Read_Binary<RW>(input,object.ambient,object.diffuse,object.specular,object.shininess);}

    static void Write(std::ostream& output,const OPENGL_MATERIAL& object)
    {char version=1;Write_Binary<RW>(output,version,object.ambient,object.diffuse,object.specular,object.shininess);}
};
}
#endif
