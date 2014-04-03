//#####################################################################
// Copyright 2009, Avi Robinson-Mosher.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class READ_WRITE_OPENGL_RIGID_BODY_HINTS
//#####################################################################
#ifndef __READ_WRITE_OPENGL_RIGID_BODY_HINTS__
#define __READ_WRITE_OPENGL_RIGID_BODY_HINTS__

#include <PhysBAM_Tools/Read_Write/Utilities/READ_WRITE.h>
#include <PhysBAM_Rendering/PhysBAM_OpenGL/OpenGL/OPENGL_RIGID_BODY_HINTS.h>
namespace PhysBAM{

template<class RW>
class Read_Write<OPENGL_RIGID_BODY_HINTS,RW>
{
public:
    static void Read(std::istream& input,OPENGL_RIGID_BODY_HINTS& object)
    {char version;Read_Binary<RW>(input,version);
    if(version!=1) throw READ_ERROR(STRING_UTILITIES::string_sprintf("Unrecognized OPENGL_RIGID_BODY_HINTS version %d",version));
    Read_Binary<RW>(input,object.material,object.include_bounding_box);}

    static void Write(std::ostream& output,const OPENGL_RIGID_BODY_HINTS& object)
    {char version=1;Write_Binary<RW>(output,version,object.material,object.include_bounding_box);}
};
}
#endif
