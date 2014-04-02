//#####################################################################
// Copyright 2009, Avi Robinson-Mosher.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class READ_WRITE_OPENGL_COLOR
//#####################################################################
#ifndef __READ_WRITE_OPENGL_COLOR__
#define __READ_WRITE_OPENGL_COLOR__

#include <PhysBAM_Tools/Read_Write/Utilities/READ_WRITE.h>
#include <PhysBAM_Rendering/PhysBAM_OpenGL/OpenGL/OPENGL_COLOR.h>
namespace PhysBAM{

template<class RW>
class Read_Write<OPENGL_COLOR,RW>
{
public:
    static void Read(std::istream& input,OPENGL_COLOR& object)
    {for(int i=0;i<4;i++){float tmp;Read_Binary<float>(input,tmp);object.rgba[i]=tmp;}} // Read as float

    static void Write(std::ostream& output,const OPENGL_COLOR& object)
    {for(int i=0;i<4;i++){Write_Binary<float>(output,(float)object.rgba[i]);}}
};
}
#endif
