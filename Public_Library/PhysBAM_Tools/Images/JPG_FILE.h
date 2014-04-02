//#####################################################################
// Copyright 2003, Eran Guendelman.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class JPG_FILE 
//
// Requires libjpeg (http://www.ijg.org/)
// NOTE: Appears to give the error message
//   JPEG parameter struct mismatch: library thinks size is 372, caller expects 376
// with gcc which goes away if you don't use -malign-double
//#####################################################################
#ifndef COMPILE_WITHOUT_READ_WRITE_SUPPORT
#ifndef _JPG_FILE_h
#define _JPG_FILE_h

#include <PhysBAM_Tools/Arrays/ARRAYS_FORWARD.h>
#include <PhysBAM_Tools/Vectors/VECTOR_FORWARD.h>
#include <string>
namespace PhysBAM{

template<class T>
class JPG_FILE
{
public:
    JPG_FILE()
    {}
    
//#####################################################################
    static void Read(const std::string& filename,ARRAY<VECTOR<T,3> ,VECTOR<int,2> >& image);
    static void Read(const std::string& filename,ARRAY<VECTOR<T,4> ,VECTOR<int,2> >& image);
    template<int d> static void Write(const std::string& filename,const ARRAY<VECTOR<T,d> ,VECTOR<int,2> >& image);
    static bool Is_Supported();
//#####################################################################
};
}
#endif
#endif
