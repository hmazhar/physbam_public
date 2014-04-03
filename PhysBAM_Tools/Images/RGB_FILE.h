//#####################################################################
// Copyright 2005, Ronald Fedkiw, Geoffrey Irving, Andrew Selle.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class RGB_FILE 
//#####################################################################
#ifndef COMPILE_WITHOUT_READ_WRITE_SUPPORT
#ifndef __RGB_FILE__
#define __RGB_FILE__

#include <PhysBAM_Tools/Arrays/ARRAYS_FORWARD.h>
#include <PhysBAM_Tools/Vectors/VECTOR_FORWARD.h>
#include <string>
namespace PhysBAM{

template<class T>
class RGB_FILE
{
public:
    RGB_FILE()
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
