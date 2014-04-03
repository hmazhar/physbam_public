//#####################################################################
// Copyright 2002-2005, Geoffrey Irving, Igor Neverov, Andrew Selle.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class EXR_FILE 
// Requires OpenEXR (http://www.openexr.com/)
//#####################################################################
#ifndef COMPILE_WITHOUT_READ_WRITE_SUPPORT
#ifndef _EXR_FILE_h
#define _EXR_FILE_h

#include <PhysBAM_Tools/Arrays/ARRAYS_FORWARD.h>
#include <PhysBAM_Tools/Log/DEBUG_UTILITIES.h>
#include <PhysBAM_Tools/Vectors/VECTOR_FORWARD.h>
#include <string>
namespace PhysBAM{

template<class T>
class EXR_FILE
{
public:
    EXR_FILE()
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
