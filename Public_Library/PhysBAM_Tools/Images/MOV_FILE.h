//#####################################################################
// Copyright 2008, Andrew Selle.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#ifndef COMPILE_WITHOUT_READ_WRITE_SUPPORT
#ifndef __MOV_FILE__
#define __MOV_FILE__

#include <PhysBAM_Tools/Arrays/ARRAY.h>
#include <cstdio>
namespace PhysBAM{
class QT_ATOM;

template<class T>
class MOV_WRITER
{
    int frames_per_second;
    int width,height;
    FILE* fp;
    QT_ATOM* current_mov;
    ARRAY<int> sample_offsets;
    ARRAY<int> sample_lengths;

public:
    MOV_WRITER(const std::string& filename,const int frames_per_second=24);
    ~MOV_WRITER();
    void Add_Frame(ARRAY<VECTOR<T,3> ,VECTOR<int,2> >& image);
    void Write_Footer();
    static bool Enabled();
//#####################################################################
};
}
#endif
#endif
