//#####################################################################
// Copyright 2007, Craig Schroeder.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class EPS_FILE
//#####################################################################
#ifndef COMPILE_WITHOUT_READ_WRITE_SUPPORT
#ifndef __EPS_FILE__
#define __EPS_FILE__

#include <PhysBAM_Tools/Arrays/ARRAY.h>
#include <PhysBAM_Tools/Math_Tools/RANGE.h>
#include <PhysBAM_Tools/Read_Write/Utilities/FILE_UTILITIES.h>
#include <iostream>
#include <string>
namespace PhysBAM{

template<class T>
class EPS_FILE
{
protected:
    typedef VECTOR<T,2> TV;
    std::ostream* stream;
public:
    RANGE<TV> bounding_box;
    RANGE<TV> output_box;

    EPS_FILE(const std::string& filename,const RANGE<TV>& box=RANGE<TV>(TV(),TV(500,500)));
    ~EPS_FILE();

    void Finish();
    void Emit(const std::string& str);
    void Emit(const TV &pt);
    void Bound(const TV& pt);
    template<int d> void Bound(const VECTOR<TV,d>& pts);
    template<class T_OBJECT> void Draw_Object_Colored(const T_OBJECT& object,const VECTOR<T,3>& color);
    void Line_Color(const VECTOR<T,3>& color);
    void Write_Head();
    void Write_Tail();
    void Compute_Transform(T& scale,TV& shift);
    void Draw_Point(const TV &pt);
    void Draw_Line(const TV &a,const TV &b);
    void Draw_Object(const RANGE<TV>& box);
//#####################################################################
};
}
#endif
#endif
