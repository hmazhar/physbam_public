//#####################################################################
// Copyright 2009.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class OPENGL_EPS_OUTPUT
//#####################################################################
#ifndef __OPENGL_EPS_OUTPUT__
#define __OPENGL_EPS_OUTPUT__
#include <PhysBAM_Tools/Arrays/ARRAY.h>
#include <PhysBAM_Tools/Vectors/VECTOR.h>
namespace PhysBAM{
template<class T>
class OPENGL_EPS_OUTPUT
{
protected:
    typedef VECTOR<T,3> TV;
    std::ostream& stream;
public:
    int glmode;
    ARRAY<TV> buffer;

    OPENGL_EPS_OUTPUT(const std::string& filename);
    virtual ~OPENGL_EPS_OUTPUT();

    void Begin(int mode);
    void End();
    template<class T2> void Vertex(const VECTOR<T2,3>& p);
    template<class T2> void Vertex(const VECTOR<T2,2>& p);
    template<class T2> void Vertex(const VECTOR<T2,1>& p);
    void Head();
    void Tail();
    void Emit(const TV& p);
    void Emit(const char* p);
    void Set_Color(const TV& color);
protected:
    void Draw_Point(const TV& p);
    void Draw_Line(const TV& a,const TV& b);
    void Draw_Polygon(int i,int n);
    void Transform_Buffer();
};
}
#endif
