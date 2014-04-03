//#####################################################################
// Copyright 2007, Jonathan Su.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class OPENGL_COMPONENT_LEVELSET_1D
//#####################################################################
#ifndef __OPENGL_COMPONENT_LEVELSET_1D__
#define __OPENGL_COMPONENT_LEVELSET_1D__

#include <PhysBAM_Rendering/PhysBAM_OpenGL/OpenGL/OPENGL_LEVELSET_1D.h>
#include <PhysBAM_Rendering/PhysBAM_OpenGL/OpenGL_Components/OPENGL_COMPONENT.h>
namespace PhysBAM{

template<class T,class RW=T>
class OPENGL_COMPONENT_LEVELSET_1D:public OPENGL_COMPONENT
{
    typedef VECTOR<T,1> TV;typedef GRID<TV> T_GRID;
    typedef typename GRID_ARRAYS_POLICY<T_GRID>::ARRAYS_SCALAR T_ARRAYS_SCALAR;
private:
    std::string levelset_filename;
    int frame_loaded;
    bool valid;
public:
    OPENGL_LEVELSET_1D<T>* opengl_levelset;

//##################################################################### 
    OPENGL_COMPONENT_LEVELSET_1D(GRID<TV> &grid,const std::string& levelset_filename_input,OPENGL_COLOR point_color,OPENGL_COLOR line_color);
    ~OPENGL_COMPONENT_LEVELSET_1D();
    bool Valid_Frame(int frame_input) const PHYSBAM_OVERRIDE;
    void Set_Frame(int frame_input) PHYSBAM_OVERRIDE;
    void Set_Draw(bool draw_input=true) PHYSBAM_OVERRIDE;
    void Display(const int in_color=1) const PHYSBAM_OVERRIDE;
    virtual RANGE<VECTOR<float,3> > Bounding_Box() const PHYSBAM_OVERRIDE;
private:
    void Reinitialize();
//##################################################################### 
};
}
#endif
