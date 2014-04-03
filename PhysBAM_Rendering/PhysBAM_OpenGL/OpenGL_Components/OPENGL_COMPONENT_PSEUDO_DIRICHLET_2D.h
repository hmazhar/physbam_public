//#####################################################################
// Copyright 2004, Eran Guendelman.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class OPENGL_COMPONENT_PSEUDO_DIRICHLET_2D
//#####################################################################
#ifndef __OPENGL_COMPONENT_PSEUDO_DIRICHLET_2D__
#define __OPENGL_COMPONENT_PSEUDO_DIRICHLET_2D__

#include <PhysBAM_Tools/Arrays/ARRAY.h>
#include <PhysBAM_Tools/Data_Structures/TRIPLE.h>
#include <PhysBAM_Tools/Grids_Uniform/GRID.h>
#include <PhysBAM_Tools/Vectors/VECTOR_2D.h>
#include <PhysBAM_Rendering/PhysBAM_OpenGL/OpenGL_Components/OPENGL_COMPONENT.h>

namespace PhysBAM
{
template<class TV> class GRID;

template<class T,class RW=T>
class OPENGL_COMPONENT_PSEUDO_DIRICHLET_2D : public OPENGL_COMPONENT
{
    typedef VECTOR<T,2> TV;
public:
    GRID<TV> mac_grid;
    ARRAY<TRIPLE<VECTOR<int,2>,VECTOR<T,2>,char> > pseudo_dirichlet_cells;
private:
    T velocity_scale;
    std::string filename;
    int frame_loaded;
    bool valid;
public:
    OPENGL_COMPONENT_PSEUDO_DIRICHLET_2D(const GRID<TV> &grid,const std::string &filename_input);
    
    bool Valid_Frame(int frame_input) const PHYSBAM_OVERRIDE;
    bool Is_Up_To_Date(int frame) const PHYSBAM_OVERRIDE { return valid && frame_loaded == frame; }
    void Set_Frame(int frame_input) PHYSBAM_OVERRIDE;
    void Set_Draw(bool draw_input = true) PHYSBAM_OVERRIDE;
    void Display(const int in_color=1) const PHYSBAM_OVERRIDE;
    bool Use_Bounding_Box() const PHYSBAM_OVERRIDE { return draw && valid; }
    virtual RANGE<VECTOR<float,3> > Bounding_Box() const PHYSBAM_OVERRIDE;

    void Set_Vector_Size(const T vector_size);

    void Increase_Vector_Size();
    void Decrease_Vector_Size();

    DEFINE_COMPONENT_CALLBACK(OPENGL_COMPONENT_PSEUDO_DIRICHLET_2D, Increase_Vector_Size, "Increase vector size");
    DEFINE_COMPONENT_CALLBACK(OPENGL_COMPONENT_PSEUDO_DIRICHLET_2D, Decrease_Vector_Size, "Decrease vector size");

private:
    void Reinitialize(bool force=false);
};
}
#endif
