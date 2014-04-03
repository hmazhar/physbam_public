#ifndef COMPILE_WITHOUT_DYADIC_SUPPORT
//#####################################################################
// Copyright 2004-2005, Frank Losasso, Andrew Selle.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class OPENGL_COMPONENT_OCTREE_NODE_VECTOR_FIELD
//#####################################################################
#ifndef __OPENGL_COMPONENT_OCTREE_NODE_VECTOR_FIELD__
#define __OPENGL_COMPONENT_OCTREE_NODE_VECTOR_FIELD__

#include <PhysBAM_Rendering/PhysBAM_OpenGL/OpenGL/OPENGL_OCTREE_NODE_VECTOR_FIELD.h>
#include <PhysBAM_Rendering/PhysBAM_OpenGL/OpenGL_Components/OPENGL_COMPONENT.h>

namespace PhysBAM
{
template<class T,class RW=T>
class OPENGL_COMPONENT_OCTREE_NODE_VECTOR_FIELD:public OPENGL_COMPONENT
{
public:
    OPENGL_OCTREE_NODE_VECTOR_FIELD<T> opengl_octree_node_vector_field;

private:
    std::string vector_field_filename;
    int frame_loaded;
    bool valid;
public:
    OPENGL_COMPONENT_OCTREE_NODE_VECTOR_FIELD(OCTREE_GRID<T> &grid,const std::string &vector_field_filename_input);
    ~OPENGL_COMPONENT_OCTREE_NODE_VECTOR_FIELD();

    bool Valid_Frame(int frame_input) const PHYSBAM_OVERRIDE;
    bool Is_Up_To_Date(int frame) const PHYSBAM_OVERRIDE { return valid && frame_loaded == frame; }

    void Set_Frame(int frame_input) PHYSBAM_OVERRIDE;
    void Set_Draw(bool draw_input = true) PHYSBAM_OVERRIDE;

    void Display(const int in_color=1) const PHYSBAM_OVERRIDE;
    bool Use_Bounding_Box() const PHYSBAM_OVERRIDE { return draw && valid; }
    virtual RANGE<VECTOR<float,3> > Bounding_Box() const PHYSBAM_OVERRIDE;

    virtual void Set_Slice(OPENGL_SLICE *slice_input) PHYSBAM_OVERRIDE {slice=slice_input;opengl_octree_node_vector_field.Set_Slice(slice_input);}
    virtual void Slice_Has_Changed() {if(draw) opengl_octree_node_vector_field.Slice_Has_Changed();}

    void Print_Selection_Info(std::ostream& output_stream,OPENGL_SELECTION* current_selection) const PHYSBAM_OVERRIDE;

    void Increase_Vector_Size();
    void Decrease_Vector_Size();
    void Toggle_Arrowhead();

    DEFINE_COMPONENT_CALLBACK(OPENGL_COMPONENT_OCTREE_NODE_VECTOR_FIELD, Increase_Vector_Size, "Increase vector size");
    DEFINE_COMPONENT_CALLBACK(OPENGL_COMPONENT_OCTREE_NODE_VECTOR_FIELD, Decrease_Vector_Size, "Decrease vector size");
    DEFINE_COMPONENT_CALLBACK(OPENGL_COMPONENT_OCTREE_NODE_VECTOR_FIELD, Toggle_Arrowhead, "Toggle arrowhead");

private:
    void Reinitialize();

};

}

#endif
#endif
