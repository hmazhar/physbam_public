#ifndef COMPILE_WITHOUT_DYADIC_SUPPORT
//#####################################################################
// Copyright 2004, Frank Losasso.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#include <PhysBAM_Tools/Read_Write/Arrays/READ_WRITE_ARRAY.h>
#include <PhysBAM_Tools/Read_Write/Utilities/FILE_UTILITIES.h>
#include <PhysBAM_Rendering/PhysBAM_OpenGL/OpenGL_Components/OPENGL_COMPONENT_QUADTREE_NODE_VECTOR_FIELD.h>
using namespace PhysBAM;

template<class T,class RW> OPENGL_COMPONENT_QUADTREE_NODE_VECTOR_FIELD<T,RW>::
OPENGL_COMPONENT_QUADTREE_NODE_VECTOR_FIELD(QUADTREE_GRID<T> &grid,const std::string &vector_field_filename)
: OPENGL_COMPONENT("Quadtree Based Vector Field"), opengl_quadtree_node_vector_field(grid,*(new ARRAY<VECTOR<T,2> >)), 
                    vector_field_filename(vector_field_filename), valid(false)
{
    is_animation = FILE_UTILITIES::Is_Animated(vector_field_filename);
    frame_loaded = -1;

    Reinitialize();
}

template<class T,class RW> OPENGL_COMPONENT_QUADTREE_NODE_VECTOR_FIELD<T,RW>::
~OPENGL_COMPONENT_QUADTREE_NODE_VECTOR_FIELD()
{
}

template<class T,class RW> bool OPENGL_COMPONENT_QUADTREE_NODE_VECTOR_FIELD<T,RW>::
Valid_Frame(int frame_input) const
{
    return FILE_UTILITIES::Frame_File_Exists(vector_field_filename, frame_input);
}

template<class T,class RW> void OPENGL_COMPONENT_QUADTREE_NODE_VECTOR_FIELD<T,RW>::
Set_Frame(int frame_input)
{
    OPENGL_COMPONENT::Set_Frame(frame_input);
    Reinitialize();
}

template<class T,class RW> void OPENGL_COMPONENT_QUADTREE_NODE_VECTOR_FIELD<T,RW>::
Set_Draw(bool draw_input)
{
    OPENGL_COMPONENT::Set_Draw(draw_input);
    Reinitialize();
}

template<class T,class RW> void OPENGL_COMPONENT_QUADTREE_NODE_VECTOR_FIELD<T,RW>::
Display(const int in_color) const
{
    if (valid && draw) opengl_quadtree_node_vector_field.Display(in_color);
}

template<class T,class RW> RANGE<VECTOR<float,3> > OPENGL_COMPONENT_QUADTREE_NODE_VECTOR_FIELD<T,RW>::
Bounding_Box() const
{
    if (valid && draw) return opengl_quadtree_node_vector_field.Bounding_Box();
    else return RANGE<VECTOR<float,3> >::Centered_Box();
}

template<class T,class RW> void OPENGL_COMPONENT_QUADTREE_NODE_VECTOR_FIELD<T,RW>::
Reinitialize()
{
    if (draw)
    {
        if (!valid ||
            (is_animation && frame_loaded != frame) ||
            (!is_animation && frame_loaded < 0))
        {
            valid = false;

            std::string tmp_filename = FILE_UTILITIES::Get_Frame_Filename(vector_field_filename, frame);
            if (FILE_UTILITIES::File_Exists(tmp_filename))
                FILE_UTILITIES::Read_From_File<RW>(tmp_filename,opengl_quadtree_node_vector_field.V);
            else
                return;

            opengl_quadtree_node_vector_field.Update();
            frame_loaded = frame;
            valid = true;
        }
    }
}

template<class T,class RW> void OPENGL_COMPONENT_QUADTREE_NODE_VECTOR_FIELD<T,RW>::
Print_Selection_Info(std::ostream& output_stream,OPENGL_SELECTION* current_selection) const
{
    if(Is_Up_To_Date(frame)){
        output_stream<<component_name<<": ";
        opengl_quadtree_node_vector_field.Print_Selection_Info(output_stream,current_selection);}
}

template<class T,class RW> void OPENGL_COMPONENT_QUADTREE_NODE_VECTOR_FIELD<T,RW>::
Increase_Vector_Size()
{
    opengl_quadtree_node_vector_field.Scale_Vector_Size(1.1);
}

template<class T,class RW> void OPENGL_COMPONENT_QUADTREE_NODE_VECTOR_FIELD<T,RW>::
Decrease_Vector_Size()
{
    opengl_quadtree_node_vector_field.Scale_Vector_Size(1/1.1);
}

template<class T,class RW> void OPENGL_COMPONENT_QUADTREE_NODE_VECTOR_FIELD<T,RW>::
Toggle_Arrowhead()
{
    opengl_quadtree_node_vector_field.draw_arrowhead = !opengl_quadtree_node_vector_field.draw_arrowhead;
}
template class OPENGL_COMPONENT_QUADTREE_NODE_VECTOR_FIELD<float,float>;
#ifndef COMPILE_WITHOUT_DOUBLE_SUPPORT
template class OPENGL_COMPONENT_QUADTREE_NODE_VECTOR_FIELD<double,double>;
#endif
#endif
