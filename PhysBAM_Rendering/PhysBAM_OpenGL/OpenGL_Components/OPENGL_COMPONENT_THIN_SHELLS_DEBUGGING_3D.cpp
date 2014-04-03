//#####################################################################
// Copyright 2004, Eran Guendelman.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#include <PhysBAM_Tools/Read_Write/Data_Structures/READ_WRITE_PAIR.h>
#include <PhysBAM_Tools/Read_Write/Grids_Uniform_Arrays/READ_WRITE_ARRAYS.h>
#include <PhysBAM_Tools/Read_Write/Math_Tools/READ_WRITE_RANGE.h>
#include <PhysBAM_Tools/Read_Write/Utilities/FILE_UTILITIES.h>
#include <PhysBAM_Rendering/PhysBAM_OpenGL/OpenGL/OPENGL_COLOR.h>
#include <PhysBAM_Rendering/PhysBAM_OpenGL/OpenGL/OPENGL_SHAPES.h>
#include <PhysBAM_Rendering/PhysBAM_OpenGL/OpenGL_Components/OPENGL_COMPONENT_THIN_SHELLS_DEBUGGING_3D.h>
using namespace PhysBAM;
//#####################################################################
// Constructor
//#####################################################################
template<class T,class RW> OPENGL_COMPONENT_THIN_SHELLS_DEBUGGING_3D<T,RW>::
OPENGL_COMPONENT_THIN_SHELLS_DEBUGGING_3D(const GRID<TV> &grid,const std::string& directory)
    :OPENGL_COMPONENT("Thin Shells Debugging"),grid(grid),
    invalid_color_map(OPENGL_COLOR::Red()),
    opengl_density_valid_mask(grid,density_valid_mask,&invalid_color_map,OPENGL_SCALAR_FIELD_3D<T,bool>::DRAW_POINTS),
    directory(directory),frame_loaded(-1),valid(false),
    draw_density_valid_mask(false),draw_node_neighbors_visible(false),draw_face_corners_visible(false)
{
    is_animation=true;
    mac_grid=grid.Get_MAC_Grid();u_grid=grid.Get_X_Face_Grid();v_grid=grid.Get_Y_Face_Grid();w_grid=grid.Get_Z_Face_Grid();
}
//#####################################################################
// Function Valid_Frame
//#####################################################################
template<class T,class RW> bool OPENGL_COMPONENT_THIN_SHELLS_DEBUGGING_3D<T,RW>::
Valid_Frame(int frame_input) const
{
    // TODO: make more accurate
    return false;
}
//#####################################################################
// Function Set_Frame
//#####################################################################
template<class T,class RW> void OPENGL_COMPONENT_THIN_SHELLS_DEBUGGING_3D<T,RW>::
Set_Frame(int frame_input)
{
    OPENGL_COMPONENT::Set_Frame(frame_input);
    Reinitialize();
}
//#####################################################################
// Function Set_Draw
//#####################################################################
template<class T,class RW> void OPENGL_COMPONENT_THIN_SHELLS_DEBUGGING_3D<T,RW>::
Set_Draw(bool draw_input)
{
    OPENGL_COMPONENT::Set_Draw(draw_input);
    Reinitialize();
}
//#####################################################################
// Function Display
//#####################################################################
template<class T,class RW> void OPENGL_COMPONENT_THIN_SHELLS_DEBUGGING_3D<T,RW>::
Display(const int in_color) const
{
    OPENGL_COLOR node_neighbor_not_visible_color=OPENGL_COLOR::Magenta(0.5,0.5);
    OPENGL_COLOR face_corners_not_visible_from_face_center_color=OPENGL_COLOR::Magenta(1);
    if(valid && draw){
        if(draw_node_neighbors_visible || draw_face_corners_visible){
            glPushAttrib(GL_ENABLE_BIT | GL_LIGHTING_BIT | GL_LINE_BIT | GL_CURRENT_BIT);

            if (slice && slice->Is_Slice_Mode()) slice->Enable_Clip_Planes();

            glDisable(GL_LIGHTING);

            if(draw_face_corners_visible){
                face_corners_not_visible_from_face_center_color.Send_To_GL_Pipeline();
                glLineWidth(1);
                OpenGL_Begin(GL_LINES);
                for(int i=face_corners_visible_from_face_center_u.domain.min_corner.x;i<=face_corners_visible_from_face_center_u.domain.max_corner.x;i++) for(int j=face_corners_visible_from_face_center_u.domain.min_corner.y;j<=face_corners_visible_from_face_center_u.domain.max_corner.y;j++) for(int k=face_corners_visible_from_face_center_u.domain.min_corner.z;k<=face_corners_visible_from_face_center_u.domain.max_corner.z;k++){
                    if(!face_corners_visible_from_face_center_u(i,j,k)(1).x){OpenGL_Line(u_grid.X(i,j,k),grid.X(i,j,k));}
                    if(!face_corners_visible_from_face_center_u(i,j,k)(2).x){OpenGL_Line(u_grid.X(i,j,k),grid.X(i,j+1,k));}
                    if(!face_corners_visible_from_face_center_u(i,j,k)(3).x){OpenGL_Line(u_grid.X(i,j,k),grid.X(i,j,k+1));}
                    if(!face_corners_visible_from_face_center_u(i,j,k)(4).x){OpenGL_Line(u_grid.X(i,j,k),grid.X(i,j+1,k+1));}
                }
                for(int i=face_corners_visible_from_face_center_v.domain.min_corner.x;i<=face_corners_visible_from_face_center_v.domain.max_corner.x;i++) for(int j=face_corners_visible_from_face_center_v.domain.min_corner.y;j<=face_corners_visible_from_face_center_v.domain.max_corner.y;j++) for(int k=face_corners_visible_from_face_center_v.domain.min_corner.z;k<=face_corners_visible_from_face_center_v.domain.max_corner.z;k++){
                    if(!face_corners_visible_from_face_center_v(i,j,k)(1).x){OpenGL_Line(v_grid.X(i,j,k),grid.X(i,j,k));}
                    if(!face_corners_visible_from_face_center_v(i,j,k)(2).x){OpenGL_Line(v_grid.X(i,j,k),grid.X(i+1,j,k));}
                    if(!face_corners_visible_from_face_center_v(i,j,k)(3).x){OpenGL_Line(v_grid.X(i,j,k),grid.X(i,j,k+1));}
                    if(!face_corners_visible_from_face_center_v(i,j,k)(4).x){OpenGL_Line(v_grid.X(i,j,k),grid.X(i+1,j,k+1));}
                }
                for(int i=face_corners_visible_from_face_center_w.domain.min_corner.x;i<=face_corners_visible_from_face_center_w.domain.max_corner.x;i++) for(int j=face_corners_visible_from_face_center_w.domain.min_corner.y;j<=face_corners_visible_from_face_center_w.domain.max_corner.y;j++) for(int k=face_corners_visible_from_face_center_w.domain.min_corner.z;k<=face_corners_visible_from_face_center_w.domain.max_corner.z;k++){
                    if(!face_corners_visible_from_face_center_w(i,j,k)(1).x){OpenGL_Line(w_grid.X(i,j,k),grid.X(i,j,k));}
                    if(!face_corners_visible_from_face_center_w(i,j,k)(2).x){OpenGL_Line(w_grid.X(i,j,k),grid.X(i+1,j,k));}
                    if(!face_corners_visible_from_face_center_w(i,j,k)(3).x){OpenGL_Line(w_grid.X(i,j,k),grid.X(i,j+1,k));}
                    if(!face_corners_visible_from_face_center_w(i,j,k)(4).x){OpenGL_Line(w_grid.X(i,j,k),grid.X(i+1,j+1,k));}
                }
                OpenGL_End();
            }

            if(draw_node_neighbors_visible){
                node_neighbor_not_visible_color.Send_To_GL_Pipeline();
                glLineWidth(5);
                OpenGL_Begin(GL_LINES);
                for(int i=node_neighbors_visible.domain.min_corner.x;i<=node_neighbors_visible.domain.max_corner.x;i++) for(int j=node_neighbors_visible.domain.min_corner.y;j<=node_neighbors_visible.domain.max_corner.y;j++) for(int k=node_neighbors_visible.domain.min_corner.z;k<=node_neighbors_visible.domain.max_corner.z;k++){
                    if(!node_neighbors_visible(i,j,k)(1)){OpenGL_Line(grid.X(i,j,k),grid.X(i+1,j,k));}
                    if(!node_neighbors_visible(i,j,k)(2)){OpenGL_Line(grid.X(i,j,k),grid.X(i,j+1,k));}
                    if(!node_neighbors_visible(i,j,k)(3)){OpenGL_Line(grid.X(i,j,k),grid.X(i,j,k+1)); }
                }
                OpenGL_End();
            }

            glPopAttrib();
        }

        if(draw_density_valid_mask && density_valid_mask.counts.x) 
            opengl_density_valid_mask.Display(in_color);
    }
}
//#####################################################################
// Function Reinitialize
//#####################################################################
template<class T,class RW> void OPENGL_COMPONENT_THIN_SHELLS_DEBUGGING_3D<T,RW>::
Reinitialize(bool force)
{
    if(force || (draw && (!valid || (is_animation && (frame_loaded != frame)) || (!is_animation && (frame_loaded < 0)))))
    {
        valid = false;

        if(draw_node_neighbors_visible || draw_face_corners_visible){
            std::string tmp_filename = FILE_UTILITIES::Get_Frame_Filename(directory+"/%d/thin_shells_grid_visibility",frame);
            if (FILE_UTILITIES::File_Exists(tmp_filename))
                FILE_UTILITIES::Read_From_File<RW>(tmp_filename,node_neighbors_visible,face_corners_visible_from_face_center_u,face_corners_visible_from_face_center_v,face_corners_visible_from_face_center_w);
        }

        if(draw_density_valid_mask){
            std::string tmp_filename = FILE_UTILITIES::Get_Frame_Filename(directory+"/%d/density_valid_mask",frame);
            if (FILE_UTILITIES::File_Exists(tmp_filename)){
                FILE_UTILITIES::Read_From_File<RW>(tmp_filename,density_valid_mask);
                for(int i=density_valid_mask.domain.min_corner.x;i<=density_valid_mask.domain.max_corner.x;i++) for(int j=density_valid_mask.domain.min_corner.y;j<=density_valid_mask.domain.max_corner.y;j++) for(int k=density_valid_mask.domain.min_corner.z;k<=density_valid_mask.domain.max_corner.z;k++)
                    density_valid_mask(i,j,k)=!density_valid_mask(i,j,k); // negate
                opengl_density_valid_mask.Update();}
            else density_valid_mask.Clean_Memory();
        }

        frame_loaded=frame;
        valid=true;
    }
}
//#####################################################################
// Function Toggle_Draw_Grid_Visibility_Mode
//#####################################################################
template<class T,class RW> void OPENGL_COMPONENT_THIN_SHELLS_DEBUGGING_3D<T,RW>::
Toggle_Draw_Grid_Visibility_Mode()
{
    int mode=((int)draw_node_neighbors_visible<<1) + ((int)draw_face_corners_visible);
    mode=(mode+1)%4;
    draw_node_neighbors_visible=(mode&0x02)!=0;
    draw_face_corners_visible=(mode&0x01)!=0;
    Reinitialize(true);
}
//#####################################################################
// Function Toggle_Draw_Density_Valid_Mask
//#####################################################################
template<class T,class RW> void OPENGL_COMPONENT_THIN_SHELLS_DEBUGGING_3D<T,RW>::
Toggle_Draw_Density_Valid_Mask()
{
    draw_density_valid_mask=!draw_density_valid_mask;
    Reinitialize(true);
}
//#####################################################################
// Function Bounding_Box
//#####################################################################
template<class T,class RW> RANGE<VECTOR<float,3> > OPENGL_COMPONENT_THIN_SHELLS_DEBUGGING_3D<T,RW>::
Bounding_Box() const
{
    if (valid && draw) return RANGE<VECTOR<float,3> >(grid.domain);
    else return RANGE<VECTOR<float,3> >::Centered_Box();
}
//#####################################################################
// Function Set_Slice
//#####################################################################
template<class T,class RW> void OPENGL_COMPONENT_THIN_SHELLS_DEBUGGING_3D<T,RW>::
Set_Slice(OPENGL_SLICE *slice_input)
{
    slice=slice_input;
    opengl_density_valid_mask.Set_Slice(slice_input);
}
//#####################################################################
// Function Slice_Has_Changed
//#####################################################################
template<class T,class RW> void OPENGL_COMPONENT_THIN_SHELLS_DEBUGGING_3D<T,RW>::
Slice_Has_Changed()
{
    opengl_density_valid_mask.Slice_Has_Changed();
}
template class OPENGL_COMPONENT_THIN_SHELLS_DEBUGGING_3D<float,float>;
#ifndef COMPILE_WITHOUT_DOUBLE_SUPPORT
template class OPENGL_COMPONENT_THIN_SHELLS_DEBUGGING_3D<double,double>;
#endif
