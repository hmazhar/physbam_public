//#####################################################################
// Copyright 2004-2009, Eran Guendelman, Andrew Selle.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#include <PhysBAM_Tools/Grids_Uniform_Arrays/ARRAYS_ND.h>
#include <PhysBAM_Tools/Grids_Uniform_Arrays/ARRAYS_ND_VIEW.h>
#include <PhysBAM_Tools/Grids_Uniform_Arrays/FACE_ARRAYS.h>
#include <PhysBAM_Tools/Math_Tools/RANGE.h>
#include <PhysBAM_Rendering/PhysBAM_OpenGL/OpenGL/OPENGL_FACE_SCALAR_FIELD_3D.h>
#include <PhysBAM_Rendering/PhysBAM_OpenGL/OpenGL/OPENGL_GRID_3D.h>
#include <PhysBAM_Rendering/PhysBAM_OpenGL/OpenGL/OPENGL_POINTS_3D.h>
#include <PhysBAM_Rendering/PhysBAM_OpenGL/OpenGL/OPENGL_UNIFORM_SLICE.h>
#include <PhysBAM_Rendering/PhysBAM_OpenGL/OpenGL/OPENGL_WORLD.h>
namespace PhysBAM{
//#####################################################################
// OPENGL_FACE_SCALAR_FIELD_3D
//#####################################################################
template<class T,class T2> OPENGL_FACE_SCALAR_FIELD_3D<T,T2>::
OPENGL_FACE_SCALAR_FIELD_3D(const GRID<TV> &grid_input,ARRAY<T2,FACE_INDEX<3> > &face_values_input,OPENGL_COLOR_MAP<T2> *color_map_input)
    :grid(grid_input),face_values(face_values_input),x_face_values(face_values_input.Component(1)),y_face_values(face_values_input.Component(2)),z_face_values(face_values_input.Component(3)),
    color_map(color_map_input), scale(1), opengl_points(*new ARRAY<VECTOR<T,3> >)
{
    PHYSBAM_ASSERT(color_map);
}
//#####################################################################
// ~OPENGL_FACE_SCALAR_FIELD_3D
//#####################################################################
template<class T,class T2> OPENGL_FACE_SCALAR_FIELD_3D<T,T2>::
~OPENGL_FACE_SCALAR_FIELD_3D()
{
    delete &opengl_points.points;
}
//#####################################################################
// Display
//#####################################################################
template<class T,class T2> void OPENGL_FACE_SCALAR_FIELD_3D<T,T2>::
Display(const int in_color) const
{
    if (x_face_values.counts.x == 0 || y_face_values.counts.y == 0 || z_face_values.counts.z == 0) return;

    glMatrixMode(GL_MODELVIEW);
    glPushMatrix();
    Send_Transform_To_GL_Pipeline();
    opengl_points.Display(in_color);
    glPopMatrix();
}
//#####################################################################
// Bounding_Box
//#####################################################################
template<class T,class T2> RANGE<VECTOR<float,3> > OPENGL_FACE_SCALAR_FIELD_3D<T,T2>::
Bounding_Box() const
{
    return World_Space_Box((RANGE<VECTOR<float,3> >)grid.domain);
}
//#####################################################################
// Update
//#####################################################################
template<class T,class T2> void OPENGL_FACE_SCALAR_FIELD_3D<T,T2>::
Update()
{
    opengl_points.points.Resize(x_face_values.counts.Product()+y_face_values.counts.Product()+z_face_values.counts.Product());
    int index=1;
    VECTOR<int,3> index_start,index_end;
    OPENGL_UNIFORM_SLICE* slice=(OPENGL_UNIFORM_SLICE*)this->slice;
    OPENGL_UNIFORM_SLICE::Get_Face_Index_Range(slice, x_face_values, 1, index_start, index_end);
    for(int i=index_start.x;i<=index_end.x;i++) for(int j=index_start.y;j<=index_end.y;j++) for(int k=index_start.z;k<=index_end.z;k++) {
        opengl_points.points(index)=grid.X_Face(i,j,k);
        opengl_points.Set_Point_Color(index,color_map->Lookup(x_face_values(i,j,k)));
        index++;
    }
    OPENGL_UNIFORM_SLICE::Get_Face_Index_Range(slice, y_face_values, 2, index_start, index_end);
    for(int i=index_start.x;i<=index_end.x;i++) for(int j=index_start.y;j<=index_end.y;j++) for(int k=index_start.z;k<=index_end.z;k++) {
        opengl_points.points(index)=grid.Y_Face(i,j,k);
        opengl_points.Set_Point_Color(index,color_map->Lookup(y_face_values(i,j,k)));
        index++;
    }
    OPENGL_UNIFORM_SLICE::Get_Face_Index_Range(slice, z_face_values, 3, index_start, index_end);
    for(int i=index_start.x;i<=index_end.x;i++) for(int j=index_start.y;j<=index_end.y;j++) for(int k=index_start.z;k<=index_end.z;k++) {
        opengl_points.points(index)=grid.Z_Face(i,j,k);
        opengl_points.Set_Point_Color(index,color_map->Lookup(z_face_values(i,j,k)));
        index++;
    }
    opengl_points.points.Resize(index-1);
}
//#####################################################################
// Update
//#####################################################################
template<> void OPENGL_FACE_SCALAR_FIELD_3D<float,bool>::
Update()
{
    OPENGL_UNIFORM_SLICE* slice=(OPENGL_UNIFORM_SLICE*)this->slice;
    opengl_points.points.Resize(x_face_values.counts.Product()+y_face_values.counts.Product()+z_face_values.counts.Product());
    opengl_points.color=color_map->Lookup(true);
    int index=1;
    VECTOR<int,3> index_start,index_end;
    OPENGL_UNIFORM_SLICE::Get_Face_Index_Range(slice, x_face_values, 1, index_start, index_end,scale);
    for(int i=index_start.x;i<=index_end.x;i++) for(int j=index_start.y;j<=index_end.y;j++) for(int k=index_start.z;k<=index_end.z;k++)
        if(x_face_values(i,j,k)) opengl_points.points(index++)=grid.X_Face(i,j,k);
    OPENGL_UNIFORM_SLICE::Get_Face_Index_Range(slice, y_face_values, 2, index_start, index_end,scale);
    for(int i=index_start.x;i<=index_end.x;i++) for(int j=index_start.y;j<=index_end.y;j++) for(int k=index_start.z;k<=index_end.z;k++)
        if(y_face_values(i,j,k)) opengl_points.points(index++)=grid.Y_Face(i,j,k);
    OPENGL_UNIFORM_SLICE::Get_Face_Index_Range(slice, z_face_values, 3, index_start, index_end,scale);
    for(int i=index_start.x;i<=index_end.x;i++) for(int j=index_start.y;j<=index_end.y;j++) for(int k=index_start.z;k<=index_end.z;k++)
        if(z_face_values(i,j,k)) opengl_points.points(index++)=grid.Z_Face(i,j,k);
    opengl_points.points.Resize(index-1);
}
//#####################################################################
// Update
//#####################################################################
template<> void OPENGL_FACE_SCALAR_FIELD_3D<double,bool>::
Update()
{
    OPENGL_UNIFORM_SLICE* slice=(OPENGL_UNIFORM_SLICE*)this->slice;
    opengl_points.points.Resize(x_face_values.counts.Product()+y_face_values.counts.Product()+z_face_values.counts.Product());
    opengl_points.color=color_map->Lookup(true);
    int index=1;
    VECTOR<int,3> index_start,index_end;
    OPENGL_UNIFORM_SLICE::Get_Face_Index_Range(slice, x_face_values, 1, index_start, index_end);
    for(int i=index_start.x;i<=index_end.x;i++) for(int j=index_start.y;j<=index_end.y;j++) for(int k=index_start.z;k<=index_end.z;k++)
        if(x_face_values(i,j,k)) opengl_points.points(index++)=grid.X_Face(i,j,k);
    OPENGL_UNIFORM_SLICE::Get_Face_Index_Range(slice, y_face_values, 2, index_start, index_end);
    for(int i=index_start.x;i<=index_end.x;i++) for(int j=index_start.y;j<=index_end.y;j++) for(int k=index_start.z;k<=index_end.z;k++)
        if(y_face_values(i,j,k)) opengl_points.points(index++)=grid.Y_Face(i,j,k);
    OPENGL_UNIFORM_SLICE::Get_Face_Index_Range(slice, z_face_values, 3, index_start, index_end);
    for(int i=index_start.x;i<=index_end.x;i++) for(int j=index_start.y;j<=index_end.y;j++) for(int k=index_start.z;k<=index_end.z;k++)
        if(z_face_values(i,j,k)) opengl_points.points(index++)=grid.Z_Face(i,j,k);
    opengl_points.points.Resize(index-1);
}
//#####################################################################
// Slice_Has_Changed
//#####################################################################
template<class T,class T2> void OPENGL_FACE_SCALAR_FIELD_3D<T,T2>::
Slice_Has_Changed()
{
    Update();
}
//#####################################################################
// Print_Selection_Info
//#####################################################################
template<class T,class T2> void OPENGL_FACE_SCALAR_FIELD_3D<T,T2>::
Print_Selection_Info(std::ostream& output_stream,OPENGL_SELECTION* selection) const
{
    // TODO: this should also interpolate to particles
    if(selection && selection->type==OPENGL_SELECTION::GRID_CELL_3D && grid.Is_MAC_Grid()){
        VECTOR<int,3> index=((OPENGL_SELECTION_GRID_CELL_3D<T>*)selection)->index;
        T2 left=x_face_values(index),right=x_face_values(index.x+1,index.y,index.z),
            bottom=y_face_values(index),top=y_face_values(index.x,index.y+1,index.z),
            back=z_face_values(index),front=z_face_values(index.x,index.y,index.z+1);
        output_stream<<"    left = "<<left<<",right = "<<right<<std::endl;
        output_stream<<"    bottom = "<<bottom<<",top = "<<top<<std::endl;
        output_stream<<"    back = "<<back<<", front = "<<front<<std::endl;}
}
//#####################################################################
template class OPENGL_FACE_SCALAR_FIELD_3D<float,int>;
template class OPENGL_FACE_SCALAR_FIELD_3D<float,bool>;
template class OPENGL_FACE_SCALAR_FIELD_3D<float,float>;
#ifndef COMPILE_WITHOUT_DOUBLE_SUPPORT
template class OPENGL_FACE_SCALAR_FIELD_3D<double,int>;
template class OPENGL_FACE_SCALAR_FIELD_3D<double,bool>;
template class OPENGL_FACE_SCALAR_FIELD_3D<double,double>;
#endif
}
