//#####################################################################
// Copyright 2004-2009, Eran Guendelman, Andrew Selle.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#include <PhysBAM_Tools/Grids_Uniform_Arrays/ARRAYS_ND.h>
#include <PhysBAM_Tools/Math_Tools/RANGE.h>
#include <PhysBAM_Rendering/PhysBAM_OpenGL/OpenGL/OPENGL_FACE_SCALAR_FIELD_2D.h>
#include <PhysBAM_Rendering/PhysBAM_OpenGL/OpenGL/OPENGL_GRID_2D.h>
#include <PhysBAM_Rendering/PhysBAM_OpenGL/OpenGL/OPENGL_POINTS_2D.h>
#include <PhysBAM_Rendering/PhysBAM_OpenGL/OpenGL/OPENGL_WORLD.h>
namespace PhysBAM{
//#####################################################################
// OPENGL_FACE_SCALAR_FIELD_2D
//#####################################################################
template<class T,class T2> OPENGL_FACE_SCALAR_FIELD_2D<T,T2>::
OPENGL_FACE_SCALAR_FIELD_2D(const GRID<TV> &grid_input,ARRAY<T2,FACE_INDEX<2> > &face_values_input,OPENGL_COLOR_MAP<T2> *color_map_input)
    :grid(grid_input),face_values(face_values_input),x_face_values(face_values.Component(1)),y_face_values(face_values.Component(2)),
    color_map(color_map_input), opengl_points(*new ARRAY<VECTOR<T,2> >)
{
    PHYSBAM_ASSERT(color_map);
}
//#####################################################################
// ~OPENGL_FACE_SCALAR_FIELD_2D
//#####################################################################
template<class T,class T2> OPENGL_FACE_SCALAR_FIELD_2D<T,T2>::
~OPENGL_FACE_SCALAR_FIELD_2D()
{
    delete &opengl_points.points;
}
//#####################################################################
// Display
//#####################################################################
template<class T,class T2> void OPENGL_FACE_SCALAR_FIELD_2D<T,T2>::
Display(const int in_color) const
{
    if (x_face_values.counts.x == 0 || y_face_values.counts.x == 0) return;

    glMatrixMode(GL_MODELVIEW);
    glPushMatrix();
    Send_Transform_To_GL_Pipeline();
    opengl_points.Display(in_color);
    glPopMatrix();
}
//#####################################################################
// Bounding_Box
//#####################################################################
template<class T,class T2> RANGE<VECTOR<float,3> > OPENGL_FACE_SCALAR_FIELD_2D<T,T2>::
Bounding_Box() const
{
    return World_Space_Box(RANGE<VECTOR<float,3> >(grid.domain.min_corner.x,grid.domain.max_corner.x,grid.domain.min_corner.y,grid.domain.max_corner.y,0,0));
}
//#####################################################################
// Update
//#####################################################################
template<class T,class T2> void OPENGL_FACE_SCALAR_FIELD_2D<T,T2>::
Update()
{
    opengl_points.points.Resize(x_face_values.counts.Product()+y_face_values.counts.Product());
    int index=1;
    VECTOR<int,2> index_start,index_end;
    index_start=VECTOR<int,2>(x_face_values.domain.min_corner.x,x_face_values.domain.min_corner.y);
    index_end=VECTOR<int,2>(x_face_values.domain.max_corner.x,x_face_values.domain.max_corner.y);
    for(int i=index_start.x;i<=index_end.x;i++) for(int j=index_start.y;j<=index_end.y;j++) {
        opengl_points.points(index)=grid.X_Face(i,j);
        opengl_points.Set_Point_Color(index,color_map->Lookup(x_face_values(i,j)));
        index++;
    }
    index_start=VECTOR<int,2>(y_face_values.domain.min_corner.x,y_face_values.domain.min_corner.y);
    index_end=VECTOR<int,2>(y_face_values.domain.max_corner.x,y_face_values.domain.max_corner.y);
    for(int i=index_start.x;i<=index_end.x;i++) for(int j=index_start.y;j<=index_end.y;j++) {
        opengl_points.points(index)=grid.Y_Face(i,j);
        opengl_points.Set_Point_Color(index,color_map->Lookup(y_face_values(i,j)));
        index++;
    }
    opengl_points.points.Resize(index-1);
}
//#####################################################################
// Print_Selection_Info
//#####################################################################
template<class T,class T2> void OPENGL_FACE_SCALAR_FIELD_2D<T,T2>::
Print_Selection_Info(std::ostream& output_stream,OPENGL_SELECTION* selection) const
{
    // TODO: this should also interpolate to particles
    if(selection && selection->type==OPENGL_SELECTION::GRID_CELL_2D && grid.Is_MAC_Grid()){
        VECTOR<int,2> index=((OPENGL_SELECTION_GRID_CELL_2D<T>*)selection)->index;
        T2 left=x_face_values(index.x,index.y);
        T2 right=x_face_values(index.x+1,index.y);
        T2 bottom=y_face_values(index.x,index.y);
        T2 top=y_face_values(index.x,index.y+1);
        output_stream<<"    left = "<<left<<",right = "<<right<<std::endl;
        output_stream<<"    bottom = "<<bottom<<",top = "<<top<<std::endl;}
}
//#####################################################################
// Update
//#####################################################################
template<> void OPENGL_FACE_SCALAR_FIELD_2D<float,bool>::
Update()
{
    opengl_points.points.Resize(x_face_values.counts.Product()+y_face_values.counts.Product());
    opengl_points.color=color_map->Lookup(true);
    int index=1;
    VECTOR<int,2> index_start,index_end;
    index_start=VECTOR<int,2>(x_face_values.domain.min_corner.x,x_face_values.domain.min_corner.y);
    index_end=VECTOR<int,2>(x_face_values.domain.max_corner.x,x_face_values.domain.max_corner.y);
    for(int i=index_start.x;i<=index_end.x;i++) for(int j=index_start.y;j<=index_end.y;j++)
        if(x_face_values(i,j)) opengl_points.points(index++)=grid.X_Face(i,j);
    index_start=VECTOR<int,2>(y_face_values.domain.min_corner.x,y_face_values.domain.min_corner.y);
    index_end=VECTOR<int,2>(y_face_values.domain.max_corner.x,y_face_values.domain.max_corner.y);
    for(int i=index_start.x;i<=index_end.x;i++) for(int j=index_start.y;j<=index_end.y;j++)
        if(y_face_values(i,j)) opengl_points.points(index++)=grid.Y_Face(i,j);
    opengl_points.points.Resize(index-1);
}
#ifndef COMPILE_WITHOUT_DOUBLE_SUPPORT
//#####################################################################
// Update
//#####################################################################
template<> void OPENGL_FACE_SCALAR_FIELD_2D<double,bool>::
Update()
{
    opengl_points.points.Resize(x_face_values.counts.Product()+y_face_values.counts.Product());
    opengl_points.color=color_map->Lookup(true);
    int index=1;
    VECTOR<int,2> index_start,index_end;
    index_start=VECTOR<int,2>(x_face_values.domain.min_corner.x,x_face_values.domain.min_corner.y);
    index_end=VECTOR<int,2>(x_face_values.domain.max_corner.x,x_face_values.domain.max_corner.y);
    for(int i=index_start.x;i<=index_end.x;i++) for(int j=index_start.y;j<=index_end.y;j++)
        if(x_face_values(i,j)) opengl_points.points(index++)=grid.X_Face(i,j);
    index_start=VECTOR<int,2>(y_face_values.domain.min_corner.x,y_face_values.domain.min_corner.y);
    index_end=VECTOR<int,2>(y_face_values.domain.max_corner.x,y_face_values.domain.max_corner.y);
    for(int i=index_start.x;i<=index_end.x;i++) for(int j=index_start.y;j<=index_end.y;j++)
        if(y_face_values(i,j)) opengl_points.points(index++)=grid.Y_Face(i,j);
    opengl_points.points.Resize(index-1);
}
//#####################################################################
#endif
template class OPENGL_FACE_SCALAR_FIELD_2D<float,int>;
template class OPENGL_FACE_SCALAR_FIELD_2D<float,bool>;
template class OPENGL_FACE_SCALAR_FIELD_2D<float,float>;
#ifndef COMPILE_WITHOUT_DOUBLE_SUPPORT
template class OPENGL_FACE_SCALAR_FIELD_2D<double,int>;
template class OPENGL_FACE_SCALAR_FIELD_2D<double,bool>;
template class OPENGL_FACE_SCALAR_FIELD_2D<double,double>;
#endif
}
