#ifndef COMPILE_WITHOUT_DYADIC_SUPPORT
//#####################################################################
// Copyright 2004, Eran Guendelman.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#include <PhysBAM_Tools/Grids_Dyadic/DYADIC_GRID_ITERATOR_FACE.h>
#include <PhysBAM_Tools/Grids_Dyadic/OCTREE_GRID.h>
#include <PhysBAM_Geometry/Basic_Geometry/ORIENTED_BOX.h>
#include <PhysBAM_Rendering/PhysBAM_OpenGL/OpenGL/OPENGL_COLOR_MAP.h>
#include <PhysBAM_Rendering/PhysBAM_OpenGL/OpenGL/OPENGL_OCTREE_FACE_SCALAR_FIELD.h>
#include <PhysBAM_Rendering/PhysBAM_OpenGL/OpenGL/OPENGL_OCTREE_GRID.h>
using namespace PhysBAM;

//#####################################################################
// Display
//#####################################################################
template<class T,class T2> void OPENGL_OCTREE_FACE_SCALAR_FIELD<T,T2>::
Display(const int in_color) const
{
    glPushAttrib(GL_LIGHTING_BIT | GL_TEXTURE_BIT | GL_POINT_BIT | GL_ENABLE_BIT);
    glPointSize(point_size);glDisable(GL_LIGHTING);glDisable(GL_TEXTURE_2D);

    if(slice && slice->Is_Slice_Mode())slice->Enable_Clip_Planes();

    if(draw_points){
        OpenGL_Begin(GL_POINTS);
        for(DYADIC_GRID_ITERATOR_FACE<OCTREE_GRID<T> > iterator(grid,grid.Map_All_Faces());iterator.Valid();iterator.Next()){
            color_map->Lookup(value(iterator.Face_Index())).Send_To_GL_Pipeline();OpenGL_Vertex(iterator.Location());}
        OpenGL_End();}
    else{
        OpenGL_Begin(GL_LINES);
        for(DYADIC_GRID_ITERATOR_FACE<OCTREE_GRID<T> > iterator(grid,grid.Map_All_Faces());iterator.Valid();iterator.Next()){
            T2 val=value(iterator.Face_Index());
            color_map->Lookup(val).Send_To_GL_Pipeline();
            OpenGL_Line(iterator.Location(),VECTOR<T,3>(iterator.Location()+val*line_size*VECTOR<T,3>::Axis_Vector(iterator.Axis())));}
        OpenGL_End();}

    glPopAttrib();
}
//#####################################################################
// Bounding_Box
//#####################################################################
template<class T,class T2> RANGE<VECTOR<float,3> > OPENGL_OCTREE_FACE_SCALAR_FIELD<T,T2>::
Bounding_Box() const
{
    return World_Space_Box(grid.uniform_grid.domain);
}
//#####################################################################
// Print_Selection_Info
//#####################################################################
template<class T,class T2> void OPENGL_OCTREE_FACE_SCALAR_FIELD<T,T2>::
Print_Selection_Info(std::ostream& output_stream,OPENGL_SELECTION* current_selection) const
{
    if(current_selection && current_selection->type==OPENGL_SELECTION::OCTREE_FACE){
        int index=((OPENGL_SELECTION_OCTREE_FACE<T>*)current_selection)->index;
        output_stream<<value(index);}
    output_stream<<std::endl;
}
//#####################################################################
template class OPENGL_OCTREE_FACE_SCALAR_FIELD<float>;
template class OPENGL_OCTREE_FACE_SCALAR_FIELD<float,int>;
template class OPENGL_OCTREE_FACE_SCALAR_FIELD<float,bool>;
#ifndef COMPILE_WITHOUT_DOUBLE_SUPPORT
template class OPENGL_OCTREE_FACE_SCALAR_FIELD<double>;
template class OPENGL_OCTREE_FACE_SCALAR_FIELD<double,int>;
template class OPENGL_OCTREE_FACE_SCALAR_FIELD<double,bool>;
#endif
#endif
