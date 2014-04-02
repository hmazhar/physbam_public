//#####################################################################
// Copyright 2010, Jon Gretarsson.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#if !COMPILE_WITHOUT_DYADIC_SUPPORT || COMPILE_WITH_BINTREE_SUPPORT
#include <PhysBAM_Tools/Arrays_Computations/ARRAY_MIN_MAX.h>
#include <PhysBAM_Tools/Grids_Dyadic/DYADIC_GRID_ITERATOR_FACE.h>
#include <PhysBAM_Geometry/Basic_Geometry/ORIENTED_BOX.h>
#include <PhysBAM_Rendering/PhysBAM_OpenGL/OpenGL/OPENGL_BINTREE_FACE_SCALAR_FIELD.h>
#include <PhysBAM_Rendering/PhysBAM_OpenGL/OpenGL/OPENGL_BINTREE_GRID.h>
using namespace PhysBAM;

//#####################################################################
// Display
//#####################################################################
namespace {
template<class T,class T2> 
void Display_Helper(const OPENGL_BINTREE_FACE_SCALAR_FIELD<T,T2>& self,const int in_color,const T2& value,const VECTOR<T,1>& location) {
    self.color_map->Lookup(value).Send_To_GL_Pipeline();OpenGL_Vertex(location.Insert((T)value,2));
}

void Display_Helper(const OPENGL_BINTREE_FACE_SCALAR_FIELD<float,bool>& self,const int in_color,const bool& value,const VECTOR<float,1>& location) {
    if(value){self.color_map->Lookup(value).Send_To_GL_Pipeline();OpenGL_Vertex(location);}
}

void Display_Helper(const OPENGL_BINTREE_FACE_SCALAR_FIELD<double,bool>& self,const int in_color,const bool& value,const VECTOR<double,1>& location) {
    if(value){self.color_map->Lookup(value).Send_To_GL_Pipeline();OpenGL_Vertex(location);}
}
};
template<class T,class T2> void OPENGL_BINTREE_FACE_SCALAR_FIELD<T,T2>::
Display(const int in_color) const
{
    glPushAttrib(GL_LIGHTING_BIT | GL_TEXTURE_BIT | GL_POINT_BIT);
    glPointSize(point_size);glDisable(GL_LIGHTING);glDisable(GL_TEXTURE_2D);

    if(draw_points) OpenGL_Begin(GL_POINTS); else OpenGL_Begin(GL_LINE_STRIP);
    for(DYADIC_GRID_ITERATOR_FACE<BINTREE_GRID<T> > iterator(grid,grid.Map_All_Faces());iterator.Valid();iterator.Next()){
        Display_Helper(*this,in_color,value(iterator.Face_Index()),iterator.Location());}
    OpenGL_End();

    glPopAttrib();
}
//#####################################################################
// Bounding_Box
//#####################################################################
template<class T,class T2> RANGE<VECTOR<float,3> > OPENGL_BINTREE_FACE_SCALAR_FIELD<T,T2>::
Bounding_Box() const
{
    return RANGE<VECTOR<float,3> >(grid.uniform_grid.domain.min_corner.x,grid.uniform_grid.domain.max_corner.x,ARRAYS_COMPUTATIONS::Min(value),ARRAYS_COMPUTATIONS::Max(value),0,0);
}
//#####################################################################
// Print_Selection_Info
//#####################################################################
template<class T,class T2> void OPENGL_BINTREE_FACE_SCALAR_FIELD<T,T2>::
Print_Selection_Info(std::ostream& output_stream,OPENGL_SELECTION* current_selection) const
{
    if(current_selection && current_selection->type==OPENGL_SELECTION::BINTREE_CELL){
        int index=((OPENGL_SELECTION_BINTREE_CELL<T>*)current_selection)->index;
        BINTREE_CELL<T>* cell=grid.Cell_Pointer_From_Index()(index);
        T u_left,u_right;
        u_left=value(cell->Face(0));
        u_right=value(cell->Face(1));
        VECTOR<T,1> center_velocity(0.5*(u_left+u_right));
        output_stream<<"    u left = "<<u_left<<",right = "<<u_right<<std::endl;
        output_stream<<"    center velocity = "<<center_velocity<<std::endl;
        T ux=(u_right-u_left)/cell->DX().x;
        output_stream<<"divergence = "<<ux<<" (ux="<<ux<<")"<<std::endl;}
}
//#####################################################################
template class OPENGL_BINTREE_FACE_SCALAR_FIELD<float>;
template class OPENGL_BINTREE_FACE_SCALAR_FIELD<float,int>;
template class OPENGL_BINTREE_FACE_SCALAR_FIELD<float,bool>;
#ifndef COMPILE_WITHOUT_DOUBLE_SUPPORT
template class OPENGL_BINTREE_FACE_SCALAR_FIELD<double>;
template class OPENGL_BINTREE_FACE_SCALAR_FIELD<double,int>;
template class OPENGL_BINTREE_FACE_SCALAR_FIELD<double,bool>;
#endif
#endif
