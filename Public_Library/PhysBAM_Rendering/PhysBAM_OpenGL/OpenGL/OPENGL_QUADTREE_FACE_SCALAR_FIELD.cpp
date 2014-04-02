#ifndef COMPILE_WITHOUT_DYADIC_SUPPORT
//#####################################################################
// Copyright 2004-2009, Frank Losasso, Andrew Selle.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#include <PhysBAM_Tools/Grids_Dyadic/DYADIC_GRID_ITERATOR_FACE.h>
#include <PhysBAM_Geometry/Basic_Geometry/ORIENTED_BOX.h>
#include <PhysBAM_Rendering/PhysBAM_OpenGL/OpenGL/OPENGL_QUADTREE_FACE_SCALAR_FIELD.h>
#include <PhysBAM_Rendering/PhysBAM_OpenGL/OpenGL/OPENGL_QUADTREE_GRID.h>
using namespace PhysBAM;

//#####################################################################
// Display
//#####################################################################
template<class T,class T2> void OPENGL_QUADTREE_FACE_SCALAR_FIELD<T,T2>::
Display(const int in_color) const
{
    glPushAttrib(GL_LIGHTING_BIT | GL_TEXTURE_BIT | GL_POINT_BIT);
    glPointSize(point_size);glDisable(GL_LIGHTING);glDisable(GL_TEXTURE_2D);

    if(draw_points){
        OpenGL_Begin(GL_POINTS);
        for(DYADIC_GRID_ITERATOR_FACE<QUADTREE_GRID<T> > iterator(grid,grid.Map_All_Faces());iterator.Valid();iterator.Next()){
            color_map->Lookup(value(iterator.Face_Index())).Send_To_GL_Pipeline();OpenGL_Vertex(iterator.Location());}
        OpenGL_End();}
    else{
        OpenGL_Begin(GL_LINES);
        for(DYADIC_GRID_ITERATOR_FACE<QUADTREE_GRID<T> > iterator(grid,grid.Map_All_Faces());iterator.Valid();iterator.Next()){
            T2 val=value(iterator.Face_Index());
            color_map->Lookup(val).Send_To_GL_Pipeline();
            OpenGL_Line(iterator.Location(),iterator.Location()+val*line_size*VECTOR<T,2>::Axis_Vector(iterator.Axis()));}
        OpenGL_End();}

    glPopAttrib();
}
//#####################################################################
// Bounding_Box
//#####################################################################
template<class T,class T2> RANGE<VECTOR<float,3> > OPENGL_QUADTREE_FACE_SCALAR_FIELD<T,T2>::
Bounding_Box() const
{
    return RANGE<VECTOR<float,3> >(grid.uniform_grid.domain.min_corner.x,grid.uniform_grid.domain.max_corner.x,grid.uniform_grid.domain.min_corner.y,grid.uniform_grid.domain.max_corner.y,0,0);
}
//#####################################################################
// Print_Selection_Info
//#####################################################################
template<class T,class T2> void OPENGL_QUADTREE_FACE_SCALAR_FIELD<T,T2>::
Print_Selection_Info(std::ostream& output_stream,OPENGL_SELECTION* current_selection) const
{
    if(current_selection && current_selection->type==OPENGL_SELECTION::QUADTREE_CELL){
        int index=((OPENGL_SELECTION_QUADTREE_CELL<T>*)current_selection)->index;
        QUADTREE_CELL<T>* cell=grid.Cell_Pointer_From_Index()(index);
        T u_left,u_right,v_bottom,v_top;
        u_left=value(cell->Face(0));
        u_right=value(cell->Face(1));
        v_bottom=value(cell->Face(2));
        v_top=value(cell->Face(3));
        VECTOR<T,2> center_velocity(0.5*(u_left+u_right),0.5*(v_bottom+v_top));
        output_stream<<"    u left = "<<u_left<<",right = "<<u_right<<std::endl;
        output_stream<<"    v bottom = "<<v_bottom<<",top = "<<v_top<<std::endl;
        output_stream<<"    center velocity = "<<center_velocity<<std::endl;
        T ux=(u_right-u_left)/cell->DX().x,vy=(v_top-v_bottom)/cell->DX().y;
        output_stream<<"divergence = "<<ux+vy<<" (ux="<<ux<<", vy="<<vy<<")"<<std::endl;}
}
//#####################################################################
template class OPENGL_QUADTREE_FACE_SCALAR_FIELD<float>;
template class OPENGL_QUADTREE_FACE_SCALAR_FIELD<float,int>;
template class OPENGL_QUADTREE_FACE_SCALAR_FIELD<float,bool>;
#ifndef COMPILE_WITHOUT_DOUBLE_SUPPORT
template class OPENGL_QUADTREE_FACE_SCALAR_FIELD<double>;
template class OPENGL_QUADTREE_FACE_SCALAR_FIELD<double,int>;
template class OPENGL_QUADTREE_FACE_SCALAR_FIELD<double,bool>;
#endif
#endif
