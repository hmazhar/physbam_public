#ifndef COMPILE_WITHOUT_RLE_SUPPORT
//#####################################################################
// Copyright 2005, Geoffrey Irving.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#include <PhysBAM_Tools/Grids_RLE/RLE_GRID_ITERATOR_CELL_2D.h>
#include <PhysBAM_Rendering/PhysBAM_OpenGL/OpenGL/OPENGL_RLE_CELL_SCALAR_FIELD_2D.h>
#include <PhysBAM_Rendering/PhysBAM_OpenGL/OpenGL/OPENGL_RLE_GRID_2D.h>
#include <PhysBAM_Rendering/PhysBAM_OpenGL/OpenGL/OPENGL_SELECTION.h>
using namespace PhysBAM;
//#####################################################################
// Function Display
//#####################################################################
template<class T,class T2> void OPENGL_RLE_CELL_SCALAR_FIELD_2D<T,T2>::
Display(const int in_color) const
{
    if(!draw_filled){
        glPushAttrib(GL_LIGHTING_BIT | GL_TEXTURE_BIT | GL_POINT_BIT);
        glPointSize(point_size);glDisable(GL_LIGHTING);glDisable(GL_TEXTURE_2D);
        OpenGL_Begin(GL_POINTS);
        if(grid.long_run_cells==1)
            for(CELL_ITERATOR cell(grid,grid.number_of_ghost_cells);cell;cell++){int c=cell.Cell();if(value(c)){
                color_map->Lookup(value(c)).Send_To_GL_Pipeline();
                OpenGL_Vertex(VECTOR<T,2>(grid.uniform_grid.Axis_X_plus_half(cell.i,1),(T).5*(grid.uniform_grid.Axis_X(cell.j,2)+grid.uniform_grid.Axis_X(cell.jmax(),2))));}}
        else if(grid.long_run_cells==2)
            for(CELL_ITERATOR cell(grid,grid.number_of_ghost_cells);cell;cell++){int c=cell.Cell();
                if(value(c)){
                    color_map->Lookup(value(c)).Send_To_GL_Pipeline();
                    OpenGL_Vertex(grid.uniform_grid.Center(cell.i,cell.j));}
                if(cell.Long() && value(c+1)){
                    color_map->Lookup(value(c+1)).Send_To_GL_Pipeline();
                    OpenGL_Vertex(grid.uniform_grid.Center(cell.i,cell.jmax()-1));}}
        else PHYSBAM_NOT_IMPLEMENTED();
        OpenGL_End();glPopAttrib();}
    else if(grid.long_run_cells==2){
        glPushAttrib(GL_LIGHTING_BIT | GL_TEXTURE_BIT | GL_DEPTH_BUFFER_BIT);
        glDisable(GL_LIGHTING);glDisable(GL_TEXTURE_2D);glDepthMask(GL_FALSE);
        OpenGL_Begin(GL_QUADS);
        for(CELL_ITERATOR cell(grid,grid.number_of_ghost_cells);cell;cell++){
            color_map->Lookup(value(cell.Cell())).Send_To_GL_Pipeline();
            OpenGL_Quad_2D(grid.uniform_grid.X(cell.i,cell.j),grid.uniform_grid.X(cell.i+1,cell.jmax()));}
        OpenGL_End();glPopAttrib();}
    else PHYSBAM_NOT_IMPLEMENTED();
}
//#####################################################################
// Function Bounding_Box
//#####################################################################
template<class T,class T2> RANGE<VECTOR<float,3> > OPENGL_RLE_CELL_SCALAR_FIELD_2D<T,T2>::
Bounding_Box() const
{
    return RANGE<VECTOR<float,3> >(grid.uniform_grid.domain.min_corner.x,grid.uniform_grid.domain.max_corner.x,grid.uniform_grid.domain.min_corner.y,grid.uniform_grid.domain.max_corner.y,0,0);
}
//#####################################################################
// Function Print_Selection_Info
//#####################################################################
template<class T,class T2> void OPENGL_RLE_CELL_SCALAR_FIELD_2D<T,T2>::
Print_Selection_Info(std::ostream& output_stream,OPENGL_SELECTION* current_selection) const
{
    if(current_selection && current_selection->type==OPENGL_SELECTION::RLE_CELL_2D){
        OPENGL_SELECTION_RLE_CELL_2D<T>* selection=(OPENGL_SELECTION_RLE_CELL_2D<T>*)current_selection;
        const RLE_RUN_2D* run=selection->run;int dj=selection->I.y-run->jmin,cell=run->cell+dj;
        assert(grid.long_run_cells==2);
        output_stream<<"value = "<<value(cell)<<std::endl;}
}
//#####################################################################
template class OPENGL_RLE_CELL_SCALAR_FIELD_2D<float>;
template class OPENGL_RLE_CELL_SCALAR_FIELD_2D<float,int>;
template class OPENGL_RLE_CELL_SCALAR_FIELD_2D<float,bool>;
#ifndef COMPILE_WITHOUT_DOUBLE_SUPPORT
template class OPENGL_RLE_CELL_SCALAR_FIELD_2D<double>;
template class OPENGL_RLE_CELL_SCALAR_FIELD_2D<double,int>;
template class OPENGL_RLE_CELL_SCALAR_FIELD_2D<double,bool>;
#endif
#endif
