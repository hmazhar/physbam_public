#ifndef COMPILE_WITHOUT_RLE_SUPPORT
//#####################################################################
// Copyright 2005, Geoffrey Irving.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#include <PhysBAM_Tools/Grids_RLE/RLE_GRID_ITERATOR_CELL_3D.h>
#include <PhysBAM_Rendering/PhysBAM_OpenGL/OpenGL/OPENGL_RLE_CELL_SCALAR_FIELD_3D.h>
#include <PhysBAM_Rendering/PhysBAM_OpenGL/OpenGL/OPENGL_RLE_GRID_3D.h>
#include <PhysBAM_Rendering/PhysBAM_OpenGL/OpenGL/OPENGL_RLE_SLICE.h>
using namespace PhysBAM;
//#####################################################################
// Function Display
//#####################################################################
template<class T,class T2> void OPENGL_RLE_CELL_SCALAR_FIELD_3D<T,T2>::
Display(const int in_color) const
{
    OPENGL_RLE_SLICE* slice=(OPENGL_RLE_SLICE*)this->slice;

    RANGE<VECTOR<int,2> > region(grid.columns.domain.min_corner.x+1,grid.columns.domain.max_corner.x-1,grid.columns.domain.min_corner.y+1,grid.columns.domain.max_corner.y-1);
    if(slice && slice->mode==OPENGL_SLICE::CELL_SLICE){
        switch(slice->axis){
            case 1: region.min_corner.x=region.max_corner.x=slice->index;break;
            case 2: break; // TODO: deal with y slices
            case 3: region.min_corner.y=region.max_corner.y=slice->index;break;}}

    if(!draw_filled){
        glPushAttrib(GL_LIGHTING_BIT | GL_TEXTURE_BIT | GL_POINT_BIT);
        glPointSize(point_size);glDisable(GL_LIGHTING);
        OpenGL_Begin(GL_POINTS);
        if(grid.long_run_cells==1)
            for(CELL_ITERATOR cell(grid,region);cell;cell++){int c=cell.Cell();if(value(c)){
                color_map->Lookup(value(c)).Send_To_GL_Pipeline();
                OpenGL_Vertex(cell.Center());}}
        else if(grid.long_run_cells==2)
            for(CELL_ITERATOR cell(grid,region);cell;cell++){int c=cell.Cell();
                if(value(c)){
                    color_map->Lookup(value(c)).Send_To_GL_Pipeline();
                    OpenGL_Vertex(grid.uniform_grid.Center(cell.i,cell.j,cell.ij));}
                if(cell.Long() && value(c+1)){
                    color_map->Lookup(value(c+1)).Send_To_GL_Pipeline();
                    OpenGL_Vertex(grid.uniform_grid.Center(cell.i,cell.jmax()-1,cell.ij));}}
        else PHYSBAM_NOT_IMPLEMENTED();
        OpenGL_End();glPopAttrib();}
    else if(grid.long_run_cells==2){
        if(slice && slice->mode==OPENGL_SLICE::CELL_SLICE && slice->axis!=2){
            glPushAttrib(GL_LIGHTING_BIT | GL_TEXTURE_BIT | GL_DEPTH_BUFFER_BIT | GL_ENABLE_BIT);
            glDisable(GL_LIGHTING);glDepthMask(GL_FALSE);glDisable(GL_CULL_FACE);
            OpenGL_Begin(GL_QUADS);
            if(slice->axis==1){
                T x=grid.uniform_grid.Axis_X_plus_half(slice->index,slice->axis);
                for(CELL_ITERATOR cell(grid,region);cell;cell++){
                    color_map->Lookup(value(cell.Cell())).Send_To_GL_Pipeline();
                    OpenGL_Quad(VECTOR<T,3>(x,grid.uniform_grid.Axis_X(cell.j,2),grid.uniform_grid.Axis_X(cell.ij,3)),VECTOR<T,3>(0,0,grid.uniform_grid.dX.z),
                        VECTOR<T,3>(0,grid.uniform_grid.Axis_X(cell.jmax(),2)-grid.uniform_grid.Axis_X(cell.j,2),0));}}
            else if(slice->axis==3){
                T z=grid.uniform_grid.Axis_X_plus_half(slice->index,slice->axis);
                for(CELL_ITERATOR cell(grid,region);cell;cell++){
                    color_map->Lookup(value(cell.Cell())).Send_To_GL_Pipeline();
                    OpenGL_Quad(VECTOR<T,3>(grid.uniform_grid.Axis_X(cell.i,1),grid.uniform_grid.Axis_X(cell.j,2),z),VECTOR<T,3>(grid.uniform_grid.dX.x,0,0),
                        VECTOR<T,3>(0,grid.uniform_grid.Axis_X(cell.jmax(),2)-grid.uniform_grid.Axis_X(cell.j,2),0));}}
            OpenGL_End();glPopAttrib();}}
    else PHYSBAM_NOT_IMPLEMENTED();
}
//#####################################################################
// Function Bounding_Box
//#####################################################################
template<class T,class T2> RANGE<VECTOR<float,3> > OPENGL_RLE_CELL_SCALAR_FIELD_3D<T,T2>::
Bounding_Box() const
{
    return (RANGE<VECTOR<float,3> >)grid.uniform_grid.domain;
}
//#####################################################################
// Print_Selection_Info
//#####################################################################
template<class T,class T2> void OPENGL_RLE_CELL_SCALAR_FIELD_3D<T,T2>::
Print_Selection_Info(std::ostream& output_stream,OPENGL_SELECTION* current_selection) const
{
    if(current_selection && current_selection->type==OPENGL_SELECTION::RLE_CELL_3D){
        OPENGL_SELECTION_RLE_CELL_3D<T>* selection=(OPENGL_SELECTION_RLE_CELL_3D<T>*)current_selection;
        const RLE_RUN_3D* run=selection->run;int dj=selection->I.y-run->jmin,cell=run->cell+dj;
        output_stream<<"value = "<<value(cell)<<std::endl;}
}
//#####################################################################
template class OPENGL_RLE_CELL_SCALAR_FIELD_3D<float>;
template class OPENGL_RLE_CELL_SCALAR_FIELD_3D<float,int>;
template class OPENGL_RLE_CELL_SCALAR_FIELD_3D<float,bool>;
#ifndef COMPILE_WITHOUT_DOUBLE_SUPPORT
template class OPENGL_RLE_CELL_SCALAR_FIELD_3D<double>;
template class OPENGL_RLE_CELL_SCALAR_FIELD_3D<double,int>;
template class OPENGL_RLE_CELL_SCALAR_FIELD_3D<double,bool>;
#endif
#endif
