#ifndef COMPILE_WITHOUT_DYADIC_SUPPORT
//#####################################################################
// Copyright 2004-2005, Geoffrey Irving, Frank Losasso.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#include <PhysBAM_Tools/Grids_Dyadic/QUADTREE_GRID.h>
#include <PhysBAM_Geometry/Basic_Geometry/ORIENTED_BOX.h>
#include <PhysBAM_Geometry/Red_Green/RED_GREEN_GRID_2D.h>
#include <PhysBAM_Rendering/PhysBAM_OpenGL/OpenGL/OPENGL_ADAPTIVE_NODE_SCALAR_FIELD.h>
#include <PhysBAM_Rendering/PhysBAM_OpenGL/OpenGL/OPENGL_QUADTREE_GRID.h>
namespace PhysBAM{
//#####################################################################
// Function Display_Quadtree
//#####################################################################
template<class T> void Display_Quadtree(ARRAY<T>& value,OPENGL_COLOR_MAP<T>& color_map,QUADTREE_GRID<T>& grid,ARRAY<QUADTREE_CELL<T>*>& cell_pointer_from_index,const bool smooth_shading)
{
    glPushAttrib(GL_LIGHTING_BIT | GL_TEXTURE_BIT | GL_DEPTH_BUFFER_BIT);
    glDisable(GL_LIGHTING);glDisable(GL_TEXTURE_2D);glDepthMask(GL_FALSE);

    int node_ordering[4]={2,3,1,0};
    if(smooth_shading) glShadeModel(GL_SMOOTH); else glShadeModel(GL_FLAT);
    OpenGL_Begin(GL_QUADS);
    if(smooth_shading){
        for(int i=1;i<=cell_pointer_from_index.m;i++)
            if(cell_pointer_from_index(i)&&!cell_pointer_from_index(i)->Has_Children()){
                for(int j=0;j<4;j++){
                    int node_index=cell_pointer_from_index(i)->Node(node_ordering[j]);
                    color_map.Lookup(value(node_index)).Send_To_GL_Pipeline();
                    OpenGL_Vertex(grid.Node_Location(node_index));}}}
    else{
        for(int i=1;i<=cell_pointer_from_index.m;i++)
            if(cell_pointer_from_index(i)&&!cell_pointer_from_index(i)->Has_Children()){
                T avg_value=0;
                for(int j=0;j<4;j++)avg_value+=value(cell_pointer_from_index(i)->Node(node_ordering[j]));
                avg_value*=.25;
                for(int j=0;j<4;j++){
                    int node_index=cell_pointer_from_index(i)->Node(node_ordering[j]);
                    color_map.Lookup(avg_value).Send_To_GL_Pipeline();
                    OpenGL_Vertex(grid.Node_Location(node_index));}}}
    OpenGL_End();
    glPopAttrib();
}
//#####################################################################
// Function Display_Red_Green
//#####################################################################
template<class T> void Display_Red_Green(ARRAY<T>& value,OPENGL_COLOR_MAP<T>& color_map,ARRAY<VECTOR<T,2> >& node_locations,ARRAY<RED_TRIANGLE<T>*>& cell_pointer_from_index)
{
    glPushAttrib(GL_LIGHTING_BIT | GL_TEXTURE_BIT);
    glDisable(GL_LIGHTING);glDisable(GL_TEXTURE_2D);

    glShadeModel(GL_SMOOTH);
    OpenGL_Begin(GL_TRIANGLES);
    for(int i=1;i<=cell_pointer_from_index.m;i++){
        RED_TRIANGLE<T>* triangle=cell_pointer_from_index(i);
        if(triangle&&!triangle->Has_Red_Children()){
            if(!triangle->Has_Green_Children())
                for(int n=2;n>=0;n--){
                    int node_index=triangle->Node(n);
                    color_map.Lookup(value(node_index)).Send_To_GL_Pipeline();
                    OpenGL_Vertex(node_locations(node_index));}
            else for(int t=1;t<=triangle->green_children->elements.m;t++)
                for(int n=3;n>=1;n--){
                    int node_index=triangle->green_children->elements(t)(n);
                    color_map.Lookup(value(node_index)).Send_To_GL_Pipeline();
                    OpenGL_Vertex(node_locations(node_index));}}}
    OpenGL_End();
    glPopAttrib();
}
//#####################################################################
// Function Display
//#####################################################################
template<> void OPENGL_ADAPTIVE_NODE_SCALAR_FIELD<QUADTREE_GRID<float> >::Display(const int in_color) const
{
    Display_Quadtree(value,*color_map,grid,grid.Cell_Pointer_From_Index(),smooth_shading);
}
template<> void OPENGL_ADAPTIVE_NODE_SCALAR_FIELD<RED_GREEN_GRID_2D<float> >::Display(const int in_color) const
{
    Display_Red_Green(value,*color_map,grid.Node_Locations(),grid.Cell_Pointer_From_Index());
}
#ifndef COMPILE_WITHOUT_DOUBLE_SUPPORT
template<> void OPENGL_ADAPTIVE_NODE_SCALAR_FIELD<QUADTREE_GRID<double> >::Display(const int in_color) const
{
    Display_Quadtree(value,*color_map,grid,grid.Cell_Pointer_From_Index(),smooth_shading);
}
template<> void OPENGL_ADAPTIVE_NODE_SCALAR_FIELD<RED_GREEN_GRID_2D<double> >::Display(const int in_color) const
{
    Display_Red_Green(value,*color_map,grid.Node_Locations(),grid.Cell_Pointer_From_Index());
}
#endif
//#####################################################################
// Function Bounding_Box
//#####################################################################
template<class T_GRID> RANGE<VECTOR<float,3> > OPENGL_ADAPTIVE_NODE_SCALAR_FIELD<T_GRID>::
Bounding_Box() const
{
    return RANGE<VECTOR<float,3> >(grid.uniform_grid.domain.min_corner.x,grid.uniform_grid.domain.max_corner.x,grid.uniform_grid.domain.min_corner.y,grid.uniform_grid.domain.max_corner.y,0,0);
}
//#####################################################################
// Function Bounding_Box
//#####################################################################
template<class T_GRID> void OPENGL_ADAPTIVE_NODE_SCALAR_FIELD<T_GRID>::
Print_Selection_Info(std::ostream& output_stream,OPENGL_SELECTION* current_selection) const
{
    if(current_selection && current_selection->type==OPENGL_SELECTION::QUADTREE_NODE){
        int index=((OPENGL_SELECTION_QUADTREE_NODE<T>*)current_selection)->index;
        if(value.Valid_Index(index))
            output_stream<<value(index);}
    output_stream<<std::endl;
}
//#####################################################################
template class OPENGL_ADAPTIVE_NODE_SCALAR_FIELD<QUADTREE_GRID<float> >;
template class OPENGL_ADAPTIVE_NODE_SCALAR_FIELD<RED_GREEN_GRID_2D<float> >;
#ifndef COMPILE_WITHOUT_DOUBLE_SUPPORT
template class OPENGL_ADAPTIVE_NODE_SCALAR_FIELD<QUADTREE_GRID<double> >;
template class OPENGL_ADAPTIVE_NODE_SCALAR_FIELD<RED_GREEN_GRID_2D<double> >;
#endif
}
#endif
