//#####################################################################
// Copyright 2010, Jon Gretarsson.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#if !COMPILE_WITHOUT_DYADIC_SUPPORT || COMPILE_WITH_BINTREE_SUPPORT
#include <PhysBAM_Tools/Grids_Dyadic/DYADIC_GRID_ITERATOR_CELL.h>
#include <PhysBAM_Tools/Grids_Uniform/UNIFORM_GRID_ITERATOR_CELL.h>
#include <PhysBAM_Geometry/Basic_Geometry/ORIENTED_BOX.h>
#include <PhysBAM_Rendering/PhysBAM_OpenGL/OpenGL/OPENGL_BINTREE_GRID.h>
#include <PhysBAM_Rendering/PhysBAM_OpenGL/OpenGL/OPENGL_PREFERENCES.h>
using namespace PhysBAM;
//#####################################################################
// Display
//#####################################################################
template<class T> void OPENGL_BINTREE_GRID<T>::
Display(const int in_color) const
{
    glMatrixMode(GL_MODELVIEW);glPushMatrix();
    Send_Transform_To_GL_Pipeline();

    glPushAttrib(GL_ENABLE_BIT);glDisable(GL_LIGHTING);glDisable(GL_DEPTH_TEST);
    color.Send_To_GL_Pipeline();
    GLint mode;glGetIntegerv(GL_RENDER_MODE, &mode);

    if(mode==GL_SELECT){
        glPushAttrib(GL_ENABLE_BIT | GL_POINT_BIT);
        glDisable(GL_CULL_FACE);

        // Draw grid cells for selection
        glPushName(1);
        glPushName(0); // overwritten per cell in Draw_Cells_For_Selection_Helper

        for(DYADIC_GRID_ITERATOR_CELL<BINTREE_GRID<T> > iterator(grid,grid.number_of_ghost_cells);iterator.Valid();iterator.Next()){
            glLoadName(iterator.Cell_Index());
            OpenGL_Begin(GL_QUADS);
            BINTREE_CELL<T>* cell=iterator.Cell_Pointer();
            OpenGL_Quad_2D(grid.Node_Location(0,cell).Insert((T)-.05,2),grid.Node_Location(1,cell).Insert((T).05,2));
            OpenGL_End();}

        glPopName();

        glLoadName(2);
        glPushName(0); // overwritten per node below
        glPointSize(OPENGL_PREFERENCES::selection_point_size);
        for(int i=1;i<=grid.number_of_nodes;i++){
            glLoadName(i);
            OpenGL_Begin(GL_POINTS);
            OpenGL_Vertex(grid.Node_Location(i));
            OpenGL_End();}
        glPopName();
        glPopName();
        glPopAttrib();
    }
    else{
        // draw the coarse uniform uniform_grid
        OpenGL_Begin(GL_LINES);
        for(UNIFORM_GRID_ITERATOR_CELL<VECTOR<T,1> > iterator(grid.uniform_grid,grid.number_of_ghost_cells);iterator.Valid();iterator.Next())
            OpenGL_Line(iterator.Bounding_Box().min_corner.Insert((T)0,2).Insert((T)0,3),iterator.Bounding_Box().max_corner.Insert((T)0,2).Insert((T)0,3));
        OpenGL_End();

        (color*.5).Send_To_GL_Pipeline();
        // draw the actual tree
        OpenGL_Begin(GL_LINES);
        for(DYADIC_GRID_ITERATOR_CELL<BINTREE_GRID<T> > iterator(grid,grid.number_of_ghost_cells);iterator.Valid();iterator.Next()){
            BINTREE_CELL<T>* cell=iterator.Cell_Pointer();
            OpenGL_Line(grid.Node_Location(0,cell).Insert((T)-.05,2),grid.Node_Location(1,cell).Insert((T)-.05,2));
            OpenGL_Line(grid.Node_Location(1,cell).Insert((T)-.05,2),grid.Node_Location(1,cell).Insert((T) .05,2));
            OpenGL_Line(grid.Node_Location(1,cell).Insert((T) .05,2),grid.Node_Location(0,cell).Insert((T) .05,2));
            OpenGL_Line(grid.Node_Location(0,cell).Insert((T) .05,2),grid.Node_Location(0,cell).Insert((T)-.05,2));}
        OpenGL_End();

        // Highlight current selection
        if(current_selection){
            if(current_selection->type==OPENGL_SELECTION::BINTREE_CELL){
                OPENGL_SELECTION_BINTREE_CELL<T>* real_selection=(OPENGL_SELECTION_BINTREE_CELL<T>*)current_selection;
                const BINTREE_CELL<T>* cell=grid.Cell_Pointer_From_Index()(real_selection->index);
                OPENGL_SELECTION::Draw_Highlighted_Quad(grid.Node_Location(0,cell).Insert((T)-.05,2),grid.Node_Location(1,cell).Insert((T).05,2));}
            else if(current_selection->type==OPENGL_SELECTION::BINTREE_NODE){
                OPENGL_SELECTION_BINTREE_NODE<T>* real_selection=(OPENGL_SELECTION_BINTREE_NODE<T>*)current_selection;
                OPENGL_SELECTION::Draw_Highlighted_Vertex(real_selection->location);}}
    }
    glPopAttrib();
    glPopMatrix();
}
//#####################################################################
// Bounding_Box
//#####################################################################
template<class T> RANGE<VECTOR<float,3> > OPENGL_BINTREE_GRID<T>::
Bounding_Box() const
{
    RANGE<VECTOR<float,3> > box(grid.uniform_grid.domain.min_corner.x,grid.uniform_grid.domain.max_corner.x,-1,1,0,0);
    return World_Space_Box(box);
}
//#####################################################################
// Print_Selection_Info
//#####################################################################
template<class T> void OPENGL_BINTREE_GRID<T>::
Print_Selection_Info(std::ostream& output_stream,OPENGL_SELECTION* selection) const
{
    if(current_selection && current_selection->type==OPENGL_SELECTION::BINTREE_NODE){
        int index=((OPENGL_SELECTION_BINTREE_NODE<T>*)current_selection)->index;
        output_stream<<"Selected node "<<index<<" ("<<grid.Node_Location(index)<<")"<<std::endl;}
    else if(current_selection && current_selection->type==OPENGL_SELECTION::BINTREE_CELL){
        int index=((OPENGL_SELECTION_BINTREE_CELL<T>*)current_selection)->index;
        BINTREE_CELL<T>* cell=grid.Cell_Pointer_From_Index()(index);
        output_stream<<"Selected cell "<<index<<" ("<<cell->Center()<<")"<<std::endl;
        output_stream<<"Faces: "<<cell->Face(0)<<", "<<cell->Face(1)<<std::endl;
        output_stream<<"DX = "<<cell->DX()<<std::endl;
        output_stream<<"depth = "<<cell->Depth_Of_This_Cell()<<std::endl;}
}
//#####################################################################
// Get_Selection
//#####################################################################
template<class T> OPENGL_SELECTION *OPENGL_BINTREE_GRID<T>::
Get_Selection(GLuint *buffer,int buffer_size)
{
    OPENGL_SELECTION *selection=0;
    if(buffer_size==2){
        if(buffer[0]==1)
            selection=new OPENGL_SELECTION_BINTREE_CELL<T>(this,buffer[1]);
        else if(buffer[0]==2)
            selection=new OPENGL_SELECTION_BINTREE_NODE<T>(this,buffer[1],grid.Node_Location(buffer[1]));}
    return selection;
}
//#####################################################################
// Highlight_Selection
//#####################################################################
template<class T> void OPENGL_BINTREE_GRID<T>::
Highlight_Selection(OPENGL_SELECTION *selection)
{
    delete current_selection;current_selection=0;
    if (selection->type==OPENGL_SELECTION::BINTREE_CELL){
        OPENGL_SELECTION_BINTREE_CELL<T> *real_selection=(OPENGL_SELECTION_BINTREE_CELL<T>*)selection;
        current_selection=new OPENGL_SELECTION_BINTREE_CELL<T>(this,real_selection->index);}
    else if (selection->type == OPENGL_SELECTION::BINTREE_NODE){
        OPENGL_SELECTION_BINTREE_NODE<T> *real_selection=(OPENGL_SELECTION_BINTREE_NODE<T>*)selection;
        current_selection=new OPENGL_SELECTION_BINTREE_NODE<T>(this,real_selection->index,real_selection->location);}
}
//#####################################################################
// Clear_highlight
//#####################################################################
template<class T> void OPENGL_BINTREE_GRID<T>::
Clear_Highlight()
{
    delete current_selection;current_selection=0;
}
//#####################################################################
// Create_Or_Destroy_Selection_After_Frame_Change
//#####################################################################
template<class T> OPENGL_SELECTION* OPENGL_BINTREE_GRID<T>::
Create_Or_Destroy_Selection_After_Frame_Change(OPENGL_SELECTION* old_selection,bool& delete_selection)
{
    if(old_selection && old_selection->object==this)
        return Get_Updated_Selection(old_selection);
    return 0;
}
//#####################################################################
// Get_Updated_Selection
//#####################################################################
template<class T> OPENGL_SELECTION* OPENGL_BINTREE_GRID<T>::
Get_Updated_Selection(OPENGL_SELECTION *selection)
{
    if (selection->type == OPENGL_SELECTION::BINTREE_NODE){
        OPENGL_SELECTION_BINTREE_NODE<T> *real_selection=(OPENGL_SELECTION_BINTREE_NODE<T>*)selection;
        BINTREE_CELL<T>* leaf=grid.Clamped_Leaf_Cell(real_selection->location,1e-5);
        T tolerance_squared=sqr(1e-5);
        int new_node_index=leaf->Nearest_Node(real_selection->location);
        VECTOR<T,1> new_node_location=grid.Node_Location(new_node_index);
        if((real_selection->location-new_node_location).Magnitude_Squared()<tolerance_squared) return new OPENGL_SELECTION_BINTREE_NODE<T>(this,new_node_index,new_node_location);}
    return 0;
}
//#####################################################################
// OPENGL_SELECTION_BINTREE_CELL
//#####################################################################
template<class T> RANGE<VECTOR<float,3> > OPENGL_SELECTION_BINTREE_CELL<T>::
Bounding_Box() const
{
    PHYSBAM_ASSERT(object);
    const BINTREE_GRID<T> &grid=((OPENGL_BINTREE_GRID<T>*)object)->grid;
    const PhysBAM::BINTREE_CELL<T> *cell=grid.Cell_Pointer_From_Index()(index);
    return object->World_Space_Box(RANGE<VECTOR<float,2> >(cell->Bounding_Box()));
}
//#####################################################################
// OPENGL_SELECTION_BINTREE_NODE
//#####################################################################
template<class T> RANGE<VECTOR<float,3> > OPENGL_SELECTION_BINTREE_NODE<T>::
Bounding_Box() const
{
    PHYSBAM_ASSERT(object);
    RANGE<VECTOR<float,1> > box((VECTOR<float,1>(location)));
    return object->World_Space_Box(box);
}
//#####################################################################
template class OPENGL_BINTREE_GRID<float>;
#ifndef COMPILE_WITHOUT_DOUBLE_SUPPORT
template class OPENGL_BINTREE_GRID<double>;
#endif
#endif
