#ifndef COMPILE_WITHOUT_RLE_SUPPORT
//#####################################################################
// Copyright 2005, Geoffrey Irving.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#include <PhysBAM_Tools/Grids_RLE/RLE_GRID_ITERATOR_CELL_2D.h>
#include <PhysBAM_Tools/Grids_RLE/RLE_GRID_ITERATOR_CELL_3D.h>
#include <PhysBAM_Geometry/Basic_Geometry/ORIENTED_BOX.h>
#include <PhysBAM_Rendering/PhysBAM_OpenGL/OpenGL/OPENGL_PREFERENCES.h>
#include <PhysBAM_Rendering/PhysBAM_OpenGL/OpenGL/OPENGL_RLE_GRID_2D.h>
using namespace PhysBAM;
//#####################################################################
// Function Display
//#####################################################################
template<class T> void OPENGL_RLE_GRID_2D<T>::
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
        for(CELL_ITERATOR cell(grid,grid.number_of_ghost_cells);cell;cell++){
            glPushName(cell.i);glPushName(cell.j);glPushName(cell.run-&grid.columns(cell.i)(1));
            VECTOR<T,2> X1=grid.uniform_grid.X(cell.i,cell.j),X2=grid.uniform_grid.X(cell.i+1,cell.jmax());
            OpenGL_Begin(GL_QUADS);OpenGL_Quad_2D(X1,X2);OpenGL_End();
            glPopName();glPopName();glPopName();}
        glPopName();
        glPopAttrib();}
    else{
        (color*.5).Send_To_GL_Pipeline();
        // draw the actual tree
        OpenGL_Begin(GL_LINES);
        for(CELL_ITERATOR cell(grid,grid.number_of_ghost_cells);cell;cell++){
            VECTOR<T,2> X1=grid.uniform_grid.X(cell.i,cell.j),X2=grid.uniform_grid.X(cell.i+1,cell.jmax());
            OpenGL_Line(VECTOR<T,2>(X1.x,X1.y),VECTOR<T,2>(X2.x,X1.y));OpenGL_Line(VECTOR<T,2>(X2.x,X1.y),VECTOR<T,2>(X2.x,X2.y));
            OpenGL_Line(VECTOR<T,2>(X2.x,X2.y),VECTOR<T,2>(X1.x,X2.y));OpenGL_Line(VECTOR<T,2>(X1.x,X2.y),VECTOR<T,2>(X1.x,X1.y));}
        OpenGL_End();

        glPushAttrib(GL_LINE_BIT);glLineWidth(4*OPENGL_PREFERENCES::line_width);OpenGL_Begin(GL_LINE_LOOP);
        OpenGL_Quad_2D(grid.uniform_grid.domain.min_corner,grid.uniform_grid.domain.max_corner);
        OpenGL_End();glPopAttrib();

        // Highlight current selection
        if(current_selection && current_selection->type==OPENGL_SELECTION::RLE_CELL_2D){
            OPENGL_SELECTION_RLE_CELL_2D<T>* real_selection=(OPENGL_SELECTION_RLE_CELL_2D<T>*)current_selection;
            VECTOR<int,2> I=real_selection->I;const RLE_RUN_2D* run=real_selection->run;
            OPENGL_SELECTION::Draw_Highlighted_Quad(grid.uniform_grid.X(I),grid.uniform_grid.X(I.x+1,run->is_long?(run+1)->jmin:I.y+1));}}
    glPopAttrib();
    glPopMatrix();
}
//#####################################################################
// Function Bounding_Box
//#####################################################################
template<class T> RANGE<VECTOR<float,3> > OPENGL_RLE_GRID_2D<T>::
Bounding_Box() const
{
    RANGE<VECTOR<float,3> > box(grid.uniform_grid.domain.min_corner.x,grid.uniform_grid.domain.max_corner.x,grid.uniform_grid.domain.min_corner.y,grid.uniform_grid.domain.max_corner.y,0,0);
    return World_Space_Box(box);
}
//#####################################################################
// Function Get_Selection
//#####################################################################
template<class T> OPENGL_SELECTION* OPENGL_RLE_GRID_2D<T>::
Get_Selection(GLuint* buffer,int buffer_size)
{
    OPENGL_SELECTION* selection=0;
    if(buffer_size==4 && buffer[0]==1) selection=new OPENGL_SELECTION_RLE_CELL_2D<T>(this,VECTOR<int,2>(buffer[1],buffer[2]),&grid.columns(buffer[1])(1)+buffer[3]);
    return selection;
}
//#####################################################################
// Function Create_Or_Destroy_Selection_After_Frame_Change
//#####################################################################
template<class T> OPENGL_SELECTION* OPENGL_RLE_GRID_2D<T>::
Create_Or_Destroy_Selection_After_Frame_Change(OPENGL_SELECTION* old_selection,bool& delete_selection)
{
    if(old_selection && old_selection->object==this){
        OPENGL_SELECTION_RLE_CELL_2D<T>* selection=(OPENGL_SELECTION_RLE_CELL_2D<T>*)current_selection;
        selection->run=grid.Clamped_Run_In_Column(selection->I);
        if(selection->run->is_long) selection->I.y=selection->run->jmin;
        Highlight_Selection(selection);
    }
    return 0;
}
//#####################################################################
// Function Highlight_Selection
//#####################################################################
template<class T> void OPENGL_RLE_GRID_2D<T>::
Highlight_Selection(OPENGL_SELECTION *selection)
{
    delete current_selection;current_selection=0;
    if(selection->type==OPENGL_SELECTION::RLE_CELL_2D){
        OPENGL_SELECTION_RLE_CELL_2D<T>* real_selection=(OPENGL_SELECTION_RLE_CELL_2D<T>*)selection;
        current_selection=new OPENGL_SELECTION_RLE_CELL_2D<T>(this,real_selection->I,real_selection->run);}
}
//#####################################################################
// Function Clear_Highlight
//#####################################################################
template<class T> void OPENGL_RLE_GRID_2D<T>::
Clear_Highlight()
{
    delete current_selection;current_selection=0;
}
//#####################################################################
// Function Print_Selection
//#####################################################################
template<class T> void OPENGL_RLE_GRID_2D<T>::
Print_Selection_Info(std::ostream& output_stream,OPENGL_SELECTION* current_selection) const
{
    if(current_selection && current_selection->type==OPENGL_SELECTION::RLE_CELL_2D){
        OPENGL_SELECTION_RLE_CELL_2D<T>* selection=(OPENGL_SELECTION_RLE_CELL_2D<T>*)current_selection;
        const RLE_RUN_2D* run=selection->run;int dj=selection->I.y-run->jmin,cell=run->cell+dj;
        int fx1=run->faces[0]+dj,fx2=run->faces[1]+dj,fy1=cell,fy2=fy1+(run->is_long?2:1);
        assert(grid.long_run_cells==2);
        bool left_face_long=false,right_face_long=false;
        if(run->is_long && grid.long_run_faces_horizontal==2){ // TODO: make more efficient (maybe store some extra data in selection object to avoid this lookup)
            const RLE_RUN_2D* left_run=grid.Clamped_Run_In_Column(selection->I.x-1,selection->I.y);if(left_run->is_long && (left_run+1)->jmin>selection->I.y+1) left_face_long=true;
            const RLE_RUN_2D* right_run=grid.Clamped_Run_In_Column(selection->I.x+1,selection->I.y);if(right_run->is_long && (right_run+1)->jmin>selection->I.y+1) right_face_long=true;}

        output_stream<<"Selected cell "<<cell<<", I = "<<selection->I<<", center = "<<typename RLE_GRID_2D<T>::CELL_ITERATOR(grid,selection->I).Center()<<std::endl;
        output_stream<<"run sign = "<<run->is_long<<", index = "<<run-&grid.columns(selection->I.x)(1)<<", jmin = "<<run->jmin<<", jmax = "<<(run+1)->jmin<<std::endl;
        output_stream<<"faces x = "<<fx1<<(left_face_long?"(long)":"")<<", "<<fx2<<(right_face_long?"(long)":"")<<", y = "<<fy1<<", "<<fy2<<std::endl;}
}
//#####################################################################
// OPENGL_SELECTION_RLE_CELL_2D
//#####################################################################
template<class T> RANGE<VECTOR<float,3> > OPENGL_SELECTION_RLE_CELL_2D<T>::
Bounding_Box() const
{
    PHYSBAM_ASSERT(object);
    const RLE_GRID_2D<T>& grid=((OPENGL_RLE_GRID_2D<T>*)object)->grid;
    return object->World_Space_Box(RANGE<VECTOR<float,2> >(VECTOR<float,2>(grid.uniform_grid.X(I)),VECTOR<float,2>(grid.uniform_grid.X(I.x+1,run->is_long?(run+1)->jmin:I.y+1))));
}
//#####################################################################
template class OPENGL_RLE_GRID_2D<float>;
#ifndef COMPILE_WITHOUT_DOUBLE_SUPPORT
template class OPENGL_RLE_GRID_2D<double>;
#endif
#endif
