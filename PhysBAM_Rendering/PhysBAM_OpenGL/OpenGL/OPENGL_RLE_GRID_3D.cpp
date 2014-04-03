#ifndef COMPILE_WITHOUT_RLE_SUPPORT
//#####################################################################
// Copyright 2005, Geoffrey Irving.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#include <PhysBAM_Tools/Grids_RLE/RLE_GRID_ITERATOR_FACE_Y.h>
#include <PhysBAM_Geometry/Basic_Geometry/ORIENTED_BOX.h>
#include <PhysBAM_Rendering/PhysBAM_OpenGL/OpenGL/OPENGL_PREFERENCES.h>
#include <PhysBAM_Rendering/PhysBAM_OpenGL/OpenGL/OPENGL_RLE_GRID_3D.h>
#include <PhysBAM_Rendering/PhysBAM_OpenGL/OpenGL/OPENGL_RLE_SLICE.h>
using namespace PhysBAM;
//#####################################################################
// Function Display
//#####################################################################
template<class T> void OPENGL_RLE_GRID_3D<T>::
Display(const int in_color) const
{
    OPENGL_RLE_SLICE* slice=(OPENGL_RLE_SLICE*)this->slice;

    glMatrixMode(GL_MODELVIEW);glPushMatrix();
    Send_Transform_To_GL_Pipeline();

    glPushAttrib(GL_ENABLE_BIT);glDisable(GL_LIGHTING);glDisable(GL_DEPTH_TEST);
    color.Send_To_GL_Pipeline();
    GLint mode;glGetIntegerv(GL_RENDER_MODE, &mode);

    RANGE<VECTOR<int,2> > region(grid.columns.domain.min_corner.x+1,grid.columns.domain.max_corner.x-1,grid.columns.domain.min_corner.y+1,grid.columns.domain.max_corner.y-1);
    if(slice && slice->mode==OPENGL_SLICE::CELL_SLICE){
        switch(slice->axis){
            case 1: region.min_corner.x=region.max_corner.x=slice->index;break;
            case 2: break; // TODO: deal with y slices
            case 3: region.min_corner.y=region.max_corner.y=slice->index;break;}}

    if(mode==GL_SELECT)
//    if(true)
    {
        if(slice && slice->mode==OPENGL_SLICE::CELL_SLICE && slice->axis!=2){
            glPushAttrib(GL_ENABLE_BIT);
            glDisable(GL_CULL_FACE);
            // Draw grid cells for selection
            glPushName(1);
            VECTOR<T,3> x_vector(grid.uniform_grid.dX.x,0,0),z_vector(0,0,grid.uniform_grid.dX.z),quad_axis,quad_offset;
            if(slice->axis==1){quad_axis=z_vector;quad_offset=x_vector;}else{quad_axis=x_vector;quad_offset=z_vector;}
            for(CELL_ITERATOR cell(grid,region);cell;cell++){
                glPushName(cell.i);glPushName(cell.j);glPushName(cell.ij);glPushName(cell.run-&grid.columns(cell.i,cell.ij)(1));
                VECTOR<T,3> X1=grid.uniform_grid.X(cell.i,cell.j,cell.ij),y_axis(0,grid.uniform_grid.Axis_X(cell.jmax(),2)-grid.uniform_grid.Axis_X(cell.j,2),0);
                OpenGL_Begin(GL_QUADS);OpenGL_Quad(X1,quad_axis,y_axis);OpenGL_Quad(X1+quad_offset,quad_axis,y_axis);OpenGL_End();
                glPopName();glPopName();glPopName();glPopName();}
            glPopName();
            glPopAttrib();}}
    else{
        T dx=grid.uniform_grid.dX.x,dz=grid.uniform_grid.dX.z;
        (color*.5).Send_To_GL_Pipeline();
        // draw the actual tree
        OpenGL_Begin(GL_LINES);
        for(int i=region.min_corner.x;i<=region.max_corner.x;i++)for(int ij=region.min_corner.y;ij<=region.max_corner.y;ij++){
            T miny=grid.uniform_grid.Axis_X(grid.columns(i,ij)(1).jmin,2),maxy=grid.uniform_grid.Axis_X(grid.columns(i,ij)(grid.columns(i,ij).m).jmin,2);
            VECTOR<T,3> X1(grid.uniform_grid.Axis_X(i,1),miny,grid.uniform_grid.Axis_X(ij,3));
            OpenGL_Line(X1,VECTOR<T,3>(X1.x,maxy,X1.z));
            OpenGL_Line(VECTOR<T,3>(X1.x+dx,X1.y,X1.z),VECTOR<T,3>(X1.x+dx,maxy,X1.z));
            OpenGL_Line(VECTOR<T,3>(X1.x+dx,X1.y,X1.z+dz),VECTOR<T,3>(X1.x+dx,maxy,X1.z+dz));
            OpenGL_Line(VECTOR<T,3>(X1.x,X1.y,X1.z+dz),VECTOR<T,3>(X1.x,maxy,X1.z+dz));}
        for(FACE_Y_ITERATOR face(grid,region,true);face;face++){
            VECTOR<T,3> X1=grid.uniform_grid.X(face.i(),face.j(),face.ij());
            OpenGL_Line(X1,VECTOR<T,3>(X1.x+dx,X1.y,X1.z));
            OpenGL_Line(VECTOR<T,3>(X1.x+dx,X1.y,X1.z),VECTOR<T,3>(X1.x+dx,X1.y,X1.z+dz));
            OpenGL_Line(VECTOR<T,3>(X1.x+dx,X1.y,X1.z+dz),VECTOR<T,3>(X1.x,X1.y,X1.z+dz));
            OpenGL_Line(VECTOR<T,3>(X1.x,X1.y,X1.z+dz),X1);}
        OpenGL_End();

        // Highlight current selection
        if(current_selection && current_selection->type==OPENGL_SELECTION::RLE_CELL_3D){
            OPENGL_SELECTION_RLE_CELL_3D<T>* real_selection=(OPENGL_SELECTION_RLE_CELL_3D<T>*)current_selection;
            VECTOR<int,3> I=real_selection->I;const RLE_RUN_3D* run=real_selection->run;
            OPENGL_SELECTION::Draw_Highlighted_Box(grid.uniform_grid.X(I),grid.uniform_grid.X(I.x+1,run->is_long?(run+1)->jmin:I.y+1,I.z+1));}}
    glPopAttrib();
    glPopMatrix();
}
//#####################################################################
// Function Bounding_Box
//#####################################################################
template<class T> RANGE<VECTOR<float,3> > OPENGL_RLE_GRID_3D<T>::
Bounding_Box() const
{
    return World_Space_Box(grid.uniform_grid.domain);
}
//#####################################################################
// Function Get_Selection
//#####################################################################
template<class T> OPENGL_SELECTION* OPENGL_RLE_GRID_3D<T>::
Get_Selection(GLuint* buffer,int buffer_size)
{
    OPENGL_SELECTION* selection=0;
    if(buffer_size==5 && buffer[0]==1) selection=new OPENGL_SELECTION_RLE_CELL_3D<T>(this,VECTOR<int,3>(buffer[1],buffer[2],buffer[3]),&grid.columns(buffer[1],buffer[3])(1)+buffer[4]);
    return selection;
}
//#####################################################################
// Function Highlight_Selection
//#####################################################################
template<class T> void OPENGL_RLE_GRID_3D<T>::
Highlight_Selection(OPENGL_SELECTION *selection)
{
    delete current_selection;current_selection=0;
    if(selection->type==OPENGL_SELECTION::RLE_CELL_3D){
        OPENGL_SELECTION_RLE_CELL_3D<T>* real_selection=(OPENGL_SELECTION_RLE_CELL_3D<T>*)selection;
        current_selection=new OPENGL_SELECTION_RLE_CELL_3D<T>(this,real_selection->I,real_selection->run);}
}
//#####################################################################
// Function Clear_Highlight
//#####################################################################
template<class T> void OPENGL_RLE_GRID_3D<T>::
Clear_Highlight()
{
    delete current_selection;current_selection=0;
}
//#####################################################################
// OPENGL_SELECTION_RLE_CELL_3D
//#####################################################################
template<class T> RANGE<VECTOR<float,3> > OPENGL_SELECTION_RLE_CELL_3D<T>::
Bounding_Box() const
{
    PHYSBAM_ASSERT(object);
    const RLE_GRID_3D<T>& grid=((OPENGL_RLE_GRID_3D<T>*)object)->grid;
    return object->World_Space_Box(RANGE<VECTOR<float,3> >(VECTOR<float,3>(grid.uniform_grid.X(I)),VECTOR<float,3>(grid.uniform_grid.X(I.x+1,run->is_long?(run+1)->jmin:I.y+1,I.z))));
}
//#####################################################################
// Function Print_Selection
//#####################################################################
template<class T> void OPENGL_RLE_GRID_3D<T>::
Print_Selection_Info(std::ostream& output_stream,OPENGL_SELECTION* current_selection) const
{
    if(current_selection && current_selection->type==OPENGL_SELECTION::RLE_CELL_3D){
        OPENGL_SELECTION_RLE_CELL_3D<T>* selection=(OPENGL_SELECTION_RLE_CELL_3D<T>*)current_selection;
        const RLE_RUN_3D* run=selection->run;int dj=selection->I.y-run->jmin,cell=run->cell+dj;
        int fx1=run->faces[0]+dj,fx2=run->faces[1]+dj,fy1=cell,fy2=fy1+(run->is_long?2:1),fz1=run->faces[4]+dj,fz2=run->faces[5]+dj;
        bool left_face_long=false,right_face_long=false,front_face_long=false,back_face_long=false;
        if(run->is_long && grid.long_run_faces_horizontal==2){ // TODO: make more efficient (maybe store some extra data in selection object to avoid this lookup)
            const RLE_RUN_3D* left_run=grid.Clamped_Run_In_Column(selection->I.x-1,selection->I.y,selection->I.z);
            if(left_run->is_long && (left_run+1)->jmin>selection->I.y+1) left_face_long=true;
            const RLE_RUN_3D* right_run=grid.Clamped_Run_In_Column(selection->I.x+1,selection->I.y,selection->I.z);
            if(right_run->is_long && (right_run+1)->jmin>selection->I.y+1) right_face_long=true;
            const RLE_RUN_3D* front_run=grid.Clamped_Run_In_Column(selection->I.x,selection->I.y,selection->I.z-1);
            if(front_run->is_long && (front_run+1)->jmin>selection->I.y+1) front_face_long=true;
            const RLE_RUN_3D* back_run=grid.Clamped_Run_In_Column(selection->I.x,selection->I.y,selection->I.z+1);
            if(back_run->is_long && (back_run+1)->jmin>selection->I.y+1) back_face_long=true;}
        output_stream<<"Selected cell "<<cell<<", I = "<<selection->I<<", center = "<<typename RLE_GRID_3D<T>::CELL_ITERATOR(grid,selection->I).Center()<<std::endl;
        output_stream<<"run sign = "<<run->is_long<<", index = "<<run-&grid.columns(selection->I.x,selection->I.z)(1)<<", jmin = "<<run->jmin<<", jmax = "<<(run+1)->jmin<<std::endl;
        output_stream<<"faces x = "<<fx1<<(left_face_long?"(long)":"")<<", "<<fx2<<(right_face_long?"(long)":"")<<", y = "<<fy1<<", "<<fy2
                     <<", z = "<<fz1<<(front_face_long?"(long)":"")<<", "<<fz2<<(back_face_long?"(long)":"")<<std::endl;}
}
//#####################################################################
// Function Create_Or_Destroy_Selection_After_Frame_Change
//#####################################################################
template<class T> OPENGL_SELECTION* OPENGL_RLE_GRID_3D<T>::
Create_Or_Destroy_Selection_After_Frame_Change(OPENGL_SELECTION* old_selection,bool& delete_selection)
{
    if(old_selection && old_selection->object==this){
        OPENGL_SELECTION_RLE_CELL_3D<T>* selection=(OPENGL_SELECTION_RLE_CELL_3D<T>*)current_selection;
        selection->run=grid.Clamped_Run_In_Column(selection->I);
        if(selection->run->is_long) selection->I.y=selection->run->jmin;
        Highlight_Selection(selection);
    }
    return 0;
}
//#####################################################################
template class OPENGL_RLE_GRID_3D<float>;
#ifndef COMPILE_WITHOUT_DOUBLE_SUPPORT
template class OPENGL_RLE_GRID_3D<double>;
#endif
#endif
