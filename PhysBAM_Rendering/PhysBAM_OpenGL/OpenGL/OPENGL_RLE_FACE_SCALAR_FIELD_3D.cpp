#ifndef COMPILE_WITHOUT_RLE_SUPPORT
//#####################################################################
// Copyright 2005, Eran Guendelman, Geoffrey Irving.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#include <PhysBAM_Tools/Grids_RLE/RLE_GRID_ITERATOR_FACE_HORIZONTAL.h>
#include <PhysBAM_Tools/Grids_RLE/RLE_GRID_ITERATOR_FACE_Y.h>
#include <PhysBAM_Rendering/PhysBAM_OpenGL/OpenGL/OPENGL_RLE_FACE_SCALAR_FIELD_3D.h>
#include <PhysBAM_Rendering/PhysBAM_OpenGL/OpenGL/OPENGL_RLE_GRID_3D.h>
#include <PhysBAM_Rendering/PhysBAM_OpenGL/OpenGL/OPENGL_RLE_SLICE.h>
using namespace PhysBAM;
//#####################################################################
// Function Display
//#####################################################################
template<class T,class T2> void Compute_Maximum_Norm(const OPENGL_RLE_FACE_SCALAR_FIELD_3D<T,T2>& field,const RANGE<VECTOR<int,2> >& region)
{
    field.max_norm=0;field.max_norm_cell=0;field.max_norm_I=VECTOR<int,3>();
    for(typename RLE_GRID_3D<T>::CELL_ITERATOR cell(field.grid,region);cell;cell++){
        T2 local_norm=0;
        for(int axis=1;axis<=RLE_GRID_3D<T>::dimension;axis++)local_norm+=maxabs(field.value(cell.First_Face_Index(axis)),field.value(cell.Second_Face_Index(axis)));
        if(field.max_norm<local_norm){field.max_norm=local_norm;field.max_norm_cell=cell.Cell();field.max_norm_I=cell.I();}}
}
template<class T> void Compute_Maximum_Norm(const OPENGL_RLE_FACE_SCALAR_FIELD_3D<T,bool>& field,const RANGE<VECTOR<int,2> >& region)
{}
template<class T,class T2> void OPENGL_RLE_FACE_SCALAR_FIELD_3D<T,T2>::
Display(const int in_color) const
{
    OPENGL_RLE_SLICE* slice=(OPENGL_RLE_SLICE*)this->slice;

    PHYSBAM_ASSERT(grid.long_run_cells==2);
    glPushAttrib(GL_LIGHTING_BIT | GL_TEXTURE_BIT | GL_POINT_BIT);
    glPointSize(point_size);glDisable(GL_LIGHTING);

    RANGE<VECTOR<int,2> > region(grid.columns.domain.min_corner.x+1,grid.columns.domain.max_corner.x-1,grid.columns.domain.min_corner.y+1,grid.columns.domain.max_corner.y-1);
    if(slice && slice->mode==OPENGL_SLICE::CELL_SLICE){
        switch(slice->axis){
            case 1: region.min_corner.x=region.max_corner.x=slice->index;break;
            case 2: break; // TODO: deal with y slices
            case 3: region.min_corner.y=region.max_corner.y=slice->index;break;}}

    Compute_Maximum_Norm(*this,region);
    RANGE<VECTOR<int,2> > xregion(region.min_corner.x,region.max_corner.x+1,region.min_corner.y,region.max_corner.y),zregion(region.min_corner.x,region.max_corner.x,region.min_corner.y,region.max_corner.y+1);

    if(draw_points){
        OpenGL_Begin(GL_POINTS);
        if(grid.long_run_faces_horizontal==1){
            for(FACE_X_ITERATOR face(grid,xregion);face;face++){int f=face.Face();
                if(value(f)){color_map->Lookup(value(f)).Send_To_GL_Pipeline();OpenGL_Vertex(face.X());}}
            for(FACE_Z_ITERATOR face(grid,zregion);face;face++){int f=face.Face();
                if(value(f)){color_map->Lookup(value(f)).Send_To_GL_Pipeline();OpenGL_Vertex(face.X());}}}
        else{
            for(FACE_X_ITERATOR face(grid,xregion);face;face++){int f=face.Face();VECTOR<T,3> X=face.X();
                if(value(f)){color_map->Lookup(value(f)).Send_To_GL_Pipeline();X.y=grid.uniform_grid.Axis_X_plus_half(face.j(),2);OpenGL_Vertex(X);}
                if(face.Long() && value(f+1)){color_map->Lookup(value(f+1)).Send_To_GL_Pipeline();X.y=grid.uniform_grid.Axis_X_plus_half(face.jmax()-1,2);OpenGL_Vertex(X);}}
            for(FACE_Z_ITERATOR face(grid,zregion);face;face++){int f=face.Face();VECTOR<T,3> X=face.X();
                if(value(f)){color_map->Lookup(value(f)).Send_To_GL_Pipeline();X.y=grid.uniform_grid.Axis_X_plus_half(face.j(),2);OpenGL_Vertex(face.X());}
                if(face.Long() && value(f+1)){color_map->Lookup(value(f+1)).Send_To_GL_Pipeline();X.y=grid.uniform_grid.Axis_X_plus_half(face.jmax()-1,2);OpenGL_Vertex(X);}}}
        for(FACE_Y_ITERATOR face(grid,region,true);face;face++){int f=face.Face();
            if(value(f)){color_map->Lookup(value(f)).Send_To_GL_Pipeline();OpenGL_Vertex(face.X());}
            if(face.Long() && value(f+1)){color_map->Lookup(value(f+1)).Send_To_GL_Pipeline();
                OpenGL_Vertex(face.cell2.Center());}}
        OpenGL_End();}
    else{
        OpenGL_Begin(GL_LINES);
        if(grid.long_run_faces_horizontal==1){
            for(FACE_X_ITERATOR face(grid,xregion);face;face++){int f=face.Face();
                if(value(f)){color_map->Lookup(value(f)).Send_To_GL_Pipeline();VECTOR<T,3> X=face.X(),Y=X;X.x+=value(f)*line_size;OpenGL_Vertex(Y,X);}}
            for(FACE_Z_ITERATOR face(grid,zregion);face;face++){int f=face.Face();
                if(value(f)){color_map->Lookup(value(f)).Send_To_GL_Pipeline();VECTOR<T,3> X=face.X(),Y=X;X.z+=value(f)*line_size;OpenGL_Vertex(Y,X);}}}
        else{
            for(FACE_X_ITERATOR face(grid,xregion);face;face++){int f=face.Face();VECTOR<T,3> X=face.X();
                if(value(f)){
                    color_map->Lookup(value(f)).Send_To_GL_Pipeline();X.y=grid.uniform_grid.Axis_X_plus_half(face.j(),2);
                    OpenGL_Line(X,VECTOR<T,3>(X.x+value(f)*line_size,X.y,X.z));}
                if(face.Long() && value(f+1)){
                    color_map->Lookup(value(f+1)).Send_To_GL_Pipeline();X.y=grid.uniform_grid.Axis_X_plus_half(face.jmax()-1,2);
                    OpenGL_Line(X,VECTOR<T,3>(X.x+value(f+1)*line_size,X.y,X.z));}}
            for(FACE_Z_ITERATOR face(grid,zregion);face;face++){int f=face.Face();VECTOR<T,3> X=face.X();
                if(value(f)){
                    color_map->Lookup(value(f)).Send_To_GL_Pipeline();X.y=grid.uniform_grid.Axis_X_plus_half(face.j(),2);
                    OpenGL_Line(X,VECTOR<T,3>(X.x,X.y,X.z+value(f)*line_size));}
                if(face.Long() && value(f+1)){
                    color_map->Lookup(value(f+1)).Send_To_GL_Pipeline();X.y=grid.uniform_grid.Axis_X_plus_half(face.jmax()-1,2);
                    OpenGL_Line(X,VECTOR<T,3>(X.x,X.y,X.z+value(f+1)*line_size));}}}
        for(FACE_Y_ITERATOR face(grid,region,true);face;face++){int f=face.Face();
            if(value(f)){color_map->Lookup(value(f)).Send_To_GL_Pipeline();VECTOR<T,3> X=face.X();OpenGL_Vertex(X);X.y+=value(f)*line_size;OpenGL_Vertex(X);}
            if(face.Long() && value(f+1)){color_map->Lookup(value(f+1)).Send_To_GL_Pipeline();
                VECTOR<T,3> X=face.cell2.Center();
                OpenGL_Vertex(X);X.y+=value(f+1)*line_size;OpenGL_Vertex(X);}}
        OpenGL_End();}
    glPopAttrib();
}
//#####################################################################
// Function Bounding_Box
//#####################################################################
template<class T,class T2> RANGE<VECTOR<float,3> > OPENGL_RLE_FACE_SCALAR_FIELD_3D<T,T2>::
Bounding_Box() const
{
    return (RANGE<VECTOR<float,3> >)grid.uniform_grid.domain;
}
//#####################################################################
// Function Print_Selection_Info
//#####################################################################
template<class T,class T2> void OPENGL_RLE_FACE_SCALAR_FIELD_3D<T,T2>::
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
        ARRAY<T2>& V=value;
        output_stream<<"u left = "<<V(fx1);if(left_face_long) output_stream<<" "<<V(fx1+1);output_stream<<", ";
        output_stream<<"right = "<<V(fx2);if(right_face_long) output_stream<<" "<<V(fx2+1);output_stream<<std::endl;
        output_stream<<"v bottom = "<<V(fy1);if(run->is_long) output_stream<<", mid = "<<V(fy1+1);output_stream<<", top = "<<V(fy2)<<std::endl;
        output_stream<<"w front = "<<V(fz1);if(front_face_long) output_stream<<" "<<V(fz1+1);output_stream<<", ";
        output_stream<<"back = "<<V(fz2);if(back_face_long) output_stream<<" "<<V(fz2+1);output_stream<<std::endl;
        if(!run->is_long){
            T ux=(V(fx2)-V(fx1))*grid.uniform_grid.one_over_dX.x,vy=(V(fy2)-V(fy1))*grid.uniform_grid.one_over_dX.y,wz=(V(fz2)-V(fz1))*grid.uniform_grid.one_over_dX.z;
            output_stream<<"divergence = "<<ux+vy+wz<<" (ux="<<ux<<", vy="<<vy<<", wz="<<wz<<")"<<std::endl;}}    
}
//#####################################################################
template class OPENGL_RLE_FACE_SCALAR_FIELD_3D<float>;
template class OPENGL_RLE_FACE_SCALAR_FIELD_3D<float,int>;
template class OPENGL_RLE_FACE_SCALAR_FIELD_3D<float,bool>;
#ifndef COMPILE_WITHOUT_DOUBLE_SUPPORT
template class OPENGL_RLE_FACE_SCALAR_FIELD_3D<double>;
template class OPENGL_RLE_FACE_SCALAR_FIELD_3D<double,int>;
template class OPENGL_RLE_FACE_SCALAR_FIELD_3D<double,bool>;
#endif
#endif
