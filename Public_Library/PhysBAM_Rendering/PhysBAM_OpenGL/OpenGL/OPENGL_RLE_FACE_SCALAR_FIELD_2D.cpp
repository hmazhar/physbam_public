#ifndef COMPILE_WITHOUT_RLE_SUPPORT
//#####################################################################
// Copyright 2005, Geoffrey Irving.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#include <PhysBAM_Tools/Grids_RLE/RLE_GRID_ITERATOR_FACE_HORIZONTAL.h>
#include <PhysBAM_Tools/Grids_RLE/RLE_GRID_ITERATOR_FACE_Y.h>
#include <PhysBAM_Rendering/PhysBAM_OpenGL/OpenGL/OPENGL_RLE_FACE_SCALAR_FIELD_2D.h>
#include <PhysBAM_Rendering/PhysBAM_OpenGL/OpenGL/OPENGL_RLE_GRID_2D.h>
using namespace PhysBAM;
//#####################################################################
// Function Display
//#####################################################################
template<class T,class T2> void OPENGL_RLE_FACE_SCALAR_FIELD_2D<T,T2>::
Display(const int in_color) const
{
    PHYSBAM_ASSERT(grid.long_run_cells==2);
    glPushAttrib(GL_LIGHTING_BIT | GL_TEXTURE_BIT | GL_POINT_BIT);
    glPointSize(point_size);glDisable(GL_LIGHTING);glDisable(GL_TEXTURE_2D);

    if(draw_points){
        OpenGL_Begin(GL_POINTS);
        if(grid.long_run_faces_horizontal==1)
            for(FACE_X_ITERATOR face(grid,grid.number_of_ghost_cells);face;face++){int f=face.Face();
                if(value(f)){color_map->Lookup(value(f)).Send_To_GL_Pipeline();OpenGL_Vertex(face.X());}}
        else
            for(FACE_X_ITERATOR face(grid,grid.number_of_ghost_cells);face;face++){int f=face.Face();VECTOR<T,2> X=face.X();
                if(value(f)){color_map->Lookup(value(f)).Send_To_GL_Pipeline();X.y=grid.uniform_grid.Axis_X_plus_half(face.j(),2);OpenGL_Vertex(X);}
                if(face.Long() && value(f+1)){color_map->Lookup(value(f+1)).Send_To_GL_Pipeline();X.y=grid.uniform_grid.Axis_X_plus_half(face.jmax()-1,2);OpenGL_Vertex(X);}}
        for(FACE_Y_ITERATOR face(grid,grid.number_of_ghost_cells,true);face;face++){int f=face.Face();
            if(value(f)){color_map->Lookup(value(f)).Send_To_GL_Pipeline();OpenGL_Vertex(face.X());}
            if(face.Long() && value(f+1)){color_map->Lookup(value(f+1)).Send_To_GL_Pipeline();
                OpenGL_Vertex(face.cell2.Center());}}
        OpenGL_End();}
    else{
        OpenGL_Begin(GL_LINES);
        if(grid.long_run_faces_horizontal==1)
            for(FACE_X_ITERATOR face(grid,grid.number_of_ghost_cells);face;face++){int f=face.Face();VECTOR<T,2> X=face.X(),Y=X;
                if(value(f)){color_map->Lookup(value(f)).Send_To_GL_Pipeline();X.x+=value(f)*line_size;OpenGL_Vertex(Y,X);}}
        else
            for(FACE_X_ITERATOR face(grid,grid.number_of_ghost_cells);face;face++){int f=face.Face();VECTOR<T,2> X=face.X();
                if(value(f)){
                    color_map->Lookup(value(f)).Send_To_GL_Pipeline();X.y=grid.uniform_grid.Axis_X_plus_half(face.j(),2);
                    OpenGL_Line(X,VECTOR<T,2>(X.x+value(f)*line_size,X.y));}
                if(face.Long() && value(f+1)){
                    color_map->Lookup(value(f+1)).Send_To_GL_Pipeline();X.y=grid.uniform_grid.Axis_X_plus_half(face.jmax()-1,2);
                    OpenGL_Line(X,VECTOR<T,2>(X.x+value(f+1)*line_size,X.y));}}
        for(FACE_Y_ITERATOR face(grid,grid.number_of_ghost_cells,true);face;face++){int f=face.Face();
            if(value(f)){color_map->Lookup(value(f)).Send_To_GL_Pipeline();VECTOR<T,2> X=face.X();OpenGL_Vertex(X);X.y+=value(f)*line_size;OpenGL_Vertex(X);}
            if(face.Long() && value(f+1)){color_map->Lookup(value(f+1)).Send_To_GL_Pipeline();
                VECTOR<T,2> X=face.cell2.Center();
                OpenGL_Vertex(X);X.y+=value(f+1)*line_size;OpenGL_Vertex(X);}}
        OpenGL_End();}
    glPopAttrib();
}
//#####################################################################
// Function Bounding_Box
//#####################################################################
template<class T,class T2> RANGE<VECTOR<float,3> > OPENGL_RLE_FACE_SCALAR_FIELD_2D<T,T2>::
Bounding_Box() const
{
    return RANGE<VECTOR<float,3> >(grid.uniform_grid.domain.min_corner.x,grid.uniform_grid.domain.max_corner.x,grid.uniform_grid.domain.min_corner.y,grid.uniform_grid.domain.max_corner.y,0,0);
}
//#####################################################################
// Function Print_Selection_Info
//#####################################################################
template<class T,class T2> void OPENGL_RLE_FACE_SCALAR_FIELD_2D<T,T2>::
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

        ARRAY<T2>& V=value;
        output_stream<<"u left = "<<V(fx1);if(left_face_long) output_stream<<" "<<V(fx1+1);output_stream<<", ";
        output_stream<<"right = "<<V(fx2);if(right_face_long) output_stream<<" "<<V(fx2+1);output_stream<<std::endl;
        output_stream<<"v bottom = "<<V(fy1);if(run->is_long) output_stream<<", mid = "<<V(fy1+1);output_stream<<", top = "<<V(fy2)<<std::endl;
        if(!run->is_long){
            T ux=(V(fx2)-V(fx1))*grid.uniform_grid.one_over_dX.x,vy=(V(fy2)-V(fy1))*grid.uniform_grid.one_over_dX.y;
            output_stream<<"divergence = "<<ux+vy<<" (ux="<<ux<<", vy="<<vy<<")"<<std::endl;}}
}
//#####################################################################
template class OPENGL_RLE_FACE_SCALAR_FIELD_2D<float>;
template class OPENGL_RLE_FACE_SCALAR_FIELD_2D<float,int>;
template class OPENGL_RLE_FACE_SCALAR_FIELD_2D<float,bool>;
#ifndef COMPILE_WITHOUT_DOUBLE_SUPPORT
template class OPENGL_RLE_FACE_SCALAR_FIELD_2D<double>;
template class OPENGL_RLE_FACE_SCALAR_FIELD_2D<double,int>;
template class OPENGL_RLE_FACE_SCALAR_FIELD_2D<double,bool>;
#endif
#endif
