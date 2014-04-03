//#####################################################################
// Copyright 2003-2009, Eran Guendelman, Michael Lentine, Andrew Selle, Jonathan Su.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#include <PhysBAM_Rendering/PhysBAM_OpenGL/OpenGL/OPENGL_GRID_3D.h>
#include <PhysBAM_Rendering/PhysBAM_OpenGL/OpenGL/OPENGL_PREFERENCES.h>
#include <PhysBAM_Rendering/PhysBAM_OpenGL/OpenGL/OPENGL_PRIMITIVES.h>
#include <PhysBAM_Rendering/PhysBAM_OpenGL/OpenGL/OPENGL_UNIFORM_SLICE.h>
using namespace PhysBAM;
//#####################################################################
// OPENGL_GRID_3D
//#####################################################################
template<class T> OPENGL_GRID_3D<T>::
OPENGL_GRID_3D(GRID<TV> &grid_input,const OPENGL_COLOR &color_input) 
    :grid(grid_input),color(color_input),draw_ghost_values(true),hide_non_selected_grid(false),scale(1),current_selection(0)
{
}
//#####################################################################
// Display
//#####################################################################
template<class T> void OPENGL_GRID_3D<T>::
Display(const int in_color) const
{
    glMatrixMode(GL_MODELVIEW);
    glPushMatrix();
    Send_Transform_To_GL_Pipeline();

    glPushAttrib(GL_ENABLE_BIT);
    glDisable(GL_LIGHTING);

    color.Send_To_GL_Pipeline();

    GLint mode;
    glGetIntegerv(GL_RENDER_MODE, &mode);

    OPENGL_UNIFORM_SLICE* slice=(OPENGL_UNIFORM_SLICE*)this->slice;

    int ghost_cells=draw_ghost_values?3:0;

    if (!slice || slice->mode==OPENGL_SLICE::NO_SLICE)
    {
        VECTOR<int,3> node_start(1,1,1),node_end(grid.numbers_of_cells+1);
        if(!hide_non_selected_grid) Draw_Subgrid(node_start,node_end);
    }
    else if (slice->mode == OPENGL_SLICE::NODE_SLICE)
    {
        VECTOR<int,3> node_start(1,1,1),node_end(grid.numbers_of_cells+1);
        if (slice->axis == 1) { node_start.x = node_end.x = (slice->index-1)/scale; }
        else if (slice->axis == 2) { node_start.y = node_end.y = (slice->index-1)/scale+1; }
        else if (slice->axis == 3) { node_start.z = node_end.z = (slice->index-1)/scale+1; }

        if (mode == GL_SELECT)
        {
            // Currently only support node selection in this mode
            glPushName(2);
            Draw_Nodes_For_Selection(node_start,node_end);
            glPopName();
        }
        else
            Draw_Subgrid(node_start,node_end);
    }
    else if (slice->mode == OPENGL_SLICE::CELL_SLICE)
    {
        if (mode == GL_SELECT)
        {
            // Currently only support cell selection in this mode
            glPushAttrib(GL_ENABLE_BIT);
            glDisable(GL_CULL_FACE);

            VECTOR<T,3> x_vector(grid.dX.x,0,0),y_vector(0,grid.dX.y,0),z_vector(0,0,grid.dX.z);

            T x, y, z;
            int i, j, k;
            int i_start, i_end, j_start, j_end, k_start, k_end;
            VECTOR<T,3> axis_1, axis_2, axis_3;
            if (slice->axis == 1) { i_start = i_end = (slice->index-1)/scale+1; axis_1 = y_vector; axis_2 = z_vector; axis_3 = x_vector; } 
            else { i_start = 1-ghost_cells; i_end = grid.numbers_of_cells.x+ghost_cells; }
            if (slice->axis == 2) { j_start = j_end = (slice->index-1)/scale+1; axis_1 = z_vector; axis_2 = x_vector; axis_3 = y_vector; } 
            else { j_start = 1-ghost_cells; j_end = grid.numbers_of_cells.y+ghost_cells; }
            if (slice->axis == 3) { k_start = k_end = (slice->index-1)/scale+1; axis_1 = x_vector; axis_2 = y_vector; axis_3 = z_vector; } 
            else { k_start = 1-ghost_cells; k_end = grid.numbers_of_cells.z+ghost_cells; }

            VECTOR<T,3> pos_start=grid.Node(i_start,j_start,k_start);

            glPushName(1);
            for (i=i_start, x=pos_start.x; i<=i_end; i++, x+=grid.dX.x)
            {
                glPushName(i);
                for (j=j_start, y=pos_start.y; j<=j_end; j++, y+=grid.dX.y)
                {
                    glPushName(j);
                    for (k=k_start, z=pos_start.z; k<=k_end; k++, z+=grid.dX.z)
                    {
                        VECTOR<T,3> min_corner(x,y,z);
                        glPushName(k);
                        OpenGL_Begin(GL_QUADS);
                        OpenGL_Quad(min_corner,axis_1,axis_2);
                        OpenGL_Quad(min_corner+axis_3,axis_1,axis_2);
                        OpenGL_End();
                        glPopName();
                    }
                    glPopName();
                }
                glPopName();
            }

            glPopName();
            glPopAttrib();
        }
        else
        {
            VECTOR<int,3> node_start(1-ghost_cells,1-ghost_cells,1-ghost_cells),node_end(grid.numbers_of_cells+ghost_cells+1);
            if (slice->axis == 1) { node_start.x = (slice->index-1)/scale+1; node_end.x = node_start.x+1; }
            else if (slice->axis == 2) { node_start.y = (slice->index-1)/scale+1; node_end.y = node_start.y+1; }
            else if (slice->axis == 3) { node_start.z = (slice->index-1)/scale+1; node_end.z = node_start.z+1; }
            Draw_Subgrid(node_start,node_end);

            // Outline boundary of real domain in wider line
            if(ghost_cells>0){
                VECTOR<T,3> x000,x111;
                if(slice->axis==1){x000=VECTOR<T,3>(grid.domain.min_corner.x+grid.dX.x*((slice->index-1)/scale+1),grid.domain.min_corner.y,grid.domain.min_corner.z);x111=VECTOR<T,3>(grid.domain.min_corner.x+grid.dX.x*(slice->index-1)/scale,grid.domain.max_corner.y,grid.domain.max_corner.z);}
                else if(slice->axis==2){x000=VECTOR<T,3>(grid.domain.min_corner.x,grid.domain.min_corner.y+grid.dX.y*((slice->index-1)/scale+1),grid.domain.min_corner.z);x111=VECTOR<T,3>(grid.domain.max_corner.x,grid.domain.min_corner.y+grid.dX.y*(slice->index-1)/scale,grid.domain.max_corner.z);}
                else if(slice->axis==3){x000=VECTOR<T,3>(grid.domain.min_corner.x,grid.domain.min_corner.y,grid.domain.min_corner.z+grid.dX.z*((slice->index-1)/scale+1));x111=VECTOR<T,3>(grid.domain.max_corner.x,grid.domain.max_corner.y,grid.domain.min_corner.z+grid.dX.z*(slice->index-1)/scale);}
                glPushAttrib(GL_LINE_BIT);glLineWidth(3*OPENGL_PREFERENCES::line_width);
                VECTOR<T,3> x001(x000.x,x000.y,x111.z),x010(x000.x,x111.y,x000.z),x011(x000.x,x111.y,x111.z),
                    x100(x111.x,x000.y,x000.z),x101(x111.x,x000.y,x111.z),x110(x111.x,x111.y,x000.z);
                
                OpenGL_Begin(GL_LINES);
                OpenGL_Line(x000,x010);
                OpenGL_Line(x010,x011);
                OpenGL_Line(x011,x001);
                OpenGL_Line(x001,x000);

                OpenGL_Line(x100,x110);
                OpenGL_Line(x110,x111);
                OpenGL_Line(x111,x101);
                OpenGL_Line(x101,x100);
                
                OpenGL_Line(x000,x100);
                OpenGL_Line(x010,x110);
                OpenGL_Line(x011,x111);
                OpenGL_Line(x001,x101);
                OpenGL_End();
                glPopAttrib();}
        }
    }

    if (current_selection)
    {
        if (current_selection->type == OPENGL_SELECTION::GRID_CELL_3D)
        {
            OPENGL_SELECTION_GRID_CELL_3D<T> *real_selection = (OPENGL_SELECTION_GRID_CELL_3D<T>*)current_selection;
            VECTOR<T,3> min_corner=grid.Node(real_selection->index),max_corner=min_corner+grid.dX;
            OPENGL_SELECTION::Draw_Highlighted_Box(min_corner,max_corner);
        }
        else if (current_selection->type == OPENGL_SELECTION::GRID_NODE_3D)
        {
            OPENGL_SELECTION_GRID_NODE_3D<T> *real_selection = (OPENGL_SELECTION_GRID_NODE_3D<T>*)current_selection;
            OPENGL_SELECTION::Draw_Highlighted_Vertex(grid.Node(real_selection->index));
        }
        else if (current_selection->type == OPENGL_SELECTION::GRID_CELL_LIST_3D)
        {
            OPENGL_SELECTION_GRID_CELL_LIST_3D<T> *real_selection = (OPENGL_SELECTION_GRID_CELL_LIST_3D<T>*)current_selection;
            for(int i=1;i<=real_selection->indicies.m;i++){
                VECTOR<T,3> min_corner=grid.Node(real_selection->indicies(i)),max_corner=min_corner+grid.dX;
                OPENGL_SELECTION::Draw_Highlighted_Box(min_corner,max_corner);}
        }
        else if (current_selection->type == OPENGL_SELECTION::GRID_NODE_LIST_3D)
        {
            OPENGL_SELECTION_GRID_NODE_LIST_3D<T> *real_selection = (OPENGL_SELECTION_GRID_NODE_LIST_3D<T>*)current_selection;
            for(int i=1;i<=real_selection->indicies.m;i++) OPENGL_SELECTION::Draw_Highlighted_Vertex(grid.Node(real_selection->indicies(i)));
        }
    }

    glPopAttrib();
    glPopMatrix();
}
//#####################################################################
// Draw_Subgrid
//#####################################################################
template<class T> void OPENGL_GRID_3D<T>::
Draw_Subgrid(const VECTOR<int,3> &node_start,const VECTOR<int,3> &node_end) const
{
    int i,j,k;
    T x,y,z;

    VECTOR<T,3> start_position=grid.Node(node_start),end_position=grid.Node(node_end);

    OpenGL_Begin(GL_LINES);

    if(node_start.z!=node_end.z)
        for (i=node_start.x, x=start_position.x; i<=node_end.x; i++, x+=grid.dX.x)
            for (j=node_start.y, y=start_position.y; j<=node_end.y; j++, y+=grid.dX.y)
            {
                OpenGL_Line(VECTOR<T,3>(x,y,start_position.z),VECTOR<T,3>(x,y,end_position.z));
            }

    if(node_start.y!=node_end.y)
        for (i=node_start.x, x=start_position.x; i<=node_end.x; i++, x+=grid.dX.x)
            for (k=node_start.z, z=start_position.z; k<=node_end.z; k++, z+=grid.dX.z)
            {
                OpenGL_Line(VECTOR<T,3>(x,start_position.y,z),VECTOR<T,3>(x,end_position.y,z));
            }

    if(node_start.x!=node_end.x)
        for (j=node_start.y, y=start_position.y; j<=node_end.y; j++, y+=grid.dX.y)
            for (k=node_start.z, z=start_position.z; k<=node_end.z; k++, z+=grid.dX.z)
            {
                OpenGL_Line(VECTOR<T,3>(start_position.x,y,z),VECTOR<T,3>(end_position.x,y,z));
            }

    OpenGL_End();
}
//#####################################################################
// Draw_Nodes_For_Selection
//#####################################################################
template<class T> void OPENGL_GRID_3D<T>::
Draw_Nodes_For_Selection(const VECTOR<int,3> &node_start,const VECTOR<int,3> &node_end) const
{
    int i,j,k;
    T x,y,z;

    VECTOR<T,3> start_position=grid.Node(node_start);

    glPushAttrib(GL_POINT_BIT);
    glPointSize(OPENGL_PREFERENCES::selection_point_size);

    for(i=node_start.x, x=start_position.x; i<=node_end.x; i++, x+=grid.dX.x)
    {
        glPushName(i);
        for(j=node_start.y, y=start_position.y; j<=node_end.y; j++, y+=grid.dX.y)
        {
            glPushName(j);
            for(k=node_start.z, z=start_position.z; k<=node_end.z; k++, z+=grid.dX.z)
            {
                glPushName(k);
                OpenGL_Begin(GL_POINTS);
                OpenGL_Vertex(VECTOR<T,3>(x,y,z));
                OpenGL_End();
                glPopName();
            }
            glPopName();
        }
        glPopName();
    }

    glPopAttrib();
}
//#####################################################################
// Function Bounding_Box
//#####################################################################
template<class T> RANGE<VECTOR<float,3> > OPENGL_GRID_3D<T>::
Bounding_Box() const
{
    return World_Space_Box(grid.domain);
}
//#####################################################################
// Function Get_Selection
//#####################################################################
template<class T> OPENGL_SELECTION *OPENGL_GRID_3D<T>::
Get_Selection(GLuint *buffer,int buffer_size)
{
    OPENGL_SELECTION *selection = 0;
    if (buffer_size == 4)
    {
        if(buffer[0] == 1)
            selection = new OPENGL_SELECTION_GRID_CELL_3D<T>(this,VECTOR<int,3>(buffer[1],buffer[2],buffer[3]));
        else if (buffer[0] == 2)
            selection = new OPENGL_SELECTION_GRID_NODE_3D<T>(this,VECTOR<int,3>(buffer[1],buffer[2],buffer[3]));
    } 
    if (buffer_size%3 == 1)
    {
        ARRAY<VECTOR<int,3> > indicies;
        for(int i=0;i<buffer_size/3;i++) indicies.Append(VECTOR<int,3>(buffer[3*i+1],buffer[3*i+2],buffer[3*i+3]));
        if (buffer[0] == 3) selection = new OPENGL_SELECTION_GRID_CELL_LIST_3D<T>(this,indicies);
        if (buffer[0] == 4) selection = new OPENGL_SELECTION_GRID_NODE_LIST_3D<T>(this,indicies);
    }

    return selection;
}
//#####################################################################
// Function Highlight_Selection
//#####################################################################
template<class T> void OPENGL_GRID_3D<T>::
Highlight_Selection(OPENGL_SELECTION *selection)
{
    delete current_selection; current_selection = 0;
    if (selection->type == OPENGL_SELECTION::GRID_CELL_3D)
    {
        OPENGL_SELECTION_GRID_CELL_3D<T> *real_selection = (OPENGL_SELECTION_GRID_CELL_3D<T>*)selection;
        current_selection = new OPENGL_SELECTION_GRID_CELL_3D<T>(this,real_selection->index);
    }
    else if (selection->type == OPENGL_SELECTION::GRID_NODE_3D)
    {
        OPENGL_SELECTION_GRID_NODE_3D<T> *real_selection = (OPENGL_SELECTION_GRID_NODE_3D<T>*)selection;
        current_selection = new OPENGL_SELECTION_GRID_NODE_3D<T>(this,real_selection->index);
    }
    else if (selection->type == OPENGL_SELECTION::GRID_CELL_LIST_3D)
    {
        OPENGL_SELECTION_GRID_CELL_LIST_3D<T> *real_selection = (OPENGL_SELECTION_GRID_CELL_LIST_3D<T>*)selection;
        current_selection = new OPENGL_SELECTION_GRID_CELL_LIST_3D<T>(this,real_selection->indicies);
    }
    else if (selection->type == OPENGL_SELECTION::GRID_NODE_LIST_3D)
    {
        OPENGL_SELECTION_GRID_NODE_LIST_3D<T> *real_selection = (OPENGL_SELECTION_GRID_NODE_LIST_3D<T>*)selection;
        current_selection = new OPENGL_SELECTION_GRID_NODE_LIST_3D<T>(this,real_selection->indicies);
    }
}
//#####################################################################
// Function Clear_Highlight
//#####################################################################
template<class T> void OPENGL_GRID_3D<T>::
Clear_Highlight()
{
    delete current_selection; current_selection = 0;
}
//#####################################################################
// Function Bounding_Box
//#####################################################################
template<class T> RANGE<VECTOR<float,3> > OPENGL_SELECTION_GRID_CELL_3D<T>::
Bounding_Box() const
{
    PHYSBAM_ASSERT(object);
    const GRID<TV> &grid=((OPENGL_GRID_3D<T> *)object)->grid;
    VECTOR<T,3> min_corner=grid.Node(index),max_corner=min_corner+grid.dX;
    RANGE<VECTOR<float,3> > box((VECTOR<float,3>)min_corner,(VECTOR<float,3>)max_corner);
    return object->World_Space_Box(box);
}
//#####################################################################
// Function Bounding_Box
//#####################################################################
template<class T> RANGE<VECTOR<float,3> > OPENGL_SELECTION_GRID_NODE_3D<T>::
Bounding_Box() const
{
    PHYSBAM_ASSERT(object);
    const GRID<TV> &grid=((OPENGL_GRID_3D<T> *)object)->grid;
    RANGE<VECTOR<float,3> > box(VECTOR<float,3>(grid.Node(index)));
    return object->World_Space_Box(box);
}
//#####################################################################
// Function Bounding_Box
//#####################################################################
template<class T> RANGE<VECTOR<float,3> > OPENGL_SELECTION_GRID_CELL_LIST_3D<T>::
Bounding_Box() const
{
    PHYSBAM_ASSERT(object);
    const GRID<TV> &grid=((OPENGL_GRID_3D<T> *)object)->grid;
    VECTOR<T,3> min_corner=grid.Node(indicies(1)),max_corner=min_corner+grid.dX;
    for(int i=2;i<=indicies.m;i++){
        min_corner=TV::Componentwise_Min(min_corner,grid.Node(indicies(i)));
        max_corner=TV::Componentwise_Max(max_corner,grid.Node(indicies(i))+grid.dX);}
    RANGE<VECTOR<float,3> > box((VECTOR<float,3>)min_corner,(VECTOR<float,3>)max_corner);
    return object->World_Space_Box(box);
}
//#####################################################################
// Function Print_Selection_Info
//#####################################################################
template<class T> void OPENGL_GRID_3D<T>::
Print_Selection_Info(std::ostream& stream,OPENGL_SELECTION* selection) const
{
    if(current_selection && current_selection->type==OPENGL_SELECTION::GRID_NODE_3D){
        VECTOR<int,3> index=((OPENGL_SELECTION_GRID_NODE_3D<T>*)current_selection)->index;
        stream<<"Selected node "<<index<<" ("<<grid.Get_Regular_Grid().X(index)<<")"<<std::endl;}
    else if(current_selection && current_selection->type==OPENGL_SELECTION::GRID_CELL_3D){
        VECTOR<int,3> index=((OPENGL_SELECTION_GRID_CELL_3D<T>*)current_selection)->index;
        stream<<"Selected cell "<<index<<" ("<<grid.Get_MAC_Grid().X(index)<<")"<<std::endl;}
}
//#####################################################################
// Function Bounding_Box
//#####################################################################
template<class T> RANGE<VECTOR<float,3> > OPENGL_SELECTION_GRID_NODE_LIST_3D<T>::
Bounding_Box() const
{
    PHYSBAM_ASSERT(object);
    const GRID<TV> &grid=((OPENGL_GRID_3D<T> *)object)->grid;
    VECTOR<T,3> min_corner=grid.Node(indicies(1)),max_corner=min_corner;
    for(int i=2;i<=indicies.m;i++){
        if(grid.Node(indicies(i)).x<min_corner.x) min_corner.x=grid.Node(indicies(i)).x;
        if(grid.Node(indicies(i)).y<min_corner.y) min_corner.y=grid.Node(indicies(i)).y;
        if(grid.Node(indicies(i)).z<min_corner.z) min_corner.z=grid.Node(indicies(i)).z;
        if(grid.Node(indicies(i)).x>max_corner.x) max_corner.x=grid.Node(indicies(i)).x;
        if(grid.Node(indicies(i)).y>max_corner.y) max_corner.y=grid.Node(indicies(i)).y;
        if(grid.Node(indicies(i)).z>max_corner.z) max_corner.z=grid.Node(indicies(i)).z;}
    RANGE<VECTOR<float,3> > box((VECTOR<float,3>)min_corner,(VECTOR<float,3>)max_corner);
    return object->World_Space_Box(box);
}
//#####################################################################
// Function Bounding_Box
//#####################################################################
template<class T> void OPENGL_GRID_3D<T>::
Toggle_Draw_Ghost_Values()
{
    draw_ghost_values=!draw_ghost_values;
}
template class OPENGL_GRID_3D<float>;
#ifndef COMPILE_WITHOUT_DOUBLE_SUPPORT
template class OPENGL_GRID_3D<double>;
#endif
