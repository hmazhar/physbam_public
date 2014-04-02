#ifndef COMPILE_WITHOUT_DYADIC_SUPPORT
//#####################################################################
// Copyright 2004-2009, Frank Losasso, Andrew Selle.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#include <PhysBAM_Tools/Grids_Dyadic/DYADIC_GRID_ITERATOR_FACE.h>
#include <PhysBAM_Geometry/Basic_Geometry/ORIENTED_BOX.h>
#include <PhysBAM_Rendering/PhysBAM_OpenGL/OpenGL/OPENGL_OCTREE_GRID.h>
#include <PhysBAM_Rendering/PhysBAM_OpenGL/OpenGL/OPENGL_OCTREE_SLICE.h>
#include <PhysBAM_Rendering/PhysBAM_OpenGL/OpenGL/OPENGL_PREFERENCES.h>
using namespace PhysBAM;


template<class T> void OPENGL_OCTREE_GRID<T>::
Display(const int in_color) const
{
    PHYSBAM_ASSERT(slice);OPENGL_OCTREE_SLICE* slice=(OPENGL_OCTREE_SLICE*)this->slice;
    glMatrixMode(GL_MODELVIEW);glPushMatrix();
    Send_Transform_To_GL_Pipeline();

    glPushAttrib(GL_ENABLE_BIT | GL_CURRENT_BIT);glDisable(GL_LIGHTING);
    color.Send_To_GL_Pipeline();
    GLint mode;glGetIntegerv(GL_RENDER_MODE, &mode);

    if(mode==GL_SELECT){ 
        if(slice && slice->Is_Slice_Mode()){
            glPushAttrib(GL_ENABLE_BIT | GL_POINT_BIT);
            glDisable(GL_CULL_FACE);

            // TODO Draw grid faces for selection
            OPENGL_OCTREE_SLICE* slice=(OPENGL_OCTREE_SLICE*)this->slice;
    
            glPushName(1);
            glPushName(0); // overwritten per cell in Draw_Cells_For_Selection_Helper
            for(DYADIC_GRID_ITERATOR_FACE<OCTREE_GRID<T> > iterator(grid,grid.Map_All_Faces());iterator.Valid();iterator.Next()){
                if(slice && iterator.Axis()+1!=slice->axis) continue;
                if(slice && abs(iterator.Location()[slice->axis]-slice->position)>.00001) continue;
                glLoadName(iterator.Face_Index());
                for(int k=0;k<4;k++) glPushName(iterator.Face_Node(k));
                glPushName(iterator.Deepest_Cell_Index());glPushName(iterator.Other_Cell_Index());
                OpenGL_Begin(GL_QUADS);
                for(int k=0;k<4;k++) OpenGL_Vertex(grid.Node_Location(iterator.Face_Node(k)));
                OpenGL_End();
                glPopName();glPopName();glPopName();glPopName();glPopName();glPopName();
            }

            glPopName();

            glLoadName(2);
            glPushName(0); // overwritten per node below
            glPointSize(OPENGL_PREFERENCES::selection_point_size);
            for(int i=1;i<=grid.number_of_nodes;i++){
                if(abs(grid.Node_Location(i)[slice->axis]-(T)slice->position)>.00001) continue;
                glLoadName(i);
                OpenGL_Begin(GL_POINTS);
                OpenGL_Vertex(grid.Node_Location(i));
                OpenGL_End();}
            glPopName();
            glPopName();
            glPopAttrib();}}
    else{
        if(!slice || !slice->Is_Slice_Mode()) Draw_Refined_Cells(); // don't use display list in non-slice mode because it's too slow
        else{
            if(reinitialize_display_list){reinitialize_display_list=false;
                glNewList(display_list_id,GL_COMPILE_AND_EXECUTE);
                Draw_Coarse_Cells_In_Slice_Mode();Draw_Refined_Cells();
                OpenGL_EndList();
            }
            else glCallList(display_list_id);}

        // Highlight current selection
        if(current_selection){
            if(current_selection->type==OPENGL_SELECTION::OCTREE_CELL){
                OPENGL_SELECTION_OCTREE_CELL<T>* real_selection=(OPENGL_SELECTION_OCTREE_CELL<T>*)current_selection;
                const PhysBAM::OCTREE_CELL<T> *cell=grid.Cell_Pointer_From_Index()(real_selection->index);
                RANGE<VECTOR<T,3> > box=cell->Bounding_Box();
                OPENGL_SELECTION::Draw_Highlighted_Box(box.Minimum_Corner(),box.Maximum_Corner());
            }
            else if(current_selection->type==OPENGL_SELECTION::OCTREE_FACE){
                OPENGL_SELECTION_OCTREE_FACE<T>* real_selection=(OPENGL_SELECTION_OCTREE_FACE<T>*)current_selection;
                OPENGL_SELECTION::Draw_Highlighted_Quad(grid.Node_Location(real_selection->node[0]),grid.Node_Location(real_selection->node[1]),grid.Node_Location(real_selection->node[2]),grid.Node_Location(real_selection->node[3]));
            }
            else if(current_selection->type==OPENGL_SELECTION::OCTREE_NODE){
                OPENGL_SELECTION_OCTREE_NODE<T>* real_selection=(OPENGL_SELECTION_OCTREE_NODE<T>*)current_selection;
                OPENGL_SELECTION::Draw_Highlighted_Vertex(real_selection->location);
            }
        }
    }
    glPopAttrib();
    glPopMatrix();
}

template<class T> void OPENGL_OCTREE_GRID<T>::
Draw_Coarse_Cells_In_Slice_Mode() const
{
    // draw the coarse uniform uniform_grid
    glPushAttrib(GL_DEPTH_BUFFER_BIT);
    glDepthMask(GL_FALSE);
    OpenGL_Begin(GL_LINES);
    OPENGL_OCTREE_SLICE* slice=(OPENGL_OCTREE_SLICE*)this->slice;
    VECTOR<int,3> index=grid.uniform_grid.Closest_Node(VECTOR<T,3>(slice->position,slice->position,slice->position));
    if(slice->axis==1){/*if(abs(grid.uniform_grid.z(index[slice->axis])-slice->position)<1e-5)*/{
        T y,z;int j,ij,n=grid.uniform_grid.numbers_of_cells.y,mn=grid.uniform_grid.numbers_of_cells.z;T x=slice->position;
        for (j=1,y=grid.uniform_grid.domain.min_corner.y;j<=n+1;j++,y+=grid.uniform_grid.dX.y){OpenGL_Line(VECTOR<T,3>(x,y,grid.uniform_grid.domain.min_corner.z),VECTOR<T,3>(x,y,grid.uniform_grid.domain.max_corner.z));}
        for (ij=1,z=grid.uniform_grid.domain.min_corner.z;ij<=mn+1;ij++,z+=grid.uniform_grid.dX.z){OpenGL_Line(VECTOR<T,3>(x,grid.uniform_grid.domain.min_corner.y,z),VECTOR<T,3>(x,grid.uniform_grid.domain.max_corner.y,z));}}}
    else if(slice->axis==2){/*if(abs(grid.uniform_grid.z(index[slice->axis])-slice->position)<1e-5)*/{
        T x,z;int i,ij,m=grid.uniform_grid.numbers_of_cells.x,mn=grid.uniform_grid.numbers_of_cells.z;T y=slice->position;
        for (i=1,x=grid.uniform_grid.domain.min_corner.x;i<=m+1;i++,x+=grid.uniform_grid.dX.x){OpenGL_Line(VECTOR<T,3>(x,y,grid.uniform_grid.domain.min_corner.z),VECTOR<T,3>(x,y,grid.uniform_grid.domain.max_corner.z));}
        for (ij=1,z=grid.uniform_grid.domain.min_corner.z;ij<=mn+1;ij++,z+=grid.uniform_grid.dX.z){OpenGL_Line(VECTOR<T,3>(grid.uniform_grid.domain.min_corner.x,y,z),VECTOR<T,3>(grid.uniform_grid.domain.max_corner.x,y,z));}}}
    else if(slice->axis==3){/*if(abs(grid.uniform_grid.z(index[slice->axis])-slice->position)<1e-5)*/{
        T x,y;int i,j,m=grid.uniform_grid.numbers_of_cells.x,n=grid.uniform_grid.numbers_of_cells.y;T z=slice->position;
        for (i=1,x=grid.uniform_grid.domain.min_corner.x;i<=m+1;i++,x+=grid.uniform_grid.dX.x){OpenGL_Line(VECTOR<T,3>(x,grid.uniform_grid.domain.min_corner.y,z),VECTOR<T,3>(x,grid.uniform_grid.domain.max_corner.y,z));}
        for (j=1,y=grid.uniform_grid.domain.min_corner.y;j<=n+1;j++,y+=grid.uniform_grid.dX.y){OpenGL_Line(VECTOR<T,3>(grid.uniform_grid.domain.min_corner.x,y,z),VECTOR<T,3>(grid.uniform_grid.domain.max_corner.x,y,z));}}}
    OpenGL_End();
    glPopAttrib();
}

template<class T> void OPENGL_OCTREE_GRID<T>::
Draw_Refined_Cells() const
{
    glPushAttrib(GL_CURRENT_BIT | GL_ENABLE_BIT | GL_POLYGON_BIT);
    (color*.5).Send_To_GL_Pipeline();
    // draw the actual tree
    OpenGL_Begin(GL_LINES);

    OPENGL_OCTREE_SLICE* slice=(OPENGL_OCTREE_SLICE*)this->slice;
    for(DYADIC_GRID_ITERATOR_FACE<OCTREE_GRID<T> > iterator(grid,grid.Map_All_Faces());iterator.Valid();iterator.Next()){
        if(slice && slice->Is_Slice_Mode() && (iterator.Axis()+1!=slice->axis || abs(iterator.Location()[slice->axis]-(T)slice->position)>.00001)) continue;
        OpenGL_Line(grid.Node_Location(iterator.Face_Node(0)),grid.Node_Location(iterator.Face_Node(1)));
        OpenGL_Line(grid.Node_Location(iterator.Face_Node(0)),grid.Node_Location(iterator.Face_Node(2)));
        OpenGL_Line(grid.Node_Location(iterator.Face_Node(2)),grid.Node_Location(iterator.Face_Node(3)));
        OpenGL_Line(grid.Node_Location(iterator.Face_Node(1)),grid.Node_Location(iterator.Face_Node(3)));
    }

//    MAP_OCTREE_MESH::Map_Faces<T,const OPENGL_OCTREE_GRID<T>*>(grid.uniform_grid,grid.cells,grid.number_of_ghost_cells,this,Draw_Faces_Wireframe_Helper);    OPENGL_OCTREE_SLICE* slice=(OPENGL_OCTREE_SLICE*)this->slice;

    OpenGL_End(); 
    glPopAttrib();
}

template<class T> RANGE<VECTOR<float,3> > OPENGL_OCTREE_GRID<T>::
Bounding_Box() const
{
    return World_Space_Box(RANGE<VECTOR<float,3> >(grid.uniform_grid.domain));
}

template<class T> OPENGL_SELECTION *OPENGL_OCTREE_GRID<T>::
Get_Selection(GLuint *buffer,int buffer_size)
{
    OPENGL_SELECTION *selection=0;
    if(buffer_size==8 && buffer[0]==1){
        selection=new OPENGL_SELECTION_OCTREE_FACE<T>(this,buffer[1],buffer[2],buffer[3],buffer[4],buffer[5],buffer[6],buffer[7]);}
    else if(buffer_size==2 && buffer[0]==2) selection=new OPENGL_SELECTION_OCTREE_NODE<T>(this,buffer[1],grid.Node_Location(buffer[1]));
    return selection;
}

template<class T> void OPENGL_OCTREE_GRID<T>::
Highlight_Selection(OPENGL_SELECTION *selection)
{
    delete current_selection;current_selection=0;
    if (selection->type == OPENGL_SELECTION::OCTREE_CELL){
        OPENGL_SELECTION_OCTREE_CELL<T> *real_selection=(OPENGL_SELECTION_OCTREE_CELL<T>*)selection;
        current_selection=new OPENGL_SELECTION_OCTREE_CELL<T>(this,real_selection->index);}
    else if (selection->type==OPENGL_SELECTION::OCTREE_FACE){
        OPENGL_SELECTION_OCTREE_FACE<T> *real_selection=(OPENGL_SELECTION_OCTREE_FACE<T>*)selection;
        current_selection=new OPENGL_SELECTION_OCTREE_FACE<T>(this,real_selection->index,real_selection->node[0],real_selection->node[1],real_selection->node[2],real_selection->node[3],real_selection->cell[0],real_selection->cell[1]);}
    else if (selection->type == OPENGL_SELECTION::OCTREE_NODE){
        OPENGL_SELECTION_OCTREE_NODE<T> *real_selection=(OPENGL_SELECTION_OCTREE_NODE<T>*)selection;
        current_selection=new OPENGL_SELECTION_OCTREE_NODE<T>(this,real_selection->index,real_selection->location);}
}

template<class T> void OPENGL_OCTREE_GRID<T>::
Clear_Highlight()
{
    delete current_selection;current_selection=0;
}

template<class T> void OPENGL_OCTREE_GRID<T>::
Update()
{
    reinitialize_display_list=true;
}

template<class T> OPENGL_SELECTION* OPENGL_OCTREE_GRID<T>::
Get_Cell_Selection(const int cell_index)
{
    return new OPENGL_SELECTION_OCTREE_CELL<T>(this,cell_index);
}

template<class T> OPENGL_SELECTION* OPENGL_OCTREE_GRID<T>::
Get_Node_Selection(const int node_index)
{
    return new OPENGL_SELECTION_OCTREE_NODE<T>(this,node_index,grid.Node_Location(node_index));
}
//#####################################################################
// Create_Or_Destroy_Selection_After_Frame_Change
//#####################################################################
template<class T> OPENGL_SELECTION* OPENGL_OCTREE_GRID<T>::
Create_Or_Destroy_Selection_After_Frame_Change(OPENGL_SELECTION* old_selection,bool& delete_selection)
{
    if(old_selection && old_selection->object==this)
        return Get_Updated_Selection(old_selection);
    return 0;
}
template<class T> OPENGL_SELECTION* OPENGL_OCTREE_GRID<T>::
Get_Updated_Selection(OPENGL_SELECTION *selection)
{
    if (selection->type == OPENGL_SELECTION::OCTREE_NODE){
        OPENGL_SELECTION_OCTREE_NODE<T> *real_selection=(OPENGL_SELECTION_OCTREE_NODE<T>*)selection;
        OCTREE_CELL<T>* leaf=grid.Clamped_Leaf_Cell(real_selection->location,1e-5);
        T tolerance_squared=sqr(1e-5);
        int new_node_index=leaf->Nearest_Node(real_selection->location);
        VECTOR<T,3> new_node_location=grid.Node_Location(new_node_index);
        if((real_selection->location-new_node_location).Magnitude_Squared()<tolerance_squared) return new OPENGL_SELECTION_OCTREE_NODE<T>(this,new_node_index,new_node_location);}
    return 0;
}
//#####################################################################
// OPENGL_SELECTION_OCTREE_CELL
//#####################################################################
template<class T> RANGE<VECTOR<float,3> > OPENGL_SELECTION_OCTREE_CELL<T>::
Bounding_Box() const
{
    PHYSBAM_ASSERT(object);
    const OCTREE_GRID<T> &grid=((OPENGL_OCTREE_GRID<T>*)object)->grid;
    const PhysBAM::OCTREE_CELL<T> *cell=grid.Cell_Pointer_From_Index()(index);
    return object->World_Space_Box(RANGE<VECTOR<float,3> >(cell->Bounding_Box()));
}
//#####################################################################
// OPENGL_SELECTION_OCTREE_FACE
//#####################################################################
template<class T> RANGE<VECTOR<float,3> > OPENGL_SELECTION_OCTREE_FACE<T>::
Bounding_Box() const
{
    PHYSBAM_ASSERT(object);
    const OCTREE_GRID<T> &grid=((OPENGL_OCTREE_GRID<T>*)object)->grid;
    RANGE<VECTOR<float,3> > box(VECTOR<float,3>(grid.Node_Location(node[0])));
    for(int i=1;i<=3;i++) box.Enlarge_Nonempty_Box_To_Include_Point(VECTOR<float,3>(grid.Node_Location(node[i])));
    return object->World_Space_Box(box);
}
//#####################################################################
// OPENGL_SELECTION_OCTREE_NODE
//#####################################################################
template<class T> RANGE<VECTOR<float,3> > OPENGL_SELECTION_OCTREE_NODE<T>::
Bounding_Box() const
{
    PHYSBAM_ASSERT(object);
    RANGE<VECTOR<float,3> > box((VECTOR<float,3>(location)));
    return object->World_Space_Box(box);
}
//#####################################################################
// Print_Selection_Info
//#####################################################################
template<class T> void OPENGL_OCTREE_GRID<T>::
Print_Selection_Info(std::ostream& output_stream,OPENGL_SELECTION* selection) const
{
    if(current_selection && current_selection->type==OPENGL_SELECTION::OCTREE_NODE){
        int index=((OPENGL_SELECTION_OCTREE_NODE<T>*)current_selection)->index;
        output_stream<<"Selected node "<<index<<" ("<<grid.Node_Location(index)<<")"<<std::endl;}
    else if(current_selection && current_selection->type==OPENGL_SELECTION::OCTREE_FACE){
        OPENGL_SELECTION_OCTREE_FACE<T>* real_selection=(OPENGL_SELECTION_OCTREE_FACE<T>*)current_selection;
        int index=real_selection->index;
        output_stream<<"Selected face "<<index<<" ("<<grid.Face_Location(index)<<")"<<std::endl;
        output_stream<<"Corner nodes "<<real_selection->node[0]<<", "<<real_selection->node[1]<<", "<<real_selection->node[2]<<", "<<real_selection->node[3]<<std::endl;
        output_stream<<"Neighbor cells: "<<real_selection->cell[0]<<", "<<real_selection->cell[1]<<std::endl;}
    else if(current_selection && current_selection->type==OPENGL_SELECTION::OCTREE_CELL){
        int index=((OPENGL_SELECTION_OCTREE_CELL<T>*)current_selection)->index;
        OCTREE_CELL<T>* cell=grid.Cell_Pointer_From_Index()(index);
        output_stream<<"Selected cell "<<index<<" ("<<cell->Center()<<")"<<std::endl;}
}

//#####################################################################
template class OPENGL_OCTREE_GRID<float>;
#ifndef COMPILE_WITHOUT_DOUBLE_SUPPORT
template class OPENGL_OCTREE_GRID<double>;
#endif
#endif
