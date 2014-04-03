//#####################################################################
// Copyright 2004-2009, Eran Guendelman, Andrew Selle.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#include <PhysBAM_Tools/Grids_Uniform_Arrays/ARRAYS_ND.h>
#include <PhysBAM_Tools/Grids_Uniform_Interpolation/LINEAR_INTERPOLATION_UNIFORM.h>
#include <PhysBAM_Tools/Math_Tools/RANGE.h>
#include <PhysBAM_Rendering/PhysBAM_OpenGL/OpenGL/OPENGL_COLOR_RAMP.h>
#include <PhysBAM_Rendering/PhysBAM_OpenGL/OpenGL/OPENGL_POINTS_3D.h>
#include <PhysBAM_Rendering/PhysBAM_OpenGL/OpenGL/OPENGL_SCALAR_FIELD_3D.h>
#include <PhysBAM_Rendering/PhysBAM_OpenGL/OpenGL/OPENGL_TEXTURED_RECT.h>
#include <PhysBAM_Rendering/PhysBAM_OpenGL/OpenGL/OPENGL_UNIFORM_SLICE.h>
#include <PhysBAM_Rendering/PhysBAM_OpenGL/OpenGL/OPENGL_WORLD.h>
#include <PhysBAM_Rendering/PhysBAM_OpenGL/OpenGL_Components/OPENGL_COMPONENT_PARTICLES_3D.h>
namespace PhysBAM{

template<class T,class T2> OPENGL_SCALAR_FIELD_3D<T,T2>::
OPENGL_SCALAR_FIELD_3D(const GRID<TV> &grid_input,ARRAY<T2,VECTOR<int,3> > &values_input,OPENGL_COLOR_MAP<T2> *color_map_input,DRAW_MODE draw_mode_input)
    :grid(grid_input),values(values_input),current_color_map(1),opengl_textured_rect(0),opengl_points(0),smooth_slice_texture(false),scale_range(false)
{
    PHYSBAM_ASSERT(color_map_input);
    Initialize_Color_Maps(color_map_input);
    Set_Draw_Mode(draw_mode_input);
}

template<class T,class T2> OPENGL_SCALAR_FIELD_3D<T,T2>::
~OPENGL_SCALAR_FIELD_3D()
{
    Delete_Textured_Rect();
    Delete_Points();
    color_maps.Delete_Pointers_And_Clean_Memory();
}

template<> void OPENGL_SCALAR_FIELD_3D<float,bool>::
Initialize_Color_Maps(OPENGL_COLOR_MAP<bool>* color_map_input)
{
    color_maps.Append(color_map_input);
}

template<> void OPENGL_SCALAR_FIELD_3D<double,bool>::
Initialize_Color_Maps(OPENGL_COLOR_MAP<bool>* color_map_input)
{
    color_maps.Append(color_map_input);
}

template<> void OPENGL_SCALAR_FIELD_3D<float,int>::
Initialize_Color_Maps(OPENGL_COLOR_MAP<int>* color_map_input)
{
    color_maps.Append(color_map_input);
}

template<> void OPENGL_SCALAR_FIELD_3D<double,int>::
Initialize_Color_Maps(OPENGL_COLOR_MAP<int>* color_map_input)
{
    color_maps.Append(color_map_input);
}

template<class T,class T2> void OPENGL_SCALAR_FIELD_3D<T,T2>::
Initialize_Color_Maps(OPENGL_COLOR_MAP<T2>* color_map_input)
{
    color_maps.Append(color_map_input);
    color_maps.Append(OPENGL_COLOR_RAMP<T2>::Matlab_Jet(8e4,1e6));
    color_maps.Append(OPENGL_COLOR_RAMP<T2>::Matlab_Hot(8e4,1e6));
}

template<class T,class T2> void OPENGL_SCALAR_FIELD_3D<T,T2>::
Delete_Textured_Rect()
{
    if(opengl_textured_rect){delete opengl_textured_rect->texture;}
    delete opengl_textured_rect;opengl_textured_rect=0;
}

template<class T,class T2> void OPENGL_SCALAR_FIELD_3D<T,T2>::
Delete_Points()
{
    if(opengl_points) delete &opengl_points->points;
    delete opengl_points;opengl_points=0;
}

template<> void OPENGL_SCALAR_FIELD_3D<float,bool>::
Set_Scale_Range(const bool range_min,const bool range_max)
{PHYSBAM_FATAL_ERROR();}

template<> void OPENGL_SCALAR_FIELD_3D<double,bool>::
Set_Scale_Range(const bool range_min,const bool range_max)
{PHYSBAM_FATAL_ERROR();}

template<> void OPENGL_SCALAR_FIELD_3D<float,int>::
Set_Scale_Range(const int range_min,const int range_max)
{PHYSBAM_FATAL_ERROR();}

template<> void OPENGL_SCALAR_FIELD_3D<double,int>::
Set_Scale_Range(const int range_min,const int range_max)
{PHYSBAM_FATAL_ERROR();}

template<class T,class T2> void OPENGL_SCALAR_FIELD_3D<T,T2>::
Set_Scale_Range(const T2 range_min,const T2 range_max)
{
    scale_range=true;
    scale_range_min=range_min;
    T2 range_length=(range_max-range_min);
    scale_range_dx=range_length>0?(T2)1/range_length:(T2)0;
}

template<class T,class T2> void OPENGL_SCALAR_FIELD_3D<T,T2>::
Reset_Scale_Range()
{
    scale_range=false;
}

template<> bool OPENGL_SCALAR_FIELD_3D<float,bool>::
Pre_Map_Value(const bool value) const
{
    return value;
}

template<> bool OPENGL_SCALAR_FIELD_3D<double,bool>::
Pre_Map_Value(const bool value) const
{
    return value;
}


template<class T,class T2> T2 OPENGL_SCALAR_FIELD_3D<T,T2>::
Pre_Map_Value(const T2 value) const
{
    if(!scale_range) return value;
    else return (value-scale_range_min)*scale_range_dx; 
}

#define DO_COLOR \
{\
    T2 value=values(i,j,k);\
    if(value==T2()) continue;\
    color_maps(current_color_map)->Lookup(Pre_Map_Value(value)).Send_To_GL_Pipeline(); \
}

template<class T,class T2> void OPENGL_SCALAR_FIELD_3D<T,T2>::
Display(const int in_color) const
{
    if(values.counts.x==0 || values.counts.y==0) return;
    OPENGL_UNIFORM_SLICE* slice=(OPENGL_UNIFORM_SLICE*)this->slice;

    glMatrixMode(GL_MODELVIEW);
    glPushMatrix();
    Send_Transform_To_GL_Pipeline();

    if(draw_mode==DRAW_TEXTURE){
        if(!slice || slice->mode==OPENGL_SLICE::NO_SLICE) Display_3D();
        else if(opengl_textured_rect) opengl_textured_rect->Display(in_color);}
    else{
        PHYSBAM_ASSERT(opengl_points);
        if(slice && slice->Is_Slice_Mode()){
            glPushAttrib(GL_ENABLE_BIT);
            slice->Enable_Clip_Planes();}
        opengl_points->Display(in_color);
        if(slice && slice->Is_Slice_Mode()){
            glPopAttrib();}}

    glPopMatrix();
}

template<class T,class T2> void OPENGL_SCALAR_FIELD_3D<T,T2>::
Display_3D() const
{
    glPushAttrib(GL_ENABLE_BIT|GL_DEPTH_BUFFER_BIT|GL_CURRENT_BIT);
    glDisable(GL_LIGHTING);
    glEnable(GL_DEPTH_TEST);
    glEnable(GL_BLEND);
    glDepthMask(GL_FALSE);

    VECTOR<float,3> view_forward,view_up,view_right;
    OPENGL_WORLD::Singleton()->Get_View_Frame(view_forward,view_up,view_right);
    int dominant_axis=1;

    if(dominant_axis==1){
        if(view_forward[1]>0){
            OpenGL_Begin(GL_QUADS);
            for(int i=grid.counts.x;i>=1;i--) for(int j=1;j<=grid.counts.y;j++) for(int k=1;k<=grid.counts.z;k++){
                DO_COLOR;
                VECTOR<T,3> pos=grid.X(i,j,k);
                glVertex3f(pos.x,pos.y-0.5*grid.dX.y,pos.z-0.5*grid.dX.z);
                glVertex3f(pos.x,pos.y-0.5*grid.dX.y,pos.z+0.5*grid.dX.z);
                glVertex3f(pos.x,pos.y+0.5*grid.dX.y,pos.z+0.5*grid.dX.z);
                glVertex3f(pos.x,pos.y+0.5*grid.dX.y,pos.z-0.5*grid.dX.z);}
            OpenGL_End();}
        else{
            OpenGL_Begin(GL_QUADS);
            for(int i=1;i<=grid.counts.x;i++) for(int j=1;j<=grid.counts.y;j++) for(int k=1;k<=grid.counts.z;k++){
                DO_COLOR;
                VECTOR<T,3> pos=grid.X(i,j,k);
                glVertex3f(pos.x,pos.y-0.5*grid.dX.y,pos.z+0.5*grid.dX.z);
                glVertex3f(pos.x,pos.y-0.5*grid.dX.y,pos.z-0.5*grid.dX.z);
                glVertex3f(pos.x,pos.y+0.5*grid.dX.y,pos.z-0.5*grid.dX.z);
                glVertex3f(pos.x,pos.y+0.5*grid.dX.y,pos.z+0.5*grid.dX.z);}
            OpenGL_End();}}
    else if(dominant_axis==2){
        if(view_forward[2]>0){
            OpenGL_Begin(GL_QUADS);
            for(int j=grid.counts.y;j>=1;j--) for(int i=1;i<=grid.counts.x;i++) for(int k=1;k<=grid.counts.z;k++){
                DO_COLOR;
                VECTOR<T,3> pos=grid.X(i,j,k);
                glVertex3f(pos.x-0.5*grid.dX.x,pos.y,pos.z-0.5*grid.dX.z);
                glVertex3f(pos.x+0.5*grid.dX.x,pos.y,pos.z-0.5*grid.dX.z);
                glVertex3f(pos.x+0.5*grid.dX.x,pos.y,pos.z+0.5*grid.dX.z);
                glVertex3f(pos.x-0.5*grid.dX.x,pos.y,pos.z+0.5*grid.dX.z);}
            OpenGL_End();}
        else{
            OpenGL_Begin(GL_QUADS);
            for(int j=1;j<=grid.counts.y;j++) for(int i=1;i<=grid.counts.x;i++) for(int k=1;k<=grid.counts.z;k++){
                DO_COLOR;
                VECTOR<T,3> pos=grid.X(i,j,k);
                glVertex3f(pos.x-0.5*grid.dX.x,pos.y,pos.z+0.5*grid.dX.z);
                glVertex3f(pos.x+0.5*grid.dX.x,pos.y,pos.z+0.5*grid.dX.z);
                glVertex3f(pos.x+0.5*grid.dX.x,pos.y,pos.z-0.5*grid.dX.z);
                glVertex3f(pos.x-0.5*grid.dX.x,pos.y,pos.z-0.5*grid.dX.z);}
            OpenGL_End();}}
    else if(dominant_axis==3){
        if(view_forward[3]>0){
            OpenGL_Begin(GL_QUADS);
            for(int k=grid.counts.z;k>=1;k--) for(int i=1;i<=grid.counts.x;i++) for(int j=1;j<=grid.counts.y;j++){
                DO_COLOR;
                VECTOR<T,3> pos=grid.X(i,j,k);
                glVertex3f(pos.x-0.5*grid.dX.x,pos.y-0.5*grid.dX.y,pos.z);
                glVertex3f(pos.x-0.5*grid.dX.x,pos.y+0.5*grid.dX.y,pos.z);
                glVertex3f(pos.x+0.5*grid.dX.x,pos.y+0.5*grid.dX.y,pos.z);
                glVertex3f(pos.x+0.5*grid.dX.x,pos.y-0.5*grid.dX.y,pos.z);}
            OpenGL_End();}
        else{
            OpenGL_Begin(GL_QUADS);
            for(int k=1;k<=grid.counts.z;k++) for(int i=1;i<=grid.counts.x;i++) for(int j=1;j<=grid.counts.y;j++){
                DO_COLOR;
                VECTOR<T,3> pos=grid.X(i,j,k);
                glVertex3f(pos.x-0.5*grid.dX.x,pos.y+0.5*grid.dX.y,pos.z);
                glVertex3f(pos.x-0.5*grid.dX.x,pos.y-0.5*grid.dX.y,pos.z);
                glVertex3f(pos.x+0.5*grid.dX.x,pos.y-0.5*grid.dX.y,pos.z);
                glVertex3f(pos.x+0.5*grid.dX.x,pos.y+0.5*grid.dX.y,pos.z);}
            OpenGL_End();}}

    glPopAttrib();
}

template<class T,class T2> RANGE<VECTOR<float,3> > OPENGL_SCALAR_FIELD_3D<T,T2>::
Bounding_Box() const
{
    if(slice && slice->Is_Slice_Mode() && opengl_textured_rect) return opengl_textured_rect->Bounding_Box();
    else return World_Space_Box((RANGE<VECTOR<float,3> >)grid.domain);
}

template<class T,class T2> void OPENGL_SCALAR_FIELD_3D<T,T2>::
Set_Draw_Mode(DRAW_MODE draw_mode_input)
{
    draw_mode=draw_mode_input;

    if(draw_mode==DRAW_TEXTURE){
        Delete_Points();
        if(!opengl_textured_rect) opengl_textured_rect=new OPENGL_TEXTURED_RECT;}
    else{
        Delete_Textured_Rect();
        if(!opengl_points) opengl_points=new OPENGL_POINTS_3D<T>(*new ARRAY<VECTOR<T,3> >);}

    Update();
}

template<class T,class T2> void OPENGL_SCALAR_FIELD_3D<T,T2>::
Update()
{
    if(draw_mode==DRAW_TEXTURE){if(slice && slice->Is_Slice_Mode()) Update_Slice();}
    else Update_Points();
}

template<class T2,class T> static void Print_Selection_Info_Helper(std::ostream& output_stream,OPENGL_SELECTION_COMPONENT_PARTICLES_3D<T>* selection,const GRID<VECTOR<T,3> >& grid,ARRAY<T2,VECTOR<int,3> >& values)
{
    output_stream<<" @ particle = "<<LINEAR_INTERPOLATION_UNIFORM<GRID<VECTOR<T,3> >,T2>().Clamped_To_Array(grid,values,selection->location);
}
// no interpolation for bool's and int's
template<class T> static void Print_Selection_Info_Helper(std::ostream& output_stream,OPENGL_SELECTION_COMPONENT_PARTICLES_3D<T>* selection,const GRID<VECTOR<T,3> >&,ARRAY<bool,VECTOR<int,3> >& values){}
template<class T> static void Print_Selection_Info_Helper(std::ostream& output_stream,OPENGL_SELECTION_COMPONENT_PARTICLES_3D<T>* selection,const GRID<VECTOR<T,3> >&,ARRAY<int,VECTOR<int,3> >& values){}

template<class T,class T2> void OPENGL_SCALAR_FIELD_3D<T,T2>::
Print_Selection_Info(std::ostream& output_stream,OPENGL_SELECTION* current_selection) const
{
    if(current_selection && current_selection->type==OPENGL_SELECTION::GRID_CELL_3D && grid.Is_MAC_Grid()){
        VECTOR<int,3> index=((OPENGL_SELECTION_GRID_CELL_3D<T>*)current_selection)->index;
        if(values.Valid_Index(index)) output_stream<<values(index);}
    if(current_selection && current_selection->type==OPENGL_SELECTION::GRID_NODE_3D && !grid.Is_MAC_Grid()){
        VECTOR<int,3> index=((OPENGL_SELECTION_GRID_NODE_3D<T>*)current_selection)->index;
        if(values.Valid_Index(index))output_stream<<values(index);}
    if(current_selection && current_selection->type==OPENGL_SELECTION::COMPONENT_PARTICLES_3D){
        OPENGL_SELECTION_COMPONENT_PARTICLES_3D<T> *selection=(OPENGL_SELECTION_COMPONENT_PARTICLES_3D<T>*)current_selection;
        Print_Selection_Info_Helper(output_stream,selection,grid,values);}
    output_stream<<std::endl;
}

namespace{
template<class T>
void Update_Slice_Helper(OPENGL_SCALAR_FIELD_3D<T,bool>* self,int tex_width,int tex_height){
    OPENGL_UNIFORM_SLICE* slice=(OPENGL_UNIFORM_SLICE*)self->slice;
    VECTOR<int,3> domain_start(self->values.domain.min_corner.x,self->values.domain.min_corner.y,self->values.domain.min_corner.z),domain_end(self->values.domain.max_corner.x,self->values.domain.max_corner.y,self->values.domain.max_corner.z);
    if(!self->opengl_textured_rect->texture || self->opengl_textured_rect->texture->width!=tex_width || self->opengl_textured_rect->texture->height!=tex_height){
        delete self->opengl_textured_rect->texture;
        self->opengl_textured_rect->texture=new OPENGL_TEXTURE();
        self->opengl_textured_rect->texture->Initialize(tex_width,tex_height);
        self->opengl_textured_rect->texture->Set_Smooth_Shading(self->smooth_slice_texture);}

    OPENGL_COLOR* bitmap=new OPENGL_COLOR[tex_width*tex_height];
    OPENGL_COLOR_MAP<bool>* color_map=self->color_maps(self->current_color_map);
    for(int i=1;i<=tex_width;i++)
        for(int j=1;j<=tex_height;j++){
            bool value=bool();
            switch (slice->axis){
                case 1: value=self->values(slice->index,domain_start.y+j-1,domain_start.z+tex_width-i);break;
                case 2: value=self->values(domain_start.x+i-1,slice->index,domain_start.z+tex_height-j);break;
                case 3: value=self->values(domain_start.x+i-1,domain_start.y+j-1,slice->index);break;}
            int idx=(j-1)*tex_width+i-1;
            bitmap[idx]=color_map->Lookup(self->Pre_Map_Value(value));}

    self->opengl_textured_rect->texture->Update_Texture(bitmap);

    delete[]bitmap;
}

template<class T>
void Update_Slice_Helper(OPENGL_SCALAR_FIELD_3D<T,int>* self,int tex_width,int tex_height){
    OPENGL_UNIFORM_SLICE* slice=(OPENGL_UNIFORM_SLICE*)self->slice;
    VECTOR<int,3> domain_start(self->values.domain.min_corner.x,self->values.domain.min_corner.y,self->values.domain.min_corner.z),domain_end(self->values.domain.max_corner.x,self->values.domain.max_corner.y,self->values.domain.max_corner.z);
    if(!self->opengl_textured_rect->texture || self->opengl_textured_rect->texture->width!=tex_width || self->opengl_textured_rect->texture->height!=tex_height){
        delete self->opengl_textured_rect->texture;
        self->opengl_textured_rect->texture=new OPENGL_TEXTURE();
        self->opengl_textured_rect->texture->Initialize(tex_width,tex_height);
        self->opengl_textured_rect->texture->Set_Smooth_Shading(self->smooth_slice_texture);}

    OPENGL_COLOR* bitmap=new OPENGL_COLOR[tex_width*tex_height];
    OPENGL_COLOR_MAP<int>* color_map=self->color_maps(self->current_color_map);
    for(int i=1;i<=tex_width;i++)
        for(int j=1;j<=tex_height;j++){
            int value=int();
            switch (slice->axis){
                case 1: value=self->values(slice->index,domain_start.y+j-1,domain_start.z+tex_width-i);break;
                case 2: value=self->values(domain_start.x+i-1,slice->index,domain_start.z+tex_height-j);break;
                case 3: value=self->values(domain_start.x+i-1,domain_start.y+j-1,slice->index);break;}
            int idx=(j-1)*tex_width+i-1;
            bitmap[idx]=color_map->Lookup(self->Pre_Map_Value(value));}

    self->opengl_textured_rect->texture->Update_Texture(bitmap);

    delete[]bitmap;
}

template<class T,class T2>
void Update_Slice_Helper(OPENGL_SCALAR_FIELD_3D<T,T2>* self,int tex_width,int tex_height){
    typedef VECTOR<T,3> TV;
    OPENGL_UNIFORM_SLICE* slice=(OPENGL_UNIFORM_SLICE*)self->slice;
    int k=1;
    tex_width=k*tex_width;
    tex_height=k*tex_height;

    if(!self->opengl_textured_rect->texture || self->opengl_textured_rect->texture->width!=tex_width || self->opengl_textured_rect->texture->height!=tex_height){
        delete self->opengl_textured_rect->texture;
        self->opengl_textured_rect->texture=new OPENGL_TEXTURE();
        self->opengl_textured_rect->texture->Initialize(tex_width,tex_height);
        self->opengl_textured_rect->texture->Set_Smooth_Shading(self->smooth_slice_texture);}

    OPENGL_COLOR* bitmap=new OPENGL_COLOR[tex_width*tex_height];
    OPENGL_COLOR_MAP<T2>* color_map=self->color_maps(self->current_color_map);

    LINEAR_INTERPOLATION_UNIFORM<GRID<TV>,T2> interpolation;
    for(int i=1;i<=tex_width;i++)
        for(int j=1;j<=tex_height;j++){
            T2 value=T2();
            TV location;
            switch (slice->axis){
                case 1: location=TV(self->grid.X(slice->index,1,1).x,
                    self->grid.domain.min_corner.y+j*(self->grid.domain.max_corner.y-self->grid.domain.min_corner.y)/tex_height,
                    self->grid.domain.max_corner.z-i*(self->grid.domain.max_corner.z-self->grid.domain.min_corner.z)/tex_width);
                    break;
                case 2: location=TV(self->grid.domain.min_corner.x+i*(self->grid.domain.max_corner.x-self->grid.domain.min_corner.x)/tex_width,
                    self->grid.X(1,slice->index,1).y,
                    self->grid.domain.max_corner.z-j*(self->grid.domain.max_corner.z-self->grid.domain.min_corner.z)/tex_height);
                    break;
                case 3: location=TV(self->grid.domain.min_corner.x+i*(self->grid.domain.max_corner.x-self->grid.domain.min_corner.x)/tex_width,
                    self->grid.domain.min_corner.y+j*(self->grid.domain.max_corner.y-self->grid.domain.min_corner.y)/tex_height,
                    self->grid.X(1,1,slice->index).z);
                    break;}
            value=interpolation.Clamped_To_Array(self->grid,self->values,location);
            int idx=(j-1)*tex_width+i-1;
            bitmap[idx]=color_map->Lookup(self->Pre_Map_Value(value));}

    self->opengl_textured_rect->texture->Update_Texture(bitmap);
    delete []bitmap;
}
}

template<class T,class T2> void OPENGL_SCALAR_FIELD_3D<T,T2>::
Update_Slice()
{
    PHYSBAM_ASSERT(this->slice);
    OPENGL_UNIFORM_SLICE* slice=(OPENGL_UNIFORM_SLICE*)this->slice;
    VECTOR<int,3> domain_start(values.domain.min_corner.x,values.domain.min_corner.y,values.domain.min_corner.z),domain_end(values.domain.max_corner.x,values.domain.max_corner.y,values.domain.max_corner.z);
    if((slice->mode==OPENGL_SLICE::CELL_SLICE && (grid.MAC_offset==0 || slice->index<domain_start[slice->axis] || slice->index>domain_end[slice->axis])) ||
        (slice->mode==OPENGL_SLICE::NODE_SLICE && (grid.MAC_offset==0.5 || slice->index<domain_start[slice->axis] || slice->index>domain_end[slice->axis]))){
        // Currently we don't draw anything ifthe slice doesn't match where the scalar field lives
        Delete_Textured_Rect();
        return;}

    if(!opengl_textured_rect) opengl_textured_rect=new OPENGL_TEXTURED_RECT();

    // Handle values arrays which are not (1,m)(1,n)
    VECTOR<T,3> half_dX=(T)0.5*grid.dX;
    RANGE<VECTOR<int,3> > domain_indices(values.Domain_Indices());
    RANGE<VECTOR<T,3> > domain(grid.X(domain_indices.min_corner)-half_dX,grid.X(domain_indices.max_corner)+half_dX);
    opengl_textured_rect->frame->t=VECTOR<float,3>(domain.Center());

    // rectangle will face you ifyou're looking down the positive axis
    int tex_width=0,tex_height=0; // texture width and height
    switch (slice->axis){
        case 1:
            opengl_textured_rect->frame->t.x=grid.Axis_X(slice->index,slice->axis);
            opengl_textured_rect->frame->r=ROTATION<VECTOR<float,3> >(0.5*pi,VECTOR<float,3>(0,1,0));
            opengl_textured_rect->width=domain.Edge_Lengths().z;
            opengl_textured_rect->height=domain.Edge_Lengths().y;
            tex_width=values.counts.z;
            tex_height=values.counts.y;
            break;

        case 2:
            opengl_textured_rect->frame->t.y=grid.Axis_X(slice->index,slice->axis);
            opengl_textured_rect->frame->r=ROTATION<VECTOR<float,3> >(-0.5*pi,VECTOR<float,3>(1,0,0));
            opengl_textured_rect->width=domain.Edge_Lengths().x;
            opengl_textured_rect->height=domain.Edge_Lengths().z;
            tex_width=values.counts.x;
            tex_height=values.counts.z;
            break;

        case 3:
            opengl_textured_rect->frame->t.z=grid.Axis_X(slice->index,slice->axis);
            opengl_textured_rect->frame->r=ROTATION<VECTOR<float,3> >();
            opengl_textured_rect->width=domain.Edge_Lengths().x;
            opengl_textured_rect->height=domain.Edge_Lengths().y;
            tex_width=values.counts.x;
            tex_height=values.counts.y;
            break;}

    Update_Slice_Helper(this,tex_width,tex_height);
}

template<class T,class T2> void OPENGL_SCALAR_FIELD_3D<T,T2>::
Slice_Has_Changed()
{
    OPENGL_UNIFORM_SLICE* slice=(OPENGL_UNIFORM_SLICE*)this->slice;
    if(!slice || slice->mode==OPENGL_SLICE::NO_SLICE) // In 3D mode now, can erase textured rect
        Delete_Textured_Rect();
    Update();
}

template<class T,class T2> void OPENGL_SCALAR_FIELD_3D<T,T2>::
Update_Points()
{
    PHYSBAM_ASSERT(opengl_points);
    opengl_points->points.Resize(values.counts.Product());
    int index=1;
    for(int i=values.domain.min_corner.x;i<=values.domain.max_corner.x;i++) for(int j=values.domain.min_corner.y;j<=values.domain.max_corner.y;j++) for(int k=values.domain.min_corner.z;k<=values.domain.max_corner.z;k++){
        opengl_points->points(index)=grid.X(i,j,k);
        opengl_points->Set_Point_Color(index,color_maps(current_color_map)->Lookup(values(i,j,k)));
        index++;}
}

//#####################################################################
// Specialization for bool scalars: only draw point if value is "true"
//#####################################################################
template<> void OPENGL_SCALAR_FIELD_3D<float,bool>::
Update_Points()
{
    PHYSBAM_ASSERT(opengl_points);
    opengl_points->color=color_maps(current_color_map)->Lookup(true);
    opengl_points->points.Resize(values.counts.Product());
    int index=1;
    for(int i=values.domain.min_corner.x;i<=values.domain.max_corner.x;i++) for(int j=values.domain.min_corner.y;j<=values.domain.max_corner.y;j++) for(int k=values.domain.min_corner.z;k<=values.domain.max_corner.z;k++)
        if(values(i,j,k)) opengl_points->points(index++)=grid.X(i,j,k);
    opengl_points->points.Resize(index-1);
}

template<> void OPENGL_SCALAR_FIELD_3D<double,bool>::
Update_Points()
{
    PHYSBAM_ASSERT(opengl_points);
    opengl_points->color=color_maps(current_color_map)->Lookup(true);
    opengl_points->points.Resize(values.counts.Product());
    int index=1;
    for(int i=values.domain.min_corner.x;i<=values.domain.max_corner.x;i++) for(int j=values.domain.min_corner.y;j<=values.domain.max_corner.y;j++) for(int k=values.domain.min_corner.z;k<=values.domain.max_corner.z;k++)
        if(values(i,j,k)) opengl_points->points(index++)=grid.X(i,j,k);
    opengl_points->points.Resize(index-1);
}

template<class T,class T2> void OPENGL_SCALAR_FIELD_3D<T,T2>::
Toggle_Draw_Mode()
{
    DRAW_MODE new_draw_mode=(DRAW_MODE)(((int)draw_mode+1)%2);
    Set_Draw_Mode(new_draw_mode);
}

template<class T,class T2> void OPENGL_SCALAR_FIELD_3D<T,T2>::
Set_Smooth_Slice_Texture(bool smooth_slice_texture_input)
{
    smooth_slice_texture=smooth_slice_texture_input;
    if(draw_mode==DRAW_TEXTURE && opengl_textured_rect && opengl_textured_rect->texture) opengl_textured_rect->texture->Set_Smooth_Shading(smooth_slice_texture);
}

template<class T,class T2> void OPENGL_SCALAR_FIELD_3D<T,T2>::
Toggle_Smooth_Slice_Texture()
{
    Set_Smooth_Slice_Texture(!smooth_slice_texture);
}

template<class T,class T2> void OPENGL_SCALAR_FIELD_3D<T,T2>::
Toggle_Color_Map()
{
    current_color_map=current_color_map%color_maps.m+1;
    Update();
}

template class OPENGL_SCALAR_FIELD_3D<float,int>;
template class OPENGL_SCALAR_FIELD_3D<float,bool>;
template class OPENGL_SCALAR_FIELD_3D<float,float>;
#ifndef COMPILE_WITHOUT_DOUBLE_SUPPORT
template class OPENGL_SCALAR_FIELD_3D<double,int>;
template class OPENGL_SCALAR_FIELD_3D<double,bool>;
template class OPENGL_SCALAR_FIELD_3D<double,double>;
#endif
}
