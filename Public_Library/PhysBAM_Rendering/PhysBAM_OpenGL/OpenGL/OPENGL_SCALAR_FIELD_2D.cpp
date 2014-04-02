//#####################################################################
// Copyright 2004-2009, Jon Gretarsson, Eran Guendelman, Nipun Kwatra, Andrew Selle.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#include <PhysBAM_Tools/Grids_Uniform_Arrays/ARRAYS_ND.h>
#include <PhysBAM_Tools/Grids_Uniform_Interpolation/LINEAR_INTERPOLATION_UNIFORM.h>
#include <PhysBAM_Geometry/Grids_Uniform_Computations/DUALCONTOUR_2D.h>
#include <PhysBAM_Geometry/Topology_Based_Geometry/SEGMENTED_CURVE_2D.h>
#include <PhysBAM_Rendering/PhysBAM_OpenGL/OpenGL/OPENGL_COLOR_RAMP.h>
#include <PhysBAM_Rendering/PhysBAM_OpenGL/OpenGL/OPENGL_GRID_2D.h>
#include <PhysBAM_Rendering/PhysBAM_OpenGL/OpenGL/OPENGL_POINTS_2D.h>
#include <PhysBAM_Rendering/PhysBAM_OpenGL/OpenGL/OPENGL_SCALAR_FIELD_2D.h>
#include <PhysBAM_Rendering/PhysBAM_OpenGL/OpenGL/OPENGL_SEGMENTED_CURVE_2D.h>
#include <PhysBAM_Rendering/PhysBAM_OpenGL/OpenGL/OPENGL_TEXTURED_RECT.h>
#include <PhysBAM_Rendering/PhysBAM_OpenGL/OpenGL_Components/OPENGL_COMPONENT_PARTICLES_2D.h>
namespace PhysBAM{

template<class T,class T2> OPENGL_SCALAR_FIELD_2D<T,T2>::
OPENGL_SCALAR_FIELD_2D(GRID<TV> &grid_input,ARRAY<T2,VECTOR<int,2> > &values_input,OPENGL_COLOR_MAP<T2>* color_map_input,DRAW_MODE draw_mode_input)
    : grid(grid_input),values(values_input),active_cells(0),draw_ghost_values(true),current_color_map(1),opengl_textured_rect(0),opengl_points(0),scale_range(false)
{
    PHYSBAM_ASSERT(color_map_input);
    Initialize_Color_Maps(color_map_input);
    Set_Draw_Mode(draw_mode_input);
}

template<class T,class T2> OPENGL_SCALAR_FIELD_2D<T,T2>::
OPENGL_SCALAR_FIELD_2D(GRID<TV> &grid_input,ARRAY<T2,VECTOR<int,2> > &values_input,OPENGL_COLOR_MAP<T2>* color_map_input,ARRAY<bool,VECTOR<int,2> >* active_cells_input,DRAW_MODE draw_mode_input)
    : grid(grid_input),values(values_input),active_cells(active_cells_input),draw_ghost_values(true),current_color_map(1),opengl_textured_rect(0),opengl_points(0),scale_range(false)
{
    PHYSBAM_ASSERT(color_map_input);
    Initialize_Color_Maps(color_map_input);
    Set_Draw_Mode(draw_mode_input);
}

template<class T,class T2> OPENGL_SCALAR_FIELD_2D<T,T2>::
~OPENGL_SCALAR_FIELD_2D()
{
    if (opengl_points) delete &opengl_points->points;
    delete opengl_points;
    if (opengl_textured_rect) delete opengl_textured_rect->texture;
    delete opengl_textured_rect;
    color_maps.Delete_Pointers_And_Clean_Memory();
    contour_curves.Delete_Pointers_And_Clean_Memory();
}

template<class T,class T2> void OPENGL_SCALAR_FIELD_2D<T,T2>::
Set_Uniform_Contour_Values(const T2 min_value,const T2 max_value,const T2 increment)
{
    contour_values.Remove_All();
    for(T2 value=min_value;value<=max_value;value+=increment) contour_values.Append(value);
    if(draw_mode==DRAW_CONTOURS) Update_Contour_Curves();
}

template<> void OPENGL_SCALAR_FIELD_2D<float,bool>::
Initialize_Color_Maps(OPENGL_COLOR_MAP<bool>* color_map_input)
{
    color_maps.Append(color_map_input);
}

template<> void OPENGL_SCALAR_FIELD_2D<double,bool>::
Initialize_Color_Maps(OPENGL_COLOR_MAP<bool>* color_map_input)
{
    color_maps.Append(color_map_input);
}

template<> void OPENGL_SCALAR_FIELD_2D<float,int>::
Initialize_Color_Maps(OPENGL_COLOR_MAP<int>* color_map_input)
{
    color_maps.Append(color_map_input);
}

template<> void OPENGL_SCALAR_FIELD_2D<double,int>::
Initialize_Color_Maps(OPENGL_COLOR_MAP<int>* color_map_input)
{
    color_maps.Append(color_map_input);
}

template<class T,class T2> void OPENGL_SCALAR_FIELD_2D<T,T2>::
Initialize_Color_Maps(OPENGL_COLOR_MAP<T2>* color_map_input)
{
    color_maps.Append(color_map_input);
    color_maps.Append(OPENGL_COLOR_RAMP<T2>::Matlab_Jet(0,1));
    color_maps.Append(OPENGL_COLOR_RAMP<T2>::Matlab_Hot(0,1));
}

template<> void OPENGL_SCALAR_FIELD_2D<float,bool>::
Set_Scale_Range(const bool range_min,const bool range_max)
{PHYSBAM_FATAL_ERROR();}

template<> void OPENGL_SCALAR_FIELD_2D<double,bool>::
Set_Scale_Range(const bool range_min,const bool range_max)
{PHYSBAM_FATAL_ERROR();}

template<> void OPENGL_SCALAR_FIELD_2D<float,int>::
Set_Scale_Range(const int range_min,const int range_max)
{PHYSBAM_FATAL_ERROR();}

template<> void OPENGL_SCALAR_FIELD_2D<double,int>::
Set_Scale_Range(const int range_min,const int range_max)
{PHYSBAM_FATAL_ERROR();}

template<class T,class T2> void OPENGL_SCALAR_FIELD_2D<T,T2>::
Set_Scale_Range(const T2 range_min,const T2 range_max)
{
    scale_range=true;
    scale_range_min=range_min;
    T2 range_length=(range_max-range_min);
    scale_range_dx=range_length>1e-10?(T2)1/range_length:(T2)0;
}

template<class T,class T2> void OPENGL_SCALAR_FIELD_2D<T,T2>::
Reset_Scale_Range()
{
    scale_range=false;
}

template<> bool OPENGL_SCALAR_FIELD_2D<float,bool>::
Pre_Map_Value(const bool value) const
{
    return value;
}

template<> bool OPENGL_SCALAR_FIELD_2D<double,bool>::
Pre_Map_Value(const bool value) const
{
    return value;
}


template<class T,class T2> T2 OPENGL_SCALAR_FIELD_2D<T,T2>::
Pre_Map_Value(const T2 value) const
{
    if(!scale_range) return value;
    else return (value-scale_range_min)*scale_range_dx; 
}

template<class T,class T2> void OPENGL_SCALAR_FIELD_2D<T,T2>::
Display(const int in_color) const
{
    if (draw_mode==DRAW_TEXTURE) {
        PHYSBAM_ASSERT(opengl_textured_rect);
        opengl_textured_rect->Display(in_color);
    } else if (draw_mode==DRAW_POINTS) {
        PHYSBAM_ASSERT(opengl_points);
        opengl_points->Display(in_color);
    } else {
        for(int i=1;i<=contour_curves.m;i++) contour_curves(i)->Display();
    }
}

template<class T,class T2> RANGE<VECTOR<float,3> > OPENGL_SCALAR_FIELD_2D<T,T2>::
Bounding_Box() const
{
    // May not be the exact bounds, but close enough...
    return World_Space_Box(RANGE<VECTOR<float,3> >(grid.domain.min_corner.x,grid.domain.max_corner.x,grid.domain.min_corner.y,grid.domain.max_corner.y,0,0));
}

template<class T,class T2> void OPENGL_SCALAR_FIELD_2D<T,T2>::
Set_Draw_Mode(DRAW_MODE draw_mode_input)
{
    draw_mode = draw_mode_input;

    if (draw_mode!=DRAW_TEXTURE && opengl_textured_rect) {
        delete opengl_textured_rect->texture; delete opengl_textured_rect; opengl_textured_rect=0;
    }
    if (draw_mode!=DRAW_POINTS && opengl_points) {
        delete &opengl_points->points; delete opengl_points; opengl_points=0;
    }

    if (draw_mode==DRAW_TEXTURE) {
        if (!opengl_textured_rect) opengl_textured_rect=new OPENGL_TEXTURED_RECT;
    } else if (draw_mode==DRAW_POINTS) {
        if (!opengl_points) opengl_points=new OPENGL_POINTS_2D<T>(*new ARRAY<VECTOR<T,2> >);
    }

    Update();
}

template<class T,class T2> void OPENGL_SCALAR_FIELD_2D<T,T2>::
Update()
{
    VECTOR<int,2> start_index(values.domain.min_corner.x,values.domain.min_corner.y),end_index(values.domain.max_corner.x,values.domain.max_corner.y);
    if(!draw_ghost_values){start_index=clamp_min(start_index,VECTOR<int,2>(1,1));end_index=clamp_max(end_index,grid.counts);}
    if (draw_mode==DRAW_TEXTURE) Update_Texture(start_index,end_index);
    else if (draw_mode==DRAW_POINTS) Update_Points(start_index,end_index);
    else Update_Contour_Curves();
}


template<class T2,class T> static void Print_Selection_Info_Helper(std::ostream& output_stream,OPENGL_SELECTION_COMPONENT_PARTICLES_2D<T>* selection,const GRID<VECTOR<T,2> >& grid,ARRAY<T2,VECTOR<int,2> >& values)
{
    output_stream<<" @ particle = "<<LINEAR_INTERPOLATION_UNIFORM<GRID<VECTOR<T,2> >,T2>().Clamped_To_Array(grid,values,selection->location);
}
// no interpolation for bool's and int's
template<class T> static void Print_Selection_Info_Helper(std::ostream& output_stream,OPENGL_SELECTION_COMPONENT_PARTICLES_2D<T>* selection,const GRID<VECTOR<T,2> >&,ARRAY<bool,VECTOR<int,2> >& values){}
template<class T> static void Print_Selection_Info_Helper(std::ostream& output_stream,OPENGL_SELECTION_COMPONENT_PARTICLES_2D<T>* selection,const GRID<VECTOR<T,2> >&,ARRAY<int,VECTOR<int,2> >& values){}

template<class T,class T2> void OPENGL_SCALAR_FIELD_2D<T,T2>::
Print_Selection_Info(std::ostream& output_stream,OPENGL_SELECTION* current_selection) const
{
    if(current_selection && current_selection->type==OPENGL_SELECTION::GRID_CELL_2D && grid.Is_MAC_Grid()){
        VECTOR<int,2> index=((OPENGL_SELECTION_GRID_CELL_2D<T>*)current_selection)->index;
        if(values.Valid_Index(index)) output_stream<<values(index);}
    if(current_selection && current_selection->type==OPENGL_SELECTION::GRID_NODE_2D && !grid.Is_MAC_Grid()){
        VECTOR<int,2> index=((OPENGL_SELECTION_GRID_NODE_2D<T>*)current_selection)->index;
        if(values.Valid_Index(index))output_stream<<values(index);}
    if(current_selection && current_selection->type==OPENGL_SELECTION::COMPONENT_PARTICLES_2D){
        OPENGL_SELECTION_COMPONENT_PARTICLES_2D<T> *selection=(OPENGL_SELECTION_COMPONENT_PARTICLES_2D<T>*)current_selection;
        Print_Selection_Info_Helper(output_stream,selection,grid,values);}
    output_stream<<std::endl;
}

template<class T,class T2> void OPENGL_SCALAR_FIELD_2D<T,T2>::
Update_Texture(const VECTOR<int,2>& start_index,const VECTOR<int,2>& end_index)
{
    PHYSBAM_ASSERT(opengl_textured_rect);

    if ((end_index.x-start_index.x+1) == 0 || (end_index.y-start_index.y+1) == 0)
    {
        delete opengl_textured_rect->texture; opengl_textured_rect->texture = 0; return;
    }

    // Handle values arrays which are not (1,m)(1,n)
    VECTOR<T,2> half_dX=(T)0.5*grid.dX;
    RANGE<VECTOR<T,2> > domain(grid.X(start_index)-half_dX,grid.X(end_index)+half_dX);

    // Set underlying OPENGL_OBJECT's transformation
    opengl_textured_rect->frame->t=VECTOR<float,3>(Convert_2d_To_3d(domain.Center()));

    // Set OPENGL_TEXTURED_RECT's dimensions
    opengl_textured_rect->width = domain.Edge_Lengths().x;
    opengl_textured_rect->height = domain.Edge_Lengths().y;

    int tex_width = end_index.x-start_index.x+1;
    int tex_height = end_index.y-start_index.y+1;

    if (!opengl_textured_rect->texture || opengl_textured_rect->texture->width != tex_width || opengl_textured_rect->texture->height != tex_height)
    {
        delete opengl_textured_rect->texture;
        opengl_textured_rect->texture = new OPENGL_TEXTURE();
        opengl_textured_rect->texture->Initialize(tex_width, tex_height);
    }

    OPENGL_COLOR* bitmap = new OPENGL_COLOR[tex_width*tex_height];
    OPENGL_COLOR_MAP<T2>* color_map=color_maps(current_color_map);
    for (int i = 1; i <= tex_width; i++)
        for (int j = 1; j <= tex_height; j++)
        {
            int idx = (j-1)*tex_width + i-1;
            T2 value=Pre_Map_Value(values(start_index.x+i-1,start_index.y+j-1));
            
            OPENGL_COLOR color_value=color_map->Lookup(value);
            if(active_cells && !(*active_cells)(start_index.x+i-1,start_index.y+j-1)) color_value.rgba[3]=0;
            bitmap[idx]=color_value;
        }

    opengl_textured_rect->texture->Update_Texture(bitmap);

    delete[] bitmap;
}

template<class T,class T2> void OPENGL_SCALAR_FIELD_2D<T,T2>::
Update_Points(const VECTOR<int,2>& start_index,const VECTOR<int,2>& end_index)
{
    PHYSBAM_ASSERT(opengl_points);
    opengl_points->points.Resize((end_index.x-start_index.x+1)*(end_index.y-start_index.y+1));
    int index=1;
    OPENGL_COLOR_MAP<T2>* color_map=color_maps(current_color_map);
    for(int i=start_index.x,i_active_cells=1;i<=end_index.x;i++,i_active_cells++) for(int j=start_index.y,j_active_cells=1;j<=end_index.y;j++,j_active_cells++) if(!active_cells || (*active_cells)(i_active_cells,j_active_cells)){
        opengl_points->points(index)=grid.X(i,j);
        opengl_points->Set_Point_Color(index,color_map->Lookup(Pre_Map_Value(values(i,j))));
        index++;}
    opengl_points->points.Resize(index-1);
}

template<class T,class T2> void OPENGL_SCALAR_FIELD_2D<T,T2>::
Update_Contour_Curves()
{
    // not supported by these types
    PHYSBAM_WARNING(std::string("Dual-contouring is not supported for scalar fields of type ")+typeid(T).name());
    contour_curves.Delete_Pointers_And_Clean_Memory();
}

template<> void OPENGL_SCALAR_FIELD_2D<float,float>::
Update_Contour_Curves()
{
    LEVELSET_2D<GRID<TV> > scalar_field_as_levelset(grid,values);
    OPENGL_COLOR_MAP<float>* color_map=color_maps(current_color_map);
    for(int i=1;i<=contour_curves.m;i++) delete contour_curves(i);
    if(!values.counts.Contains(0)){
        contour_curves.Resize(contour_values.m);
        for(int i=1;i<=contour_values.m;i++){
            contour_curves(i)=new OPENGL_SEGMENTED_CURVE_2D<float>(*DUALCONTOUR_2D<float>::Create_Segmented_Curve_From_Levelset(scalar_field_as_levelset,contour_values(i),false),color_map->Lookup(Pre_Map_Value(contour_values(i))));}}
    else contour_curves.Clean_Memory();
}

#ifndef COMPILE_WITHOUT_DOUBLE_SUPPORT
template<> void OPENGL_SCALAR_FIELD_2D<double,double>::
Update_Contour_Curves()
{
    LEVELSET_2D<GRID<TV> > scalar_field_as_levelset(grid,values);
    OPENGL_COLOR_MAP<double>* color_map=color_maps(current_color_map);
    for(int i=1;i<=contour_curves.m;i++) delete contour_curves(i);
    if(!values.counts.Contains(0)){
        contour_curves.Resize(contour_values.m);
        for(int i=1;i<=contour_values.m;i++){
            contour_curves(i)=new OPENGL_SEGMENTED_CURVE_2D<double>(*DUALCONTOUR_2D<double>::Create_Segmented_Curve_From_Levelset(scalar_field_as_levelset,contour_values(i),false),color_map->Lookup(Pre_Map_Value(contour_values(i))));}}
    else contour_curves.Clean_Memory();
}
#endif

//#####################################################################
// Specialization for bool scalars: only draw point if value is "true"
//#####################################################################
template<> void OPENGL_SCALAR_FIELD_2D<float,bool>::
Update_Points(const VECTOR<int,2>& start_index,const VECTOR<int,2>& end_index)
{
    PHYSBAM_ASSERT(opengl_points);
    OPENGL_COLOR_MAP<bool>* color_map=color_maps(current_color_map);
    opengl_points->color=color_map->Lookup(true);
    opengl_points->points.Resize((end_index.x-start_index.x+1)*(end_index.y-start_index.y+1));
    int index=1;
    for(int i=start_index.x;i<=end_index.x;i++) for(int j=start_index.y;j<=end_index.y;j++)
        if(values(i,j)) opengl_points->points(index++)=grid.X(i,j);
    opengl_points->points.Resize(index-1);
}

template<> void OPENGL_SCALAR_FIELD_2D<double,bool>::
Update_Points(const VECTOR<int,2>& start_index,const VECTOR<int,2>& end_index)
{
    PHYSBAM_ASSERT(opengl_points);
    OPENGL_COLOR_MAP<bool>* color_map=color_maps(current_color_map);
    opengl_points->color=color_map->Lookup(true);
    opengl_points->points.Resize((end_index.x-start_index.x+1)*(end_index.y-start_index.y+1));
    int index=1;
    for(int i=start_index.x;i<=end_index.x;i++) for(int j=start_index.y;j<=end_index.y;j++)
        if(values(i,j)) opengl_points->points(index++)=grid.X(i,j);
    opengl_points->points.Resize(index-1);
}

template<> void OPENGL_SCALAR_FIELD_2D<float,bool>::
Set_Uniform_Contour_Values(const bool min_value,const bool max_value,const bool increment)
{}

template<> void OPENGL_SCALAR_FIELD_2D<double,bool>::
Set_Uniform_Contour_Values(const bool min_value,const bool max_value,const bool increment)
{}

template<class T,class T2> void OPENGL_SCALAR_FIELD_2D<T,T2>::
Toggle_Draw_Mode()
{
    DRAW_MODE new_draw_mode=(DRAW_MODE)(((int)draw_mode+1)%3);
    Set_Draw_Mode(new_draw_mode);
}

template<class T,class T2> void OPENGL_SCALAR_FIELD_2D<T,T2>::
Toggle_Smooth_Texture()
{
    if(draw_mode==DRAW_TEXTURE && opengl_textured_rect && opengl_textured_rect->texture)
        opengl_textured_rect->texture->Toggle_Smooth_Shading();
}

template<class T,class T2> void OPENGL_SCALAR_FIELD_2D<T,T2>::
Toggle_Draw_Ghost_Values()
{
    draw_ghost_values=!draw_ghost_values;
    Update();
}

template<class T,class T2> void OPENGL_SCALAR_FIELD_2D<T,T2>::
Toggle_Color_Map()
{
    current_color_map=current_color_map%color_maps.m+1;
    Update();
}

template class OPENGL_SCALAR_FIELD_2D<float,int>;
template class OPENGL_SCALAR_FIELD_2D<float,bool>;
template class OPENGL_SCALAR_FIELD_2D<float,float>;
#ifndef COMPILE_WITHOUT_DOUBLE_SUPPORT
template class OPENGL_SCALAR_FIELD_2D<double,int>;
template class OPENGL_SCALAR_FIELD_2D<double,bool>;
template class OPENGL_SCALAR_FIELD_2D<double,double>;
#endif
}
