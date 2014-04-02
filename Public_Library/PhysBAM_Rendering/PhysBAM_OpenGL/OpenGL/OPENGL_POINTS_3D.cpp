//#####################################################################
// Copyright 2004-2008, Eran Guendelman, Craig Schroeder, Tamar Shinar.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#include <PhysBAM_Tools/Arrays/INDIRECT_ARRAY.h>
#include <PhysBAM_Rendering/PhysBAM_OpenGL/OpenGL/OPENGL_POINTS_3D.h>
using namespace PhysBAM;
//#####################################################################
// Constructor
//#####################################################################
template<class T,class T_ARRAY> OPENGL_POINTS_3D<T,T_ARRAY>::
OPENGL_POINTS_3D(T_ARRAY& points_input,const OPENGL_COLOR& color_input,const T point_size_input)
    :points(points_input),color(color_input),point_size(point_size_input),draw_point_numbers(false),point_colors(0),point_ids(0),draw_mask(0)
{}
//#####################################################################
// Destructor
//#####################################################################
template<class T,class T_ARRAY> OPENGL_POINTS_3D<T,T_ARRAY>::
~OPENGL_POINTS_3D()
{
    delete point_colors;
    delete point_ids;
    delete draw_mask;
}
//#####################################################################
// Function Bounding_Box
//#####################################################################
template<class T,class T_ARRAY> RANGE<VECTOR<float,3> > OPENGL_POINTS_3D<T,T_ARRAY>::
Bounding_Box() const
{
    if(!points.Size()) return RANGE<VECTOR<float,3> >::Empty_Box();
    return World_Space_Box(RANGE<VECTOR<T,3> >::Bounding_Box(points));
}
//#####################################################################
// Function Display
//#####################################################################
template<class T,class T_ARRAY> void OPENGL_POINTS_3D<T,T_ARRAY>::
Display(const int in_color) const
{
    glMatrixMode(GL_MODELVIEW);
    glPushMatrix();
    Send_Transform_To_GL_Pipeline();

    glPushAttrib(GL_ENABLE_BIT | GL_POINT_BIT);
    glPointSize(point_size);
    glDisable(GL_LIGHTING);

    GLint mode;
    glGetIntegerv(GL_RENDER_MODE,&mode);

    OPENGL_COLOR current_color;
    glGetFloatv(GL_CURRENT_COLOR,current_color.rgba);

    if(mode==GL_SELECT){
        glPushName(0);
        color.Send_To_GL_Pipeline();
        for(int i=1;i<=points.Size();i++)if(!draw_mask || (*draw_mask)(i)){
            glLoadName(i);
            OpenGL_Begin(GL_POINTS);
            OpenGL_Vertex(points(i));
            OpenGL_End();}
        glPopName();}
    else{
        color.Send_To_GL_Pipeline(); // set color outside OpenGL_Begin/OpenGL_End to work around apparent driver bug
        OpenGL_Begin(GL_POINTS);
        if(draw_mask){
            if(point_colors){for(int i=1;i<=points.Size();i++)if((*draw_mask)(i)){(*point_colors)(i).Send_To_GL_Pipeline();OpenGL_Vertex(points(i));}}
            else{color.Send_To_GL_Pipeline();for(int i=1;i<=points.Size();i++)if((*draw_mask)(i)) OpenGL_Vertex(points(i));}}
        else if(point_colors) for(int i=1;i<=points.Size();i++){(*point_colors)(i).Send_To_GL_Pipeline();OpenGL_Vertex(points(i));}
        else{color.Send_To_GL_Pipeline();for(int i=1;i<=points.Size();i++)OpenGL_Vertex(points(i));}
        OpenGL_End();

        if(draw_point_numbers)
            for(int i=1;i<=points.Size();i++)if(!draw_mask || (*draw_mask)(i)){
                OPENGL_COLOR label_color=(point_colors)?((*point_colors)(i)*0.8):(color*0.8);
                label_color.Send_To_GL_Pipeline();
                OpenGL_String(points(i),point_ids?STRING_UTILITIES::string_sprintf("%d [id=%d] [%g %g]",i,(*point_ids)(i),points(i).x,points(i).y):STRING_UTILITIES::string_sprintf("%d",i));}}

    current_color.Send_To_GL_Pipeline();

    glPopAttrib();
    glPopMatrix();
}
//#####################################################################
// Function Set_Point_Color
//#####################################################################
template<class T,class T_ARRAY> void OPENGL_POINTS_3D<T,T_ARRAY>::
Set_Point_Color(int index,const OPENGL_COLOR &point_color)
{
    Store_Point_Colors(true);
    (*point_colors)(index)=point_color;
}
//#####################################################################
// Function Set_Point_Colors
//#####################################################################
template<class T,class T_ARRAY> void OPENGL_POINTS_3D<T,T_ARRAY>::
Set_Point_Colors(const ARRAY<int> &indices, const OPENGL_COLOR &point_color)
{
    Store_Point_Colors(true);
    for (int i=1;i<=indices.m;i++) (*point_colors)(indices(i))=point_color;
}
//#####################################################################
// Function Reset_Point_Colors
//#####################################################################
template<class T,class T_ARRAY> void OPENGL_POINTS_3D<T,T_ARRAY>::
Reset_Point_Colors()
{
    if(point_colors) ARRAYS_COMPUTATIONS::Fill(*point_colors,color);
}    
//#####################################################################
// Function Store_Point_Colors
//#####################################################################
template<class T,class T_ARRAY> void OPENGL_POINTS_3D<T,T_ARRAY>::
Store_Point_Colors(const bool store_point_colors)
{
    if(store_point_colors){
        if(!point_colors){point_colors=new ARRAY<OPENGL_COLOR>(points.Size(),false);ARRAYS_COMPUTATIONS::Fill(*point_colors,color);}
        else if(point_colors->m!=points.Size()){point_colors->Resize(points.Size());ARRAYS_COMPUTATIONS::Fill(*point_colors,color);}}
    else{delete point_colors;point_colors=0;}
}
//#####################################################################
// Function Store_Point_Ids
//#####################################################################
template<class T,class T_ARRAY> void OPENGL_POINTS_3D<T,T_ARRAY>::
Store_Point_Ids(bool store_ids)
{
    if(store_ids){if(!point_ids) point_ids=new ARRAY<int>();point_ids->Resize(points.Size());}
    else{delete point_ids;point_ids=0;}
}
//#####################################################################
// Function Store_Draw_Mask
//#####################################################################
template<class T,class T_ARRAY> void OPENGL_POINTS_3D<T,T_ARRAY>::
Store_Draw_Mask(bool store_draw_mask)
{
    if(store_draw_mask){if(!draw_mask) draw_mask=new ARRAY<bool>();draw_mask->Resize(points.Size());}
    else{delete draw_mask;draw_mask=0;}
}
//#####################################################################
// Function Filter_By_Slice_Planes
//#####################################################################
template<class T,class T_ARRAY> void OPENGL_POINTS_3D<T,T_ARRAY>::
Filter_By_Slice_Planes(const PLANE<T> &leftplane, const PLANE<T> &rightplane)
{
    Store_Draw_Mask(true);
    for(int i=1;i<=points.Size();i++) (*draw_mask)(i)=leftplane.Lazy_Inside(points(i)) && rightplane.Lazy_Inside(points(i));
}
//#####################################################################
// Function Clear_Filter
//#####################################################################
template<class T,class T_ARRAY> void OPENGL_POINTS_3D<T,T_ARRAY>::
Clear_Filter()
{
    Store_Draw_Mask(false);
}
//#####################################################################
// Function Select_Point
//#####################################################################
template<class T,class T_ARRAY> void OPENGL_POINTS_3D<T,T_ARRAY>::
Select_Point(int index)
{
    Set_Point_Color(index,OPENGL_COLOR::Yellow());
}
//#####################################################################
// Function Select_Points
//#####################################################################
template<class T,class T_ARRAY> void OPENGL_POINTS_3D<T,T_ARRAY>::
Select_Points(const ARRAY<int> &indices)
{
    Set_Point_Colors(indices,OPENGL_COLOR::Yellow());
}
//#####################################################################
// Function Clear_Selection
//#####################################################################
template<class T,class T_ARRAY> void OPENGL_POINTS_3D<T,T_ARRAY>::
Clear_Selection()
{
    Store_Point_Colors(false);
}
//#####################################################################
// Function Get_Selection
//#####################################################################
template<class T,class T_ARRAY> OPENGL_SELECTION *OPENGL_POINTS_3D<T,T_ARRAY>::
Get_Selection(GLuint *buffer,int buffer_size)
{
    if(buffer_size==1){
        OPENGL_SELECTION_POINTS_3D<T> *selection=new OPENGL_SELECTION_POINTS_3D<T>(this);
        selection->index=buffer[0];
        if(point_ids){selection->has_id=true;selection->id=(*point_ids)(buffer[0]);}else selection->has_id=false;
        return selection;}
    return 0;
}
//#####################################################################
// Function Highlight_Selection
//#####################################################################
template<class T,class T_ARRAY> void OPENGL_POINTS_3D<T,T_ARRAY>::
Highlight_Selection(OPENGL_SELECTION *selection)
{
    if (selection->type!=OPENGL_SELECTION::POINTS_3D) return;
    OPENGL_SELECTION_POINTS_3D<T> *real_selection=(OPENGL_SELECTION_POINTS_3D<T>*)selection;
    Select_Point(real_selection->index);
}
//#####################################################################
// Function Clear_Highlight
//#####################################################################
template<class T,class T_ARRAY> void OPENGL_POINTS_3D<T,T_ARRAY>::
Clear_Highlight()
{
    Clear_Selection();
}
//#####################################################################
// Function Bounding_Box
//#####################################################################
template<class T> RANGE<VECTOR<float,3> > OPENGL_SELECTION_POINTS_3D<T>::
Bounding_Box() const
{
    PHYSBAM_ASSERT(object);
    if(OPENGL_POINTS_3D<T,ARRAY<VECTOR<T,3> > >* opengl_points=dynamic_cast<OPENGL_POINTS_3D<T,ARRAY<VECTOR<T,3> > >*>(object))
        return object->World_Space_Box(RANGE<VECTOR<float,3> >(VECTOR<float,3>(opengl_points->points(index))));
    else if(OPENGL_POINTS_3D<T,INDIRECT_ARRAY<ARRAY<VECTOR<T,3> > > >* opengl_points=dynamic_cast<OPENGL_POINTS_3D<T,INDIRECT_ARRAY<ARRAY<VECTOR<T,3> > > >*>(object))
        return object->World_Space_Box(RANGE<VECTOR<float,3> >(VECTOR<float,3>(opengl_points->points(index))));
    else PHYSBAM_NOT_IMPLEMENTED();
}
//#####################################################################
// Function Print_Selection_Info
//#####################################################################
template<class T,class T_ARRAY> void OPENGL_POINTS_3D<T,T_ARRAY>::
Print_Selection_Info(std::ostream &output_stream,OPENGL_SELECTION *selection) const
{
    if(selection->type!=OPENGL_SELECTION::POINTS_3D) return;
    int particle_index=Particle_Index(dynamic_cast<OPENGL_SELECTION_POINTS_3D<T>*>(selection)->index);
    output_stream<<"Free particle "<<particle_index<<std::endl;
}
//#####################################################################
template class OPENGL_POINTS_3D<float,ARRAY<VECTOR<float,3> > >;
template class OPENGL_POINTS_3D<float,INDIRECT_ARRAY<ARRAY<VECTOR<float,3> > > >;
template class OPENGL_POINTS_3D<float,INDIRECT_ARRAY<ARRAY_VIEW<VECTOR<float,3> > > >;
#ifndef COMPILE_WITHOUT_DOUBLE_SUPPORT
template class OPENGL_POINTS_3D<double,ARRAY<VECTOR<double,3> > >;
template class OPENGL_POINTS_3D<double,INDIRECT_ARRAY<ARRAY<VECTOR<double,3> > > >;
template class OPENGL_POINTS_3D<double,INDIRECT_ARRAY<ARRAY_VIEW<VECTOR<double,3> > > >;
#endif
