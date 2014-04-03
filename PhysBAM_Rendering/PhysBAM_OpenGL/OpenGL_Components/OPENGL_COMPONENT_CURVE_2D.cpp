//#####################################################################
// Copyright 2004, Eran Guendelman.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#include <PhysBAM_Tools/Math_Tools/sign.h>
#include <PhysBAM_Tools/Read_Write/Grids_Uniform_Arrays/READ_WRITE_ARRAYS.h>
#include <PhysBAM_Tools/Read_Write/Utilities/FILE_UTILITIES.h>
#include <PhysBAM_Rendering/PhysBAM_OpenGL/OpenGL/OPENGL_PRIMITIVES.h>
#include <PhysBAM_Rendering/PhysBAM_OpenGL/OpenGL/OPENGL_SHAPES.h>
#include <PhysBAM_Rendering/PhysBAM_OpenGL/OpenGL_Components/OPENGL_COMPONENT_CURVE_2D.h>
using namespace PhysBAM;
//#####################################################################
// Constructor
//#####################################################################
template<class T,class RW> OPENGL_COMPONENT_CURVE_2D<T,RW>::
OPENGL_COMPONENT_CURVE_2D(const std::string& x_filename_input,const std::string& u_filename_input,const std::string& flux_filename_input,const T domain_xmin_input,const T domain_xmax_input)
    :OPENGL_COMPONENT("Curve 2D"),x_filename(x_filename_input),u_filename(u_filename_input),flux_filename(flux_filename_input),domain_xmin(domain_xmin_input),domain_xmax(domain_xmax_input),
    flux_scale(1000),draw_flux(false),draw_du(false),draw_piecewise_constant(false),ghost_nodes(1),
    frame_loaded(-1),valid(false),current_selection(0)
{
    is_animation=true;
    Reinitialize();
}
//#####################################################################
// Destructor
//#####################################################################
template<class T,class RW> OPENGL_COMPONENT_CURVE_2D<T,RW>::
~OPENGL_COMPONENT_CURVE_2D()
{
}
//#####################################################################
// Valid_Frame
//#####################################################################
template<class T,class RW> bool OPENGL_COMPONENT_CURVE_2D<T,RW>::
Valid_Frame(int frame_input) const
{
    return FILE_UTILITIES::Frame_File_Exists(x_filename,frame_input) && FILE_UTILITIES::Frame_File_Exists(u_filename,frame_input);
}
//#####################################################################
// Set_Frame
//#####################################################################
template<class T,class RW> void OPENGL_COMPONENT_CURVE_2D<T,RW>::
Set_Frame(int frame_input)
{
    OPENGL_COMPONENT::Set_Frame(frame_input);
    Reinitialize();
}
//#####################################################################
// Set_Draw
//#####################################################################
template<class T,class RW> void OPENGL_COMPONENT_CURVE_2D<T,RW>::
Set_Draw(bool draw_input)
{
    OPENGL_COMPONENT::Set_Draw(draw_input);
    Reinitialize();
}
//#####################################################################
// Display
//#####################################################################
template<class T,class RW> void OPENGL_COMPONENT_CURVE_2D<T,RW>::
Display(const int in_color) const
{
    glPushAttrib(GL_ENABLE_BIT | GL_POINT_BIT | GL_LINE_BIT);
    glDisable(GL_LIGHTING);
    glShadeModel(GL_FLAT);
    OPENGL_COLOR positive_segment_color=OPENGL_COLOR::Green(0.5);
    OPENGL_COLOR negative_segment_color=OPENGL_COLOR::Red(0.5);
    OPENGL_COLOR positive_point_color=OPENGL_COLOR::Green();
    OPENGL_COLOR negative_point_color=OPENGL_COLOR::Red();
    OPENGL_COLOR flux_color=OPENGL_COLOR::Magenta();
    OPENGL_COLOR du_color=OPENGL_COLOR::Yellow();

    GLint mode;
    glGetIntegerv(GL_RENDER_MODE, &mode);
    if(mode == GL_SELECT){
        glPointSize(OPENGL_PREFERENCES::selection_point_size);
        glPushName(1); // vertices
        glPushName(0);
        for(int i=x.domain.min_corner.x;i<=x.domain.max_corner.x;i++){
            glLoadName(i);
            OpenGL_Begin(GL_POINTS);
            OpenGL_Vertex(VECTOR<T,2>(x(i),u(i)));
            OpenGL_End();}
        glPopName();

        glLoadName(2); // segments
        glPushName(0);
        for(int i=x.domain.min_corner.x;i<=x.domain.max_corner.x-1;i++){
            glLoadName(i);
            OpenGL_Begin(GL_LINES);
            OpenGL_Line(VECTOR<T,2>(x(i),u(i)),VECTOR<T,2>(x(i+1),u(i+1)));
            OpenGL_End();}
        glPopName();
        glPopName();
    }
    else{
        OPENGL_COLOR::Gray(0.7).Send_To_GL_Pipeline();
        OpenGL_Begin(GL_LINES);
        OpenGL_Line(VECTOR<T,2>(domain_xmin,-1000),VECTOR<T,2>(domain_xmin,1000));
        OpenGL_Line(VECTOR<T,2>(domain_xmax,-1000),VECTOR<T,2>(domain_xmax,1000));
        OpenGL_End();
        OPENGL_COLOR::Gray(0.3).Send_To_GL_Pipeline();
        OpenGL_Begin(GL_LINES);
        OpenGL_Line(VECTOR<T,2>(domain_xmin,0),VECTOR<T,2>(domain_xmax,0));
        OpenGL_End();

        if(draw_piecewise_constant){
            OpenGL_Begin(GL_LINES);
            for(int i=x.domain.min_corner.x+1;i<=x.domain.max_corner.x-1;i++){
                T dx=(T).5*(x(i+1)-x(i-1));
                if(sign(dx)>0) positive_segment_color.Send_To_GL_Pipeline(); else negative_segment_color.Send_To_GL_Pipeline();
                OpenGL_Line(VECTOR<T,2>((T).5*(x(i-1)+x(i)),u(i)),VECTOR<T,2>((T).5*(x(i)+x(i+1)),u(i)));}
            OpenGL_End();

            glPointSize(5);
            OpenGL_Begin(GL_POINTS);
            for(int i=x.domain.min_corner.x+1;i<=x.domain.max_corner.x-1;i++){
                T dx=(T).5*(x(i+1)-x(i-1));
                if(sign(dx)>0) positive_point_color.Send_To_GL_Pipeline(); else negative_point_color.Send_To_GL_Pipeline();
                OpenGL_Vertex(VECTOR<T,2>(x(i),u(i)));}
            OpenGL_End();

            if(draw_flux && flux.counts.x){
                flux_color.Send_To_GL_Pipeline();
                OpenGL_Begin(GL_LINES);
                for(int i=x.domain.min_corner.x;i<=x.domain.max_corner.x-1;i++){
                    T dx=x(i+1)-x(i);
                    VECTOR<T,2> start_point((T).5*(x(i)+x(i+1)),(T).5*(u(i)+u(i+1)));
                    OPENGL_SHAPES::Draw_Arrow(start_point,start_point+VECTOR<T,2>(sign(dx)*flux_scale*flux(i),0));}
                OpenGL_End();

                // draw flux on segment
                glPushAttrib(GL_LINE_BIT);
                glLineWidth(2);
                OpenGL_Begin(GL_LINES);
                T arrowhead_size=.3,sin_angle=sin(0.5),cos_angle=cos(0.5);
                for(int i=x.domain.min_corner.x;i<=x.domain.max_corner.x-1;i++){
                    T dx=x(i+1)-x(i);T alpha=sign(flux(i)*dx)>0?(T).9:(T).1;
                    VECTOR<T,2> startpt((T).5*(x(i)+x(i+1)),(T).5*(u(i)+u(i+1)));
                    VECTOR<T,2> endpt(x(i)+alpha*(x(i+1)-x(i)),u(i)+alpha*(u(i+1)-u(i)));
                    VECTOR<T,2> direction=endpt-startpt,p=endpt-(cos_angle*arrowhead_size)*direction,dp=(sin_angle*arrowhead_size)*direction.Rotate_Clockwise_90();
                    OpenGL_Line(p+dp,endpt);
                    OpenGL_Line(p-dp,endpt);}
                OpenGL_End();
                glPopAttrib();
            }
        }
        else{
            OpenGL_Begin(GL_LINE_STRIP);
            for(int i=x.domain.min_corner.x;i<=x.domain.max_corner.x;i++){
                if(i>x.domain.min_corner.x && x(i)<x(i-1)) negative_segment_color.Send_To_GL_Pipeline();
                else positive_segment_color.Send_To_GL_Pipeline();
                OpenGL_Vertex(VECTOR<T,2>(x(i),u(i)));}
            OpenGL_End();

            glPointSize(5);
            OpenGL_Begin(GL_POINTS);
            for(int i=x.domain.min_corner.x;i<=x.domain.max_corner.x;i++){
                T dx=(i==x.domain.min_corner.x)?x(i+1)-x(i):(i==x.domain.max_corner.x)?x(i)-x(i-1):(T).5*(x(i+1)-x(i-1));
                if(sign(dx)>0) positive_point_color.Send_To_GL_Pipeline();
                else negative_point_color.Send_To_GL_Pipeline();
                OpenGL_Vertex(VECTOR<T,2>(x(i),u(i)));
            }
            OpenGL_End();

            if(draw_flux && flux.counts.x){
                flux_color.Send_To_GL_Pipeline();
                OpenGL_Begin(GL_LINES);
                for(int i=x.domain.min_corner.x;i<=x.domain.max_corner.x-1;i++){
                    T dx=x(i+1)-x(i);
                    VECTOR<T,2> start_point((T).5*(x(i)+x(i+1)),(T).5*(u(i)+u(i+1)));
                    OPENGL_SHAPES::Draw_Arrow(start_point,start_point+VECTOR<T,2>(sign(dx)*flux_scale*flux(i),0));}
                OpenGL_End();

                // draw flux on segment
                glPushAttrib(GL_LINE_BIT);
                glLineWidth(2);
                OpenGL_Begin(GL_LINES);
                T arrowhead_size=.3,sin_angle=sin(0.5),cos_angle=cos(0.5);
                for(int i=x.domain.min_corner.x;i<=x.domain.max_corner.x-1;i++){
                    T dx=x(i+1)-x(i);T alpha=sign(flux(i)*dx)>0?(T).9:(T).1;
                    VECTOR<T,2> startpt((T).5*(x(i)+x(i+1)),(T).5*(u(i)+u(i+1)));
                    VECTOR<T,2> endpt(x(i)+alpha*(x(i+1)-x(i)),u(i)+alpha*(u(i+1)-u(i)));
                    VECTOR<T,2> direction=endpt-startpt,p=endpt-(cos_angle*arrowhead_size)*direction,dp=(sin_angle*arrowhead_size)*direction.Rotate_Clockwise_90();
                    OpenGL_Line(p+dp,endpt);
                    OpenGL_Line(p-dp,endpt);}
                OpenGL_End();
                glPopAttrib();
            }

            if(draw_du && flux.counts.x){
                du_color.Send_To_GL_Pipeline();
                OpenGL_Begin(GL_LINES);
                for(int i=x.domain.min_corner.x+1;i<=x.domain.max_corner.x-1;i++){
                    T dx=x(i+1)-x(i),du=(flux(i-1)-flux(i))/dx;
                    VECTOR<T,2> start_point(x(i),u(i));
                    OPENGL_SHAPES::Draw_Arrow(start_point,start_point+VECTOR<T,2>(0,flux_scale*du));}
                OpenGL_End();
            }
        }

        if(current_selection){
            if(current_selection->type == OPENGL_SELECTION::COMPONENT_CURVE_VERTEX_2D){
                int index=((OPENGL_SELECTION_COMPONENT_CURVE_VERTEX_2D<T,RW> *)current_selection)->index;
                OPENGL_SELECTION::Draw_Highlighted_Vertex(VECTOR<T,3>(x(index),u(index),0));} 
            else if(current_selection->type == OPENGL_SELECTION::COMPONENT_CURVE_SEGMENT_2D){
                int index=((OPENGL_SELECTION_COMPONENT_CURVE_SEGMENT_2D<T,RW> *)current_selection)->index;
                OPENGL_SELECTION::Draw_Highlighted_Segment(VECTOR<T,3>(x(index),u(index),0),VECTOR<T,3>(x(index+1),u(index+1),0));}
        }
    }
    glPopAttrib();
}
//#####################################################################
// Bounding_Box
//#####################################################################
template<class T,class RW> RANGE<VECTOR<float,3> > OPENGL_COMPONENT_CURVE_2D<T,RW>::
Bounding_Box() const
{
    if(valid&&draw){
        T xmin=x.Min(),xmax=x.Max();
        T umin=u.Min(),umax=u.Max();
        return World_Space_Box(RANGE<VECTOR<float,3> >(xmin,xmax,umin,umax,0,0));}
    else return RANGE<VECTOR<float,3> >::Centered_Box();
}
//#####################################################################
// Function Get_Selection
//#####################################################################
template<class T,class RW> OPENGL_SELECTION *OPENGL_COMPONENT_CURVE_2D<T,RW>::
Get_Selection(GLuint *buffer,int buffer_size)
{
    OPENGL_SELECTION *selection = 0;
    if(buffer_size == 2){
        if(buffer[0] == 1) selection=Get_Vertex_Selection(buffer[1]);
        else if(buffer[0] == 2) selection=Get_Segment_Selection(buffer[1]);
    }
    return selection;
}
//#####################################################################
// Function Highlight_Selection
//#####################################################################
template<class T,class RW> void OPENGL_COMPONENT_CURVE_2D<T,RW>::
Highlight_Selection(OPENGL_SELECTION *selection)
{
    delete current_selection;current_selection=0;
    // Make a copy of selection
    if(selection->type == OPENGL_SELECTION::COMPONENT_CURVE_VERTEX_2D)
        current_selection=new OPENGL_SELECTION_COMPONENT_CURVE_VERTEX_2D<T,RW>(this,((OPENGL_SELECTION_COMPONENT_CURVE_VERTEX_2D<T,RW> *)selection)->index);
    else if(selection->type == OPENGL_SELECTION::COMPONENT_CURVE_SEGMENT_2D)
        current_selection=new OPENGL_SELECTION_COMPONENT_CURVE_SEGMENT_2D<T,RW>(this,((OPENGL_SELECTION_COMPONENT_CURVE_SEGMENT_2D<T,RW> *)selection)->index);
}
//#####################################################################
// Function Print_Selection_Info
//#####################################################################
template<class T,class RW> void OPENGL_COMPONENT_CURVE_2D<T,RW>::
Print_Selection_Info(std::ostream &output_stream, OPENGL_SELECTION *selection) const
{
    if (selection->type == OPENGL_SELECTION::COMPONENT_CURVE_VERTEX_2D) {
        int index=((OPENGL_SELECTION_COMPONENT_CURVE_VERTEX_2D<T,RW> *)selection)->index;
        output_stream << "vertex " << index << std::endl;
        output_stream << "position " << VECTOR<T,2>(x(index),u(index)) << std::endl;
        output_stream << "dx=";
        if(index>x.domain.min_corner.x&&index<x.domain.max_corner.x){output_stream<<(T).5*(x(index+1)-x(index-1));}else{output_stream<<"n/a";}
        output_stream << " left=";
        if(index>x.domain.min_corner.x){output_stream<<x(index)-x(index-1);}else{output_stream<<"n/a";}
        output_stream << " right=";
        if(index<x.domain.max_corner.x){output_stream<<x(index+1)-x(index);}else{output_stream<<"n/a";}
        output_stream << std::endl;
    }
    else if (selection->type == OPENGL_SELECTION::COMPONENT_CURVE_SEGMENT_2D) {
        int index=((OPENGL_SELECTION_COMPONENT_CURVE_SEGMENT_2D<T,RW> *)selection)->index;
        output_stream << "segment " << index << std::endl;
        output_stream << "dx=" << x(index+1)-x(index) << std::endl;
        if(index>x.domain.min_corner.x) output_stream << "dx left node=" << (T).5*(x(index+1)-x(index-1)) << std::endl;
        if(index<x.domain.max_corner.x) output_stream << "dx right node=" << (T).5*(x(index+2)-x(index)) << std::endl;
        if(flux.counts.x) output_stream << "flux=" << flux(index) << std::endl;
    }
}
//#####################################################################
// Function Clear_Highlight
//#####################################################################
template<class T,class RW> void OPENGL_COMPONENT_CURVE_2D<T,RW>::
Clear_Highlight()
{
    delete current_selection;current_selection=0;
}
//#####################################################################
// Function Get_Vertex_Selection
//#####################################################################
template<class T,class RW> OPENGL_SELECTION *OPENGL_COMPONENT_CURVE_2D<T,RW>::
Get_Vertex_Selection(int index)
{
    return new OPENGL_SELECTION_COMPONENT_CURVE_VERTEX_2D<T,RW>(this,index);
}
//#####################################################################
// Function Get_Segment_Selection
//#####################################################################
template<class T,class RW> OPENGL_SELECTION *OPENGL_COMPONENT_CURVE_2D<T,RW>::
Get_Segment_Selection(int index)
{
    return new OPENGL_SELECTION_COMPONENT_CURVE_SEGMENT_2D<T,RW>(this,index);
}
//#####################################################################
// Reinitialize
//#####################################################################
template<class T,class RW> void OPENGL_COMPONENT_CURVE_2D<T,RW>::
Reinitialize(bool force)
{
    if(draw){
        if(force || (is_animation && frame_loaded != frame) || (!is_animation && frame_loaded < 0)){
            valid=false;
            FILE_UTILITIES::Read_From_File<RW>(FILE_UTILITIES::Get_Frame_Filename(x_filename,frame),x);
            FILE_UTILITIES::Read_From_File<RW>(FILE_UTILITIES::Get_Frame_Filename(u_filename,frame),u);
            if(FILE_UTILITIES::Frame_File_Exists(flux_filename,frame)) FILE_UTILITIES::Read_From_File<RW>(FILE_UTILITIES::Get_Frame_Filename(flux_filename,frame),flux);
            else flux.Clean_Memory();
            valid=(x.counts.x && ARRAY<T,VECTOR<int,1> >::Equal_Dimensions(x,u));
            frame_loaded=frame;}}
}
//#####################################################################
// Area_Under_Curve
//#####################################################################
template<class T,class RW> T OPENGL_COMPONENT_CURVE_2D<T,RW>::
Area_Under_Curve() const
{
    T area=0;
    for(int i=0;i<=x.domain.max_corner.x-1-ghost_nodes;i++) area+=0.5*(u(i)+u(i+1))*(x(i+1)-x(i));
    return area;
}
//#####################################################################
// Callbacks
//#####################################################################
template<class T,class RW> void OPENGL_COMPONENT_CURVE_2D<T,RW>::
Increase_Vector_Size()
{
    flux_scale*=1.1;
}
template<class T,class RW> void OPENGL_COMPONENT_CURVE_2D<T,RW>::
Toggle_Draw_Flux()
{
    draw_flux=!draw_flux;
}
template<class T,class RW> void OPENGL_COMPONENT_CURVE_2D<T,RW>::
Toggle_Draw_Du()
{
    draw_du=!draw_du;
}
template<class T,class RW> void OPENGL_COMPONENT_CURVE_2D<T,RW>::
Toggle_Draw_Piecewise_Constant()
{
    draw_piecewise_constant=!draw_piecewise_constant;
}
//#####################################################################
// Decrease_Vector_Size
//#####################################################################
template<class T,class RW> void OPENGL_COMPONENT_CURVE_2D<T,RW>::
Decrease_Vector_Size()
{
    flux_scale*=1/1.1;
}
//#####################################################################
// Bounding_Box
//#####################################################################
template<class T,class RW>
RANGE<VECTOR<float,3> > OPENGL_SELECTION_COMPONENT_CURVE_VERTEX_2D<T,RW>::
Bounding_Box() const
{
    PHYSBAM_ASSERT(object);
    const OPENGL_COMPONENT_CURVE_2D<T,RW> &curve=*(OPENGL_COMPONENT_CURVE_2D<T,RW> *)object;
    RANGE<VECTOR<float,3> > box(VECTOR<float,3>(curve.x(index),curve.u(index),0));
    return object->World_Space_Box(box);
}
//#####################################################################
// Bounding_Box
//#####################################################################
template<class T,class RW>
RANGE<VECTOR<float,3> > OPENGL_SELECTION_COMPONENT_CURVE_SEGMENT_2D<T,RW>::
Bounding_Box() const
{
    PHYSBAM_ASSERT(object);
    const OPENGL_COMPONENT_CURVE_2D<T,RW> &curve=*(OPENGL_COMPONENT_CURVE_2D<T,RW> *)object;
    RANGE<VECTOR<float,3> > box=RANGE<VECTOR<float,3> >::Bounding_Box(VECTOR<float,3>(curve.x(index),curve.u(index),0),VECTOR<float,3>(curve.x(index+1),curve.u(index+1),0));
    return object->World_Space_Box(box);
}
template class OPENGL_COMPONENT_CURVE_2D<float,float>;
#ifndef COMPILE_WITHOUT_DOUBLE_SUPPORT
template class OPENGL_COMPONENT_CURVE_2D<double,double>;
#endif
