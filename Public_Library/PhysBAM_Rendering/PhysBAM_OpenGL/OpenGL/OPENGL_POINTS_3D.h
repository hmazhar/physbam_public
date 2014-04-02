//#####################################################################
// Copyright 2004-2007, Eran Guendelman, Geoffrey Irving, Craig Schroeder, Tamar Shinar.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class OPENGL_POINTS_3D
//#####################################################################
#ifndef __OPENGL_POINTS_3D__
#define __OPENGL_POINTS_3D__

#include <PhysBAM_Tools/Arrays/ARRAY.h>
#include <PhysBAM_Tools/Arrays/INDIRECT_ARRAY.h>
#include <PhysBAM_Tools/Vectors/VECTOR_3D.h>
#include <PhysBAM_Geometry/Basic_Geometry/PLANE.h>
#include <PhysBAM_Rendering/PhysBAM_OpenGL/OpenGL/OPENGL_COLOR.h>
#include <PhysBAM_Rendering/PhysBAM_OpenGL/OpenGL/OPENGL_OBJECT.h>
#include <PhysBAM_Rendering/PhysBAM_OpenGL/OpenGL/OPENGL_SELECTION.h>
namespace PhysBAM{

template<class T,class T_ARRAY=ARRAY<VECTOR<T,3> > >
class OPENGL_POINTS_3D:public OPENGL_OBJECT
{
    typedef VECTOR<T,3> TV;
public:
    T_ARRAY& points;
    OPENGL_COLOR color;
    T point_size;
    bool draw_point_numbers; 
    ARRAY<OPENGL_COLOR>* point_colors;
    ARRAY<int>* point_ids;
    ARRAY<bool>* draw_mask;

    OPENGL_POINTS_3D(T_ARRAY& points_input,const OPENGL_COLOR& color_input=OPENGL_COLOR::White(),const T point_size=5);
    ~OPENGL_POINTS_3D();

    template<class T_PARTICLES> void Set_Points_From_Particles(const T_PARTICLES& particles,const bool keep_colors=true,const bool use_ids=true,const bool use_draw_mask=false)
    {points.Resize(particles.array_collection->Size());
    //Store_Point_Ids(use_ids && particles.store_id);
    Store_Point_Colors(false);
    if(use_draw_mask) Store_Draw_Mask(true);
    for(int i=1;i<=particles.array_collection->Size();i++){
        if(draw_mask) (*draw_mask)(i)=true;
        points(i)=particles.X(i);}}

    bool Use_Bounding_Box() const PHYSBAM_OVERRIDE {return points.Size()>0;}
    virtual int Particle_Index(const int index) const {return index;}
    virtual RANGE<VECTOR<float,3> > Bounding_Box() const PHYSBAM_OVERRIDE;
    void Display(const int in_color=1) const PHYSBAM_OVERRIDE;

    virtual OPENGL_SELECTION* Get_Selection(GLuint* buffer,int buffer_size);
    void Highlight_Selection(OPENGL_SELECTION* selection) PHYSBAM_OVERRIDE;
    void Clear_Highlight() PHYSBAM_OVERRIDE;
    void Print_Selection_Info(std::ostream &output_stream,OPENGL_SELECTION *selection) const PHYSBAM_OVERRIDE;

    void Store_Point_Colors(const bool store_point_colors=true);
    void Store_Point_Ids(bool store_ids=true);
    void Store_Draw_Mask(bool store_draw_mask=true);
    
    void Set_Point_Color(int index,const OPENGL_COLOR& point_color);
    void Set_Point_Colors(const ARRAY<int>& indices,const OPENGL_COLOR& point_color);
    void Reset_Point_Colors();

    void Filter_By_Slice_Planes(const PLANE<T>& leftplane,const PLANE<T>& rightplane);
    void Clear_Filter();

    void Select_Point(int index);
    void Select_Points(const ARRAY<int>& indices);
    void Clear_Selection();
};

template<class T>
class OPENGL_SELECTION_POINTS_3D:public OPENGL_SELECTION
{
public:
    int index;
    bool has_id;
    int id;

    OPENGL_SELECTION_POINTS_3D(OPENGL_OBJECT* object):OPENGL_SELECTION(OPENGL_SELECTION::POINTS_3D,object){}

    RANGE<VECTOR<float,3> > Bounding_Box() const PHYSBAM_OVERRIDE;
};

}
#endif
