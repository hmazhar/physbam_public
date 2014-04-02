//#####################################################################
// Copyright 2009, Jon Gretarsson, Nipun Kwatra.
// This file is part of PhysBAM whose distribution is governed by the license 
// contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class OPENGL_POINT_SIMPLICES_1D
//##################################################################### 
#ifndef __OPENGL_POINT_SIMPLICES_1D__
#define __OPENGL_POINT_SIMPLICES_1D__

#include <PhysBAM_Geometry/Topology_Based_Geometry/POINT_SIMPLICES_1D.h>
#include <PhysBAM_Rendering/PhysBAM_OpenGL/OpenGL/OPENGL_COLOR.h>
#include <PhysBAM_Rendering/PhysBAM_OpenGL/OpenGL/OPENGL_OBJECT.h>
#include <PhysBAM_Rendering/PhysBAM_OpenGL/OpenGL/OPENGL_SELECTION.h>
namespace PhysBAM{

class OPENGL_SELECTION;

template<class T>
class OPENGL_POINT_SIMPLICES_1D: public OPENGL_OBJECT
{
public:
    const POINT_SIMPLICES_1D<T>& simplices;
    OPENGL_COLOR color,color_gray;
    OPENGL_COLOR vertex_color,segment_color,vertex_position_color,velocity_color;
    bool draw_vertices,draw_vertex_positions;

    OPENGL_POINT_SIMPLICES_1D(const POINT_SIMPLICES_1D<T>& simplices_input,const OPENGL_COLOR &color_input=OPENGL_COLOR::Cyan());

//#####################################################################
    void Display(const int in_color=1) const PHYSBAM_OVERRIDE;
    virtual RANGE<VECTOR<float,3> > Bounding_Box() const PHYSBAM_OVERRIDE;
//#####################################################################
};
}
#endif
