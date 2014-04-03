#ifndef COMPILE_WITHOUT_DYADIC_SUPPORT
//#####################################################################
// Copyright 2004-2005, Eran Guendelman, Frank Losasso, Andrew Selle.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class OPENGL_OCTREE_SLICE
//##################################################################### 
#include <PhysBAM_Rendering/PhysBAM_OpenGL/OpenGL/OPENGL_OCTREE_SLICE.h>
#include <PhysBAM_Rendering/PhysBAM_OpenGL/OpenGL/OPENGL_WORLD.h>
using namespace PhysBAM;
//#####################################################################
// Function Print_Slice_Info
//#####################################################################
void OPENGL_OCTREE_SLICE::
Print_Slice_Info(std::ostream& output_stream)
{
    const char* axis_name[]={"","x","y","z"};
    if(mode==CELL_SLICE){
        output_stream<<"Slice: "<<axis_name[axis]<<" from "<<position<<" to "<<+position+Get_Slice_Increment();
        RANGE<VECTOR<float,3> > domain=octree_grid->uniform_grid.domain;
        if(position+.5*Get_Slice_Increment()<domain.Minimum_Corner()[axis] || position+.5*Get_Slice_Increment()>domain.Maximum_Corner()[axis]) output_stream<<" (ghost)";
        output_stream<<std::endl;
    }
    if(mode==NODE_SLICE){
        output_stream<<"Slice: "<<axis_name[axis]<<"="<<position;
        RANGE<VECTOR<float,3> > domain=octree_grid->uniform_grid.domain;
        if(position<domain.Minimum_Corner()[axis] || position>domain.Maximum_Corner()[axis]) output_stream<<" (ghost)";
        output_stream<<std::endl;
    }
}
//#####################################################################
// Function Update_Clip_Planes
//#####################################################################
void OPENGL_OCTREE_SLICE::
Update_Clip_Planes()
{
    if (mode==NO_SLICE) {
        if (clip_plane_id1!=0) world.Remove_Clipping_Plane(clip_plane_id1); 
        if (clip_plane_id2!=0) world.Remove_Clipping_Plane(clip_plane_id2);
        clip_plane_id1=clip_plane_id2=0;
    }
    else if(mode==CELL_SLICE){
        PLANE<float> plane1(VECTOR<float,3>(0,0,0),VECTOR<float,3>(0,0,0)),plane2(VECTOR<float,3>(0,0,0),VECTOR<float,3>(0,0,0));
        float dx=Get_Slice_Increment();
        plane1.normal[axis]=1;plane1.x1[axis]=position-0.0001;
        plane2.normal[axis]=-1;plane2.x1[axis]=position+dx+0.0001;
        if (clip_plane_id1==0) clip_plane_id1=world.Add_Clipping_Plane(plane1);
        else world.Set_Clipping_Plane(clip_plane_id1,plane1);
        if (clip_plane_id2==0) clip_plane_id2=world.Add_Clipping_Plane(plane2);
        else world.Set_Clipping_Plane(clip_plane_id2,plane2);
    }
    else if(mode==NODE_SLICE){
        PLANE<float> plane1(VECTOR<float,3>(0,0,0),VECTOR<float,3>(0,0,0)),plane2(VECTOR<float,3>(0,0,0),VECTOR<float,3>(0,0,0));
        float dx=Get_Slice_Increment();
        plane1.normal[axis]=1;plane1.x1[axis]=position-dx/1.999;
        plane2.normal[axis]=-1;plane2.x1[axis]=position+dx/1.999;
        if (clip_plane_id1==0) clip_plane_id1=world.Add_Clipping_Plane(plane1);
        else world.Set_Clipping_Plane(clip_plane_id1,plane1);
        if (clip_plane_id2==0) clip_plane_id2=world.Add_Clipping_Plane(plane2);
        else world.Set_Clipping_Plane(clip_plane_id2,plane2);
    }
}
#endif
