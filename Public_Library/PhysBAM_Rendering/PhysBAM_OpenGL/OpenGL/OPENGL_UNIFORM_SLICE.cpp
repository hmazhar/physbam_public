//#####################################################################
// Copyright 2004, Eran Guendelman.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class OPENGL_UNIFORM_SLICE
//##################################################################### 
#include <PhysBAM_Rendering/PhysBAM_OpenGL/OpenGL/OPENGL_UNIFORM_SLICE.h>
#include <PhysBAM_Rendering/PhysBAM_OpenGL/OpenGL/OPENGL_WORLD.h>
using namespace PhysBAM;
//#####################################################################
// Function Print_Slice_Info
//#####################################################################
void OPENGL_UNIFORM_SLICE::
Print_Slice_Info(std::ostream& output_stream)
{
    const char* axis_name[]={"","x","y","z"};
    if(mode==CELL_SLICE){output_stream<<"Slice: cell="<<index<<", "<<axis_name[axis]<<"="<<grid.Center(index,index,index)[axis]<<std::endl;}
    else if(mode==NODE_SLICE){output_stream<<"Slice: node="<<index<<", "<<axis_name[axis]<<"="<<grid.Node(index,index,index)[axis]<<std::endl;}
}
//#####################################################################
// Function Update_Clip_Planes
//#####################################################################
void OPENGL_UNIFORM_SLICE::
Update_Clip_Planes()
{
    if (mode==NO_SLICE) {
        if (clip_plane_id1!=0) world.Remove_Clipping_Plane(clip_plane_id1); 
        if (clip_plane_id2!=0) world.Remove_Clipping_Plane(clip_plane_id2);
        clip_plane_id1=clip_plane_id2=0;
    }
    else {
        float pos=(mode==NODE_SLICE)?grid.Node(index,index,index)[axis]:grid.Center(index,index,index)[axis];
        PLANE<float> plane1(VECTOR<float,3>(0,0,0),VECTOR<float,3>(0,0,0)),plane2(VECTOR<float,3>(0,0,0),VECTOR<float,3>(0,0,0));
        plane1.normal[axis]=1;plane1.x1[axis]=pos-grid.dX[axis]/1.9;
        plane2.normal[axis]=-1;plane2.x1[axis]=pos+grid.dX[axis]/1.9;
        if (clip_plane_id1==0) clip_plane_id1=world.Add_Clipping_Plane(plane1);
        else world.Set_Clipping_Plane(clip_plane_id1,plane1);
        if (clip_plane_id2==0) clip_plane_id2=world.Add_Clipping_Plane(plane2);
        else world.Set_Clipping_Plane(clip_plane_id2,plane2);
    }
}
//#####################################################################
