//#####################################################################
// Copyright 2009.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#ifndef __TRIANGULATED_SURFACE_UPDATE_VERTEX_NORMALS__
#define __TRIANGULATED_SURFACE_UPDATE_VERTEX_NORMALS__

#include <PhysBAM_Tools/Arrays/ARRAY.h>
namespace PhysBAM{
template<class T> class TRIANGULATED_SURFACE;
template<class T> class TRIANGLE_3D;

namespace TOPOLOGY_BASED_GEOMETRY_COMPUTATIONS
{
template<class T> void
Update_Vertex_Normals(TRIANGULATED_SURFACE<T>& ts)
{  
    typedef VECTOR<T,3> TV;
    bool incident_elements_defined=ts.mesh.incident_elements!=0;if(!incident_elements_defined) ts.mesh.Initialize_Incident_Elements();
    bool triangle_list_defined=ts.triangle_list!=0;if(!triangle_list_defined) ts.Update_Triangle_List();

    if(ts.avoid_normal_interpolation_across_sharp_edges){
        delete ts.vertex_normals;ts.vertex_normals=0;
        if(!ts.face_vertex_normals) ts.face_vertex_normals=new ARRAY<VECTOR<TV,3> >(ts.mesh.elements.m);else ts.vertex_normals->Resize(ts.mesh.elements.m);
        TV face_normal,normal,zero_vector(0,0,0);
        for(int t=1;t<=ts.mesh.elements.m;t++){
            face_normal=ts.Face_Normal(t);
            for(int i=1;i<=3;i++){
                normal=zero_vector;int node=ts.mesh.elements(t)(i); 
                for(int k=1;k<=(*ts.mesh.incident_elements)(node).m;k++){
                    TRIANGLE_3D<T>& triangle=(*ts.triangle_list)((*ts.mesh.incident_elements)(node)(k));
                    TV local_normal=triangle.normal;
                    if((face_normal-local_normal).Magnitude_Squared() < ts.normal_variance_threshold) normal+=triangle.Area()*local_normal;}
                if(normal != zero_vector) normal.Normalize();
                (*ts.face_vertex_normals)(t)(i)=normal;}}}
    else{
        delete ts.face_vertex_normals;ts.face_vertex_normals=0;
        if(!ts.vertex_normals) ts.vertex_normals=new ARRAY<TV>(ts.mesh.number_nodes);
        else ts.vertex_normals->Resize(ts.mesh.number_nodes);
        for(int k=1;k<=ts.vertex_normals->m;k++){
            (*ts.vertex_normals)(k)=TV(); // initialize to zero
            for(int kk=1;kk<=(*ts.mesh.incident_elements)(k).m;kk++) 
                (*ts.vertex_normals)(k)+=(*ts.triangle_list)((*ts.mesh.incident_elements)(k)(kk)).Area()*(*ts.triangle_list)((*ts.mesh.incident_elements)(k)(kk)).normal;
            if((*ts.vertex_normals)(k) != TV()) (*ts.vertex_normals)(k).Normalize();}}

    if(!triangle_list_defined){delete ts.triangle_list;ts.triangle_list=0;}
    if(!incident_elements_defined){delete ts.mesh.incident_elements;ts.mesh.incident_elements=0;}
}
}
}
#endif
