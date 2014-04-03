//#####################################################################
// Copyright 2009.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#ifndef __TRIANGLE_SUBDIVISION_LOOP__
#define __TRIANGLE_SUBDIVISION_LOOP__

#include <PhysBAM_Tools/Vectors/VECTOR_3D.h>
#include <PhysBAM_Geometry/Topology/SEGMENT_MESH.h>
#include <PhysBAM_Geometry/Topology/TRIANGLE_MESH.h>
#include <PhysBAM_Geometry/Topology_Based_Geometry/TRIANGLE_SUBDIVISION.h>
namespace PhysBAM{

namespace TOPOLOGY_BASED_GEOMETRY_COMPUTATIONS
{
//#####################################################################
// Function Apply_Loop_Subdivision
//#####################################################################
template<class TV>
void Apply_Loop_Subdivision(TRIANGLE_SUBDIVISION& ts,ARRAY_VIEW<const TV> base_values,ARRAY_VIEW<TV> subdivided_values)
{
    typedef typename TV::SCALAR T;
    PHYSBAM_ASSERT(subdivided_values.Size()==ts.start_index_for_new_nodes-1+ts.triangle_mesh.segment_mesh->elements.m);
    if(!ts.triangle_mesh.topologically_sorted_neighbor_nodes){ts.delete_topologically_sorted_neighbor_nodes=true;ts.triangle_mesh.Initialize_Topologically_Sorted_Neighbor_Nodes();}
    if(!ts.triangle_mesh.boundary_mesh){ts.delete_boundary_mesh=true;ts.triangle_mesh.Initialize_Boundary_Mesh();}
    if(!ts.triangle_mesh.boundary_mesh->neighbor_nodes) ts.triangle_mesh.boundary_mesh->Initialize_Neighbor_Nodes();
    ARRAY<ARRAY<int> > &neighbors=*ts.triangle_mesh.topologically_sorted_neighbor_nodes,&boundary_neighbors=*ts.triangle_mesh.boundary_mesh->neighbor_nodes;
    // node values
    for(int i=1;i<=ts.triangle_mesh.number_nodes;i++)if(neighbors(i).m){
        if(ts.Node_Is_A_Corner(i) || boundary_neighbors(i).m>2 || boundary_neighbors(i).m==1)subdivided_values(i)=base_values(i);
        else if(boundary_neighbors(i).m==2) // if this is a regular boundary node
            subdivided_values(i)=(T).75*base_values(i)+(T).125*(base_values(boundary_neighbors(i)(1))+base_values(boundary_neighbors(i)(2)));
        else{ // interior node
            T alpha;switch(neighbors(i).m){
                case 3:alpha=(T).4375;break;case 4:alpha=(T).5;break;case 5:alpha=(T).54546609462891;break;case 6:alpha=(T).625;break;
                default:{T lambda=(T).375+(T).25*cos(T(2*pi)/neighbors(i).m);alpha=1-lambda*(4+lambda*(5*lambda-8))/(2*(1-lambda))+sqr(lambda);}}
            VECTOR<T,3> neighbor_sum=base_values(neighbors(i)(1));for(int j=2;j<=neighbors(i).m;j++)neighbor_sum+=base_values(neighbors(i)(j));
            subdivided_values(i)=alpha*base_values(i)+(1-alpha)/neighbors(i).m*neighbor_sum;}}
    // edge values
    for(int i=1;i<=ts.triangle_mesh.segment_mesh->elements.m;i++){
        int index=ts.start_index_for_new_nodes-1+i;
        int j,end1,end2;ts.triangle_mesh.segment_mesh->elements(i).Get(end1,end2);
        if(boundary_neighbors(end1).m && boundary_neighbors(end2).m && boundary_neighbors(end1).Find(end2,j)) // if boundary edge
            subdivided_values(index)=(T).5*(base_values(end1)+base_values(end2));
        else if(neighbors(end1).m==6 && neighbors(end2).m==6){ // if next to regular vertices (the most common situation)
            j=0;neighbors(end1).Find(end2,j);assert(j);
            int common1=neighbors(end1)(j==1?neighbors(end1).m:j-1),common2=neighbors(end1)(j==neighbors(end1).m?1:j+1);
            subdivided_values(index)=(T).375*(base_values(end1)+base_values(end2))+(T).125*(base_values(common1)+base_values(common2));}
        else{
            if(neighbors(end1).m != 6){
                T beta;switch(neighbors(end1).m){
                    case 3:beta=(T).625;break;case 4:beta=(T).640625;break;case 5:beta=(T).65906781074217;break;
                    default:{T lambda=(T).375+(T).25*cos((T)(2*pi)/neighbors(end1).m);beta=lambda*(4+lambda*(5*lambda-8))/(2*(1-lambda));}}
                subdivided_values(index)=(1-beta)*base_values(end1);
                int loop_start=0;neighbors(end1).Find(end2,loop_start);assert(loop_start);
                int j=loop_start;
                do{
                    T u=cos((T)(2*pi)*(j-loop_start)/neighbors(end1).m);
                    T w;switch(neighbors(end1).m){
                        case 3:w=(T)one_sixth*((T)1.25+u);break;case 4:w=(T).5*sqr((T).5+(T).375*u);break;
                        case 5:w=(T).16362712429686842801278667714785*sqr((T).55278640450004206071816526625374+u);break;
                        default:{T lambda=(T).375+(T).25*cos((T)(2*pi)/neighbors(end1).m);
                            w=2*cube(lambda)/(neighbors(end1).m*(1-lambda))*(1+u)*sqr(1/lambda-(T)1.5+u);}}
                    subdivided_values(index)+=w*base_values(neighbors(end1)(j));
                    j++;if(j>neighbors(end1).m)j=1;}
                while(j!=loop_start);}
            if(neighbors(end2).m != 6){
                T beta;switch(neighbors(end2).m){
                    case 3:beta=(T).625;break;case 4:beta=(T).640625;break;case 5:beta=(T).65906781074217;break;
                    default:{T lambda=(T).375+(T).25*cos((T)(2*pi)/neighbors(end2).m);beta=lambda*(4+lambda*(5*lambda-8))/(2*(1-lambda));}}
                if(neighbors(end1).m==6) subdivided_values(index)=(1-beta)*base_values(end2);
                else{subdivided_values(index)*=(T).5;subdivided_values(index)+=(T).5*(1-beta)*base_values(end2);} // need to average in this case
                int loop_start=0;neighbors(end2).Find(end1,loop_start);assert(loop_start);
                int j=loop_start;
                do{
                    T u=cos((T)(2*pi)*(j-loop_start)/neighbors(end2).m);
                    T w;switch(neighbors(end2).m){
                        case 3:w=(T)one_sixth*((T)1.25+u);break;case 4:w=(T).5*sqr((T).5+(T).375*u);break;
                        case 5:w=(T).16362712429686842801278667714785*sqr((T).55278640450004206071816526625374+u);break;
                        default:{T lambda=(T).375+(T).25*cos((T)(2*pi)/neighbors(end2).m);
                            w=2*cube(lambda)/(neighbors(end2).m*(1-lambda))*(1+u)*sqr(1/lambda-(T)1.5+u);}}
                    if(neighbors(end1).m!=6)w*=(T).5; // need to average in this case
                    subdivided_values(index)+=w*base_values(neighbors(end2)(j));
                    j++;if(j>neighbors(end2).m)j=1;}
                while(j!=loop_start);}}}
}
}
}
#endif
