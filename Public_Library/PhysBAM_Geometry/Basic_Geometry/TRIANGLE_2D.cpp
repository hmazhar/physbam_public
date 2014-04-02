//#####################################################################
// Copyright 2006, Avi Robinson-Mosher.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class TRIANGLE_2D
//##################################################################### 
#include <PhysBAM_Tools/Arrays/ARRAY.h>
#include <PhysBAM_Tools/Interpolation/LINEAR_INTERPOLATION.h>
#include <PhysBAM_Geometry/Basic_Geometry/LINE_2D.h>
#include <PhysBAM_Geometry/Basic_Geometry/SEGMENT_2D.h>
#include <PhysBAM_Geometry/Basic_Geometry/TRIANGLE_2D.h>
#include <PhysBAM_Geometry/Level_Sets/LEVELSET_UTILITIES.h>
using namespace PhysBAM;
//#####################################################################
// Function Negative_Material
//#####################################################################
template<class T> T TRIANGLE_2D<T>::
Negative_Material(const ARRAY<TV>& X,const ARRAY<T>& phis,const VECTOR<int,3>& indices)
{
    int positive_count=0;
    VECTOR<T,3> local_phi;
    T area=Signed_Area(X(indices[1]),X(indices[2]),X(indices[3])); // doesn't matter what's inside; also avoids duplicate node case (which messes up theta)
    // special case: if two phis are zero, that's an interface 
    if(!area) return 0;
    for(int i=1;i<=3;i++) positive_count+=((local_phi[i]=phis(indices[i]))>0);
    switch(positive_count){
        case 0: return area;
        case 1:
            // draw positive triangle. has correct positive/negative area based on whether triangle is backwards or not
            for(int i=1;i<=3;i++)if(local_phi[i]>0){
                VECTOR<TV,2> interface_locations;int index=i%3+1;
                for(int j=1;j<=2;j++,index=index%3+1)
                    interface_locations[j]=LINEAR_INTERPOLATION<T,TV>::Linear(X(indices[i]),X(indices[index]),LEVELSET_UTILITIES<T>::Theta(local_phi[i],local_phi[index]));
                return area-TRIANGLE_2D<T>::Signed_Area(X(indices[i]),interface_locations[1],interface_locations[2]);}
        case 2:
            // draw negative triangle
            for(int i=1;i<=3;i++)if(local_phi[i]<=0){
                VECTOR<TV,2> interface_locations;int index=i%3+1;
                for(int j=1;j<=2;j++,index=index%3+1)
                    interface_locations[j]=LINEAR_INTERPOLATION<T,TV>::Linear(X(indices[i]),X(indices[index]),LEVELSET_UTILITIES<T>::Theta(local_phi[i],local_phi[index]));
                return TRIANGLE_2D<T>::Signed_Area(X(indices[i]),interface_locations[1],interface_locations[2]);}
        case 3: return (T)0;}

    PHYSBAM_FATAL_ERROR();
}
//#####################################################################
// Function Cut_With_Hyperplane
//#####################################################################
template<class T> void TRIANGLE_2D<T>::
Cut_With_Hyperplane(ARRAY<TV>& X,const LINE_2D<T>& cutting_plane,const VECTOR<int,3>& indices,ARRAY<VECTOR<int,3> >& left_tris,ARRAY<VECTOR<int,3> >& right_tris)
{
    VECTOR<T,3> phis;
    VECTOR<TV,3> X_nodes;
    for(int i=1;i<=3;i++){X_nodes[i]=X(indices[i]);phis[i]=cutting_plane.Signed_Distance(X_nodes[i]);}
    Cut_Simplex(X,indices,X_nodes,phis,left_tris,right_tris);
}
//#####################################################################
// Function Cut_Simplex
//#####################################################################
template<class T> static inline void Add_Points_As_Triangle(const ARRAY<VECTOR<T,2> >& X,ARRAY<VECTOR<int,3> >& tris,const int x1,const int x2,const int x3)
{
    tris.Append(VECTOR<int,3>(x1,x2,x3));
}
template<class T> void TRIANGLE_2D<T>::
Cut_Simplex(ARRAY<TV>& X,const VECTOR<int,3>& indices,const VECTOR<TV,3>& X_nodes,const VECTOR<T,3>& phi_nodes,ARRAY<VECTOR<int,3> >& left_tris,ARRAY<VECTOR<int,3> >& right_tris)
{
    int positive_count=0;
    for(int i=1;i<=3;i++) if(phi_nodes[i]>0) positive_count++;
    switch(positive_count){
      case 0: // place in left list
        Add_Points_As_Triangle(X,left_tris,indices[1],indices[2],indices[3]);break;
      case 1:
        // draw positive triangle. has correct positive/negative area based on whether triangle is backwards or not
        for(int i=1;i<=3;i++)if(phi_nodes[i]>0){
            VECTOR<int,2> interface_locations;int index=i%3+1;
            VECTOR<int,2> other_locations;
            for(int j=1;j<=2;j++,index=index%3+1){
                other_locations[j]=indices[index];
                interface_locations[j]=X.Append(LINEAR_INTERPOLATION<T,TV>::Linear(X_nodes[i],X_nodes[index],LEVELSET_UTILITIES<T>::Theta(phi_nodes[i],phi_nodes[index])));}
            // add triangle to right tris
            Add_Points_As_Triangle(X,right_tris,indices[i],interface_locations[1],interface_locations[2]);
            // add two triangles to left tris
            Add_Points_As_Triangle(X,left_tris,interface_locations[1],other_locations[1],other_locations[2]);
            Add_Points_As_Triangle(X,left_tris,interface_locations[1],other_locations[2],interface_locations[2]);
            return;}
      case 2:
        // draw negative triangle
        for(int i=1;i<=3;i++)if(phi_nodes[i]<=0){
            VECTOR<int,2> interface_locations;int index=i%3+1;
            VECTOR<int,2> other_locations;
            for(int j=1;j<=2;j++,index=index%3+1){
                other_locations[j]=indices[index];
                interface_locations[j]=X.Append(LINEAR_INTERPOLATION<T,TV>::Linear(X_nodes[i],X_nodes[index],LEVELSET_UTILITIES<T>::Theta(phi_nodes[i],phi_nodes[index])));}
            // add triangle to left tris
            Add_Points_As_Triangle(X,left_tris,indices[i],interface_locations[1],interface_locations[2]);
            // add two triangles to right tris
            Add_Points_As_Triangle(X,right_tris,interface_locations[1],other_locations[1],other_locations[2]);
            Add_Points_As_Triangle(X,right_tris,interface_locations[1],other_locations[2],interface_locations[2]);
            return;}
      case 3: // place in right list
        Add_Points_As_Triangle(X,right_tris,indices[1],indices[2],indices[3]);break;}
}
//#####################################################################
// Function Clip_To_Box
//#####################################################################
template<class T> void TRIANGLE_2D<T>::
Clip_To_Box(const RANGE<TV>& box,ARRAY<TRIANGLE_2D<T> >& clipped_simplices) const
{
    // cut with all sides of box
    clipped_simplices.Remove_All();
    clipped_simplices.Append(*this);
    for(int axis=1;axis<=TV::dimension;axis++){
        for(int i=clipped_simplices.m;i>=1;i--){
            Cut_With_Hyperplane_And_Discard_Outside_Simplices(clipped_simplices(i),LINE_2D<T>(-TV::Axis_Vector(axis),box.min_corner),clipped_simplices);
            clipped_simplices.Remove_Index_Lazy(i);}
        for(int i=clipped_simplices.m;i>=1;i--){
            // TODO: make this more efficient by not removing a triangle that is fully inside
            Cut_With_Hyperplane_And_Discard_Outside_Simplices(clipped_simplices(i),LINE_2D<T>(TV::Axis_Vector(axis),box.max_corner),clipped_simplices);
            clipped_simplices.Remove_Index_Lazy(i);}}
}
//#####################################################################
// Function Cut_With_Hyperplane_And_Discard_Outside_Simplices
//#####################################################################
template<class T> void TRIANGLE_2D<T>::
Cut_With_Hyperplane_And_Discard_Outside_Simplices(const TRIANGLE_2D<T>& triangle,const LINE_2D<T>& cutting_plane,ARRAY<TRIANGLE_2D<T> >& negative_triangles) const
{
    VECTOR<T,3> phi_nodes;
    VECTOR<VECTOR<T,2>,3> X_nodes;X_nodes(1)=triangle.X[1];X_nodes(2)=triangle.X[2];X_nodes(3)=triangle.X[3];
    for(int i=1;i<=3;i++){phi_nodes(i)=cutting_plane.Signed_Distance(X_nodes(i));}

    // left simplices are in the negative halfspace, right simplices in the positive halfspace
    int positive_count=0,single_node_sign;
    for(int i=1;i<=3;i++) if(phi_nodes(i)>0) positive_count++;
    switch(positive_count){
        case 0: // in negative halfspace
            negative_triangles.Append(triangle);break;
        case 1:
        case 2:
            single_node_sign=positive_count==1?1:-1;
            // draw positive triangle. has correct positive/negative area based on whether triangle is backwards or not
            for(int i=1;i<=3;i++)if(LEVELSET_UTILITIES<T>::Sign(phi_nodes(i))==single_node_sign){
                VECTOR<VECTOR<T,2>,2> interface_locations;int index=i%3+1;
                VECTOR<int,2> other_locations;
                for(int j=1;j<=2;j++,index=index%3+1){
                    other_locations(j)=index;
                    interface_locations(j)=LINEAR_INTERPOLATION<T,VECTOR<T,2> >::Linear(X_nodes(i),X_nodes(index),LEVELSET_UTILITIES<T>::Theta(phi_nodes(i),phi_nodes(index)));}
                if(positive_count==1){ // add two triangles to negative triangles
                    negative_triangles.Append(TRIANGLE_2D<T>(interface_locations(1),X_nodes(other_locations(1)),X_nodes(other_locations(2))));
                    negative_triangles.Append(TRIANGLE_2D<T>(interface_locations(1),X_nodes(other_locations(2)),interface_locations(2)));}
                else // add triangle to negative_triangles
                    negative_triangles.Append(TRIANGLE_2D<T>(X_nodes(i),interface_locations(1),interface_locations(2)));
                return;}
        case 3: // in positive halfspace
            break;}
}
//#####################################################################
#define INSTANTIATION_HELPER(T) \
    template T TRIANGLE_2D<T>::Negative_Material(const ARRAY<VECTOR<T,2> >& X,const ARRAY<T>& phis,const VECTOR<int,3>& indices);\
    template void TRIANGLE_2D<T>::Cut_With_Hyperplane(ARRAY<VECTOR<T,2> >& X,const LINE_2D<T>& cutting_plane,const VECTOR<int,3>& indices,ARRAY<VECTOR<int,3> >& left_tris, \
        ARRAY<VECTOR<int,3> >& right_tris); \
    template void TRIANGLE_2D<T>::Cut_Simplex(ARRAY<VECTOR<T,2> >& X,const VECTOR<int,3>& indices,const VECTOR<VECTOR<T,2>,3>& X_nodes,const VECTOR<T,3>& phi_nodes, \
        ARRAY<VECTOR<int,3> >& left_tris,ARRAY<VECTOR<int,3> >& right_tris); \
    template void TRIANGLE_2D<T>::Clip_To_Box(const RANGE<VECTOR<T,2> >& box,ARRAY<TRIANGLE_2D<T> >& clipped_simplices) const; \
    template void TRIANGLE_2D<T>::Cut_With_Hyperplane_And_Discard_Outside_Simplices(const TRIANGLE_2D<T>& triangle,const LINE_2D<T>& cutting_plane,ARRAY<TRIANGLE_2D<T> >& negative_triangles) const;

INSTANTIATION_HELPER(float)
#ifndef COMPILE_WITHOUT_DOUBLE_SUPPORT
INSTANTIATION_HELPER(double)
#endif
