//#####################################################################
// Copyright 2002-2006, Robert Bridson, Ronald Fedkiw, Eran Guendelman, Geoffrey Irving, Neil Molino, Avi Robinson-Mosher, Andrew Selle, Joseph Teran.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class TETRAHEDRON
//##################################################################### 
#include <PhysBAM_Tools/Arrays/ARRAY.h>
#include <PhysBAM_Tools/Arrays/ARRAY_VIEW.h>
#include <PhysBAM_Tools/Arrays/INDIRECT_ARRAY.h>
#include <PhysBAM_Tools/Grids_Uniform/GRID.h>
#include <PhysBAM_Tools/Interpolation/LINEAR_INTERPOLATION.h>
#include <PhysBAM_Tools/Vectors/VECTOR_3D.h>
#include <PhysBAM_Geometry/Basic_Geometry/TETRAHEDRON.h>
#include <PhysBAM_Geometry/Basic_Geometry/TRIANGLE_3D.h>
#include <PhysBAM_Geometry/Level_Sets/LEVELSET_UTILITIES.h>
using namespace PhysBAM;
//#####################################################################
// Constructor
//#####################################################################
template<class T> TETRAHEDRON<T>::
TETRAHEDRON()
    :X(TV(),TV(0,1,0),TV(1,0,0),TV(0,0,-1))
{
    Create_Triangles();
}
//#####################################################################
// Constructor
//#####################################################################
template<class T> TETRAHEDRON<T>::
TETRAHEDRON(const TV& x1_input,const TV& x2_input,const TV& x3_input,const TV& x4_input)
    :X(x1_input,x2_input,x3_input,x4_input)
{
    Create_Triangles();
}
//#####################################################################
// Function Create_Triangles
//#####################################################################
template<class T> void TETRAHEDRON<T>::
Create_Triangles()
{
    if(TV::Dot_Product(TV::Cross_Product(X[2]-X[1],X[3]-X[1]),X[4]-X[1]) <= 0){
        triangle1.Specify_Three_Points(X[1],X[2],X[3]);triangle2.Specify_Three_Points(X[1],X[4],X[2]);
        triangle3.Specify_Three_Points(X[1],X[3],X[4]);triangle4.Specify_Three_Points(X[2],X[4],X[3]);}
    else{
        triangle1.Specify_Three_Points(X[1],X[3],X[2]);triangle2.Specify_Three_Points(X[1],X[2],X[4]);
        triangle3.Specify_Three_Points(X[1],X[4],X[3]);triangle4.Specify_Three_Points(X[2],X[3],X[4]);}
}
//#####################################################################
// Function Normal
//#####################################################################
template<class T> VECTOR<T,3> TETRAHEDRON<T>::
Normal(const VECTOR<T,3>& location,const int aggregate) const
{
    assert(aggregate >= 1 && aggregate <= 4);
    if(aggregate == 1) return triangle1.normal;
    else if(aggregate == 2) return triangle2.normal;
    else if(aggregate == 3) return triangle3.normal;
    else return triangle4.normal;
}
//#####################################################################
// Function Inside
//#####################################################################
template<class T> bool TETRAHEDRON<T>::
Inside(const VECTOR<T,3>& location,const T thickness_over_two) const
{
    if(triangle1.Inside(location,thickness_over_two) && triangle2.Inside(location,thickness_over_two) && triangle3.Inside(location,thickness_over_two) && triangle4.Inside(location,thickness_over_two)) 
        return true;
    else return false;
}
//#####################################################################
// Function Outside
//#####################################################################
template<class T> bool TETRAHEDRON<T>::
Outside(const VECTOR<T,3>& location,const T thickness_over_two) const
{
    return triangle1.Outside(location,thickness_over_two) || triangle2.Outside(location,thickness_over_two) || triangle3.Outside(location,thickness_over_two)
        || triangle4.Outside(location,thickness_over_two);
}
//#####################################################################
// Function Boundary
//#####################################################################
template<class T> bool TETRAHEDRON<T>::
Boundary(const VECTOR<T,3>& location,const T thickness_over_two) const
{
    return !Inside(location,thickness_over_two) && !Outside(location,thickness_over_two);
}
//#####################################################################
// Function Thickened
//#####################################################################
template<class T> TETRAHEDRON<T> TETRAHEDRON<T>::
Thickened(const T thickness_over_two) const
{
    assert(Signed_Volume()>0);
    T m1=thickness_over_two/triangle4.Signed_Distance(X[1]),m2=thickness_over_two/triangle3.Signed_Distance(X[2]),m3=thickness_over_two/triangle2.Signed_Distance(X[3]),m4=thickness_over_two/triangle1.Signed_Distance(X[4]);
    return TETRAHEDRON<T>(Point_From_Barycentric_Coordinates(VECTOR<T,3>(m1,m2,m3)),Point_From_Barycentric_Coordinates(VECTOR<T,3>((T)1-m2-m3-(T)root_three*m4,m2,m3)),
                          Point_From_Barycentric_Coordinates(VECTOR<T,3>(m1,(T)1-m1-m3-(T)root_three*m4,m3)),Point_From_Barycentric_Coordinates(VECTOR<T,3>(m1,m2,(T)1-m1-m2-(T)root_three*m4)));
}
//#####################################################################
// Function Bounding_Box
//#####################################################################
template<class T> RANGE<VECTOR<T,3> > TETRAHEDRON<T>::
Bounding_Box() const
{      
    return RANGE<TV>::Bounding_Box(X);
}
//#####################################################################
// Function Surface
//#####################################################################
template<class T> VECTOR<T,3> TETRAHEDRON<T>::
Surface(const VECTOR<T,3>& location) const
{      
    if(Inside(location)){
        int triangle=1;T distance=VECTOR<T,3>::Dot_Product(triangle1.x1-location,triangle1.normal);
        T distance_temp=VECTOR<T,3>::Dot_Product(triangle2.x1-location,triangle2.normal);
        if(distance_temp < distance){triangle=2;distance=distance_temp;}
        distance_temp=VECTOR<T,3>::Dot_Product(triangle3.x1-location,triangle3.normal);
        if(distance_temp < distance){triangle=3;distance=distance_temp;}
        distance_temp=VECTOR<T,3>::Dot_Product(triangle4.x1-location,triangle4.normal);
        if(distance_temp < distance){triangle=4;distance=distance_temp;}
        switch(triangle){
            case 1:return location+distance*triangle1.normal;
            case 2:return location+distance*triangle2.normal;
            case 3:return location+distance*triangle3.normal;
            case 4:return location+distance*triangle4.normal;
            default:return VECTOR<T,3>(0,0,0);}} // should never be called!
    else{
        VECTOR<T,3> surface_point(location);
        if(triangle1.Lazy_Outside(surface_point)) surface_point-=VECTOR<T,3>::Dot_Product(location-triangle1.x1,triangle1.normal)*triangle1.normal;
        if(triangle2.Lazy_Outside(surface_point)) surface_point-=VECTOR<T,3>::Dot_Product(location-triangle2.x1,triangle2.normal)*triangle2.normal;
        if(triangle3.Lazy_Outside(surface_point)) surface_point-=VECTOR<T,3>::Dot_Product(location-triangle3.x1,triangle3.normal)*triangle3.normal;
        if(triangle4.Lazy_Outside(surface_point)) surface_point-=VECTOR<T,3>::Dot_Product(location-triangle4.x1,triangle4.normal)*triangle4.normal;
        return surface_point;}
}
//#####################################################################
// Function Closest_Point
//#####################################################################
template<class T> VECTOR<T,3> TETRAHEDRON<T>::
Closest_Point(const VECTOR<T,3>& location,VECTOR<T,3>& weights) const
{      
    if(Inside(location)){weights=First_Three_Barycentric_Coordinates(location);return location;}
    
    VECTOR<T,3> triangle_weights,triangle_closest_point=TRIANGLE_3D<T>(X[1],X[2],X[3]).Closest_Point(location,triangle_weights),closest_point=triangle_closest_point;
    weights=VECTOR<T,3>(triangle_weights.x,triangle_weights.y,triangle_weights.z);
    T triangle_distance_squared=(closest_point-location).Magnitude_Squared(),distance_squared=triangle_distance_squared;
    
    triangle_closest_point=TRIANGLE_3D<T>(X[1],X[2],X[4]).Closest_Point(location,triangle_weights);triangle_distance_squared=(triangle_closest_point-location).Magnitude_Squared();
    if(triangle_distance_squared<distance_squared){
        weights=VECTOR<T,3>(triangle_weights.x,triangle_weights.y,0);closest_point=triangle_closest_point;distance_squared=triangle_distance_squared;}
    
    triangle_closest_point=TRIANGLE_3D<T>(X[1],X[3],X[4]).Closest_Point(location,triangle_weights);triangle_distance_squared=(triangle_closest_point-location).Magnitude_Squared();
    if(triangle_distance_squared<distance_squared){
        weights=VECTOR<T,3>(triangle_weights.x,0,triangle_weights.y);closest_point=triangle_closest_point;distance_squared=triangle_distance_squared;}
    
    triangle_closest_point=TRIANGLE_3D<T>(X[2],X[3],X[4]).Closest_Point(location,triangle_weights);triangle_distance_squared=(triangle_closest_point-location).Magnitude_Squared();
    if(triangle_distance_squared<distance_squared){
        weights=VECTOR<T,3>(0,triangle_weights.x,triangle_weights.y);return triangle_closest_point;}
        
    return closest_point;
}
//#####################################################################
// Function Volume
//#####################################################################
template<class T> T TETRAHEDRON<T>::
Volume() const
{
    return (T)one_sixth*abs(VECTOR<T,3>::Dot_Product(VECTOR<T,3>::Cross_Product(X[2]-X[1],X[3]-X[1]),X[4]-X[1]));
}
//#####################################################################
// Function Signed_Volume
//#####################################################################
template<class T> T TETRAHEDRON<T>::
Signed_Volume() const
{  
    return (T)one_sixth*VECTOR<T,3>::Dot_Product(VECTOR<T,3>::Cross_Product(X[2]-X[1],X[3]-X[1]),X[4]-X[1]); 
}
//#####################################################################
// Function Minimum_Angle
//#####################################################################
template<class T> T TETRAHEDRON<T>::
Minimum_Angle() const
{  
    return min(triangle1.Minimum_Angle(),triangle2.Minimum_Angle(),triangle3.Minimum_Angle(),triangle4.Minimum_Angle());
}
//#####################################################################
// Function Maximum_Angle
//#####################################################################
template<class T> T TETRAHEDRON<T>::
Maximum_Angle() const
{  
    return max(triangle1.Maximum_Angle(),triangle2.Maximum_Angle(),triangle3.Maximum_Angle(),triangle4.Maximum_Angle());
}
//#####################################################################
// Function Minimum_Altitude
//#####################################################################
template<class T> T TETRAHEDRON<T>::
Minimum_Altitude() const
{  
    return min(abs(triangle1.Signed_Distance(X[4])),abs(triangle2.Signed_Distance(X[3])),abs(triangle3.Signed_Distance(X[2])),abs(triangle4.Signed_Distance(X[1])));
}
//#####################################################################
// Function Aspect_Ratio
//#####################################################################
template<class T> T TETRAHEDRON<T>::
Aspect_Ratio() const
{  
    return Maximum_Edge_Length()/Minimum_Altitude();
}
//#####################################################################
// Function Minimum_Dihedral_Angle
//#####################################################################
template<class T> T TETRAHEDRON<T>::
Minimum_Dihedral_Angle() const
{  
    return acos(-min(VECTOR<T,3>::Dot_Product(triangle1.normal,triangle2.normal),VECTOR<T,3>::Dot_Product(triangle1.normal,triangle3.normal),
                               VECTOR<T,3>::Dot_Product(triangle1.normal,triangle4.normal),VECTOR<T,3>::Dot_Product(triangle2.normal,triangle3.normal),
                               VECTOR<T,3>::Dot_Product(triangle2.normal,triangle4.normal),VECTOR<T,3>::Dot_Product(triangle3.normal,triangle4.normal))); 
}
//#####################################################################
// Function Maximum_Dihedral_Angle
//#####################################################################
template<class T> T TETRAHEDRON<T>::
Maximum_Dihedral_Angle() const
{  
    return acos(-max(VECTOR<T,3>::Dot_Product(triangle1.normal,triangle2.normal),VECTOR<T,3>::Dot_Product(triangle1.normal,triangle3.normal),
                               VECTOR<T,3>::Dot_Product(triangle1.normal,triangle4.normal),VECTOR<T,3>::Dot_Product(triangle2.normal,triangle3.normal),
                               VECTOR<T,3>::Dot_Product(triangle2.normal,triangle4.normal),VECTOR<T,3>::Dot_Product(triangle3.normal,triangle4.normal)));
}
//#####################################################################
// Function Signed_Reciprocal_Aspect_Ratio
//#####################################################################
template<class T> T TETRAHEDRON<T>::
Signed_Reciprocal_Aspect_Ratio(const VECTOR<T,3>& x1,const VECTOR<T,3>& x2,const VECTOR<T,3>& x3,const VECTOR<T,3>& x4)
{
    return min(VECTOR<T,3>::Dot_Product(x4-x1,PLANE<T>::Normal(x1,x2,x3)),VECTOR<T,3>::Dot_Product(x3-x1,PLANE<T>::Normal(x1,x4,x2)),
                     VECTOR<T,3>::Dot_Product(x2-x1,PLANE<T>::Normal(x1,x3,x4)),VECTOR<T,3>::Dot_Product(x1-x2,PLANE<T>::Normal(x2,x4,x3)))/
              max((x1-x2).Magnitude(),(x1-x3).Magnitude(),(x1-x4).Magnitude(),(x2-x3).Magnitude(),(x2-x4).Magnitude(),(x3-x4).Magnitude());
}
//#####################################################################
// Function Negative_Material
//#####################################################################
template<class T> T TETRAHEDRON<T>::
Negative_Material(const ARRAY<VECTOR<T,3> >& X,const ARRAY<T>& phis,const VECTOR<int,4>& indices)
{
    VECTOR<T,4> local_phi;
    int positive_count=0;for(int i=1;i<=4;i++) positive_count+=((local_phi[i]=phis(indices[i]))>0);
    switch(positive_count){
      case 0: return Signed_Volume(X(indices[1]),X(indices[2]),X(indices[3]),X(indices[4]));
      case 1:
        for(int i=1;i<=4;i++)if(local_phi[i]>0){
            VECTOR<VECTOR<T,3>,3> interface_locations;int index=i%4+1;
            for(int j=1;j<=3;j++,index=index%4+1)
                interface_locations[j]=LINEAR_INTERPOLATION<T,VECTOR<T,3> >::Linear(X(indices[i]),X(indices[index]),LEVELSET_UTILITIES<T>::Theta(local_phi[i],local_phi[index]));
            if(i%2 == 1) exchange(interface_locations[1],interface_locations[3]);
            return Signed_Volume(X(indices[1]),X(indices[2]),X(indices[3]),X(indices[4]))+TETRAHEDRON<T>::Signed_Volume(X(indices[i]),interface_locations[1],interface_locations[2],interface_locations[3]);}
      case 2:{
        positive_count=0;int negative_count=0;int negative_indices[2],positive_indices[2];
        for(int i=1;i<=4;i++){if(local_phi[i]<=0) negative_indices[negative_count++]=i;else positive_indices[positive_count++]=i;}
        if((negative_indices[1]-negative_indices[0])%2 == 1) exchange(positive_indices[0],positive_indices[1]);  // odd wrong, even right (odd=swap)
        VECTOR<T,3> interface_locations[2][2];
        for(int j=0;j<=1;j++)for(int k=0;k<=1;k++){
            int n=negative_indices[j],p=positive_indices[k];
            interface_locations[j][k]=LINEAR_INTERPOLATION<T,VECTOR<T,3> >::Linear(X(indices[n]),X(indices[p]),LEVELSET_UTILITIES<T>::Theta(local_phi[n],local_phi[p]));}
        // the tets are: (1-,1+2-,2+2-,1-2+) and (1-,1-1+,1+2-,1-2+), and (1-,1+2-,2-,2-2+)
        VECTOR<T,3> n0=X(indices[negative_indices[0]]);
        return TETRAHEDRON<T>::Signed_Volume(n0,interface_locations[1][0],X(indices[negative_indices[1]]),interface_locations[1][1])+
            TETRAHEDRON<T>::Signed_Volume(n0,interface_locations[1][0],interface_locations[1][1],interface_locations[0][1])+
            TETRAHEDRON<T>::Signed_Volume(n0,interface_locations[0][0],interface_locations[1][0],interface_locations[0][1]);}
      case 3:
        for(int i=1;i<=4;i++)if(local_phi[i]<=0){
            VECTOR<VECTOR<T,3>,3> interface_locations;int index=i%4+1;
            for(int j=1;j<=3;j++,index=index%4+1)
                interface_locations[j]=LINEAR_INTERPOLATION<T,VECTOR<T,3> >::Linear(X(indices[i]),X(indices[index]),LEVELSET_UTILITIES<T>::Theta(local_phi[i],local_phi[index]));
            if(i%2 == 1) exchange(interface_locations[1],interface_locations[3]);
            return -TETRAHEDRON<T>::Signed_Volume(X(indices[i]),interface_locations[1],interface_locations[2],interface_locations[3]);}
      default:return 0;}

    PHYSBAM_FATAL_ERROR();
}
//#####################################################################
// Function Cut_With_Hyperplane
//#####################################################################
template<class T> void TETRAHEDRON<T>::
Cut_With_Hyperplane(ARRAY<VECTOR<T,3> >& X,const PLANE<T>& cutting_plane,const VECTOR<int,4>& indices,ARRAY<VECTOR<int,4> >& left_simplices,ARRAY<VECTOR<int,4> >& right_simplices)
{
    VECTOR<T,4> phis;
    VECTOR<VECTOR<T,3>,4> X_nodes;
    for(int i=1;i<=4;i++){X_nodes[i]=X(indices[i]);phis[i]=cutting_plane.Signed_Distance(X_nodes[i]);}
    Cut_Simplex(X,indices,X_nodes,phis,left_simplices,right_simplices);
}
//#####################################################################
// Function Clip_To_Box
//#####################################################################
template<class T> void TETRAHEDRON<T>::
Clip_To_Box(const RANGE<TV>& box,ARRAY<TETRAHEDRON<T> >& clipped_simplices) const
{
    // cut with all sides of box
    clipped_simplices.Remove_All();
    clipped_simplices.Append(*this);
    for(int axis=1;axis<=TV::dimension;axis++){
        for(int i=clipped_simplices.m;i>=1;i--){
            Cut_With_Hyperplane_And_Discard_Outside_Simplices(clipped_simplices(i),PLANE<T>(-TV::Axis_Vector(axis),box.min_corner),clipped_simplices);
            clipped_simplices.Remove_Index_Lazy(i);}
        for(int i=clipped_simplices.m;i>=1;i--){
            // TODO: make this more efficient by not removing a triangle that is fully inside
            Cut_With_Hyperplane_And_Discard_Outside_Simplices(clipped_simplices(i),PLANE<T>(TV::Axis_Vector(axis),box.max_corner),clipped_simplices);
            clipped_simplices.Remove_Index_Lazy(i);}}
}
//#####################################################################
// Function Cut_With_Hyperplane_And_Discard_Outside_Simplices
//#####################################################################
template<class T> void TETRAHEDRON<T>::
Cut_With_Hyperplane_And_Discard_Outside_Simplices(const TETRAHEDRON<T>& tetrahedron,const PLANE<T>& cutting_plane,ARRAY<TETRAHEDRON<T> >& negative_tetrahedra) const
{
    // left simplices are in the negative halfspace, right simplices in the positive halfspace
    VECTOR<T,4> phi_nodes;
    VECTOR<VECTOR<T,3>,4> X_nodes;X_nodes(1)=tetrahedron.X[1];X_nodes(2)=tetrahedron.X[2];X_nodes(3)=tetrahedron.X[3];X_nodes(4)=tetrahedron.X[4];
    for(int i=1;i<=4;i++){phi_nodes(i)=cutting_plane.Signed_Distance(X_nodes(i));}

    int positive_count=0;
    for(int i=1;i<=4;i++) if(phi_nodes[i]>0) positive_count++;
    switch(positive_count){
        case 0: // we are in the negative halfspace
            negative_tetrahedra.Append(tetrahedron);
            break;
        case 1: // tet in positive halfspace, three tets in negative
            for(int i=1;i<=4;i++)if(phi_nodes[i]>0){
                VECTOR<VECTOR<T,3>,3> interface_locations;int index=i%4+1;
                VECTOR<int,3> other_indices;
                for(int j=1;j<=3;j++,index=index%4+1){
                    other_indices[j]=index;
                    interface_locations[j]=LINEAR_INTERPOLATION<T,VECTOR<T,3> >::Linear(X_nodes[i],X_nodes[index],LEVELSET_UTILITIES<T>::Theta(phi_nodes[i],phi_nodes[index]));}
                if(i%2 == 0){exchange(interface_locations[1],interface_locations[3]);exchange(other_indices[1],other_indices[3]);}
                // (i1,o1,o2,o3), (i2,i1,o2,o3), (i3,i1,i2,o3)
                negative_tetrahedra.Append(TETRAHEDRON<T>(interface_locations[1],X_nodes(other_indices[1]),X_nodes(other_indices[2]),X_nodes(other_indices[3])));
                negative_tetrahedra.Append(TETRAHEDRON<T>(interface_locations[2],interface_locations[1],X_nodes(other_indices[2]),X_nodes(other_indices[3])));
                negative_tetrahedra.Append(TETRAHEDRON<T>(interface_locations[3],interface_locations[1],interface_locations[2],X_nodes(other_indices[3])));
                return;}
        case 2:{ // three tets in each halfspace
            positive_count=0;int negative_count=0;int negative_indices[2],positive_indices[2];
            for(int i=1;i<=4;i++){if(phi_nodes[i]<=0) negative_indices[negative_count++]=i;else positive_indices[positive_count++]=i;}
            if((negative_indices[1]-negative_indices[0])%2 == 1) exchange(positive_indices[0],positive_indices[1]);  // odd wrong, even right (odd=swap)
            VECTOR<T,3> interface_locations[2][2];
            for(int j=0;j<=1;j++)for(int k=0;k<=1;k++){
                int n=negative_indices[j],p=positive_indices[k];
                interface_locations[j][k]=LINEAR_INTERPOLATION<T,VECTOR<T,3> >::Linear(X_nodes[n],X_nodes[p],LEVELSET_UTILITIES<T>::Theta(phi_nodes[n],phi_nodes[p]));}
            // the tets are: (1-,1+2-,2+2-,1-2+) and (1-,1-1+,1+2-,1-2+), and (1-,1+2-,2-,2-2+)
            VECTOR<T,3> n0=X_nodes(negative_indices[0]);
            negative_tetrahedra.Append(TETRAHEDRON<T>(n0,interface_locations[1][0],X_nodes(negative_indices[1]),interface_locations[1][1]));
            negative_tetrahedra.Append(TETRAHEDRON<T>(n0,interface_locations[1][0],interface_locations[1][1],interface_locations[0][1]));
            negative_tetrahedra.Append(TETRAHEDRON<T>(n0,interface_locations[0][0],interface_locations[1][0],interface_locations[0][1]));
            break;}
        case 3: // tet in negative halfspace, three tets in positive
            for(int i=1;i<=4;i++)if(phi_nodes[i]<=0){
                VECTOR<VECTOR<T,3>,3> interface_locations;int index=i%4+1;
                VECTOR<int,3> other_indices;
                for(int j=1;j<=3;j++,index=index%4+1){
                    other_indices[j]=index;
                    interface_locations[j]=LINEAR_INTERPOLATION<T,VECTOR<T,3> >::Linear(X_nodes[i],X_nodes[index],LEVELSET_UTILITIES<T>::Theta(phi_nodes[i],phi_nodes[index]));}
                if(i%2 == 0){exchange(interface_locations[1],interface_locations[3]);exchange(other_indices[1],other_indices[3]);}
                negative_tetrahedra.Append(TETRAHEDRON<T>(X_nodes(i),interface_locations[1],interface_locations[2],interface_locations[3]));
                return;}
        case 4:
            break;}
}
//#####################################################################
// Function Cut_Simplex
//#####################################################################
template<class T> static inline void Add_Points_As_Tetrahedron(const ARRAY<VECTOR<T,3> >& X,ARRAY<VECTOR<int,4> >& tets,const int x1,const int x2,const int x3,const int x4)
{
    /*if(TETRAHEDRON<T>::Signed_Volume(X(x1),X(x2),X(x3),X(x4)))*/ tets.Append(VECTOR<int,4>(x1,x2,x3,x4));
}
template<class T> void TETRAHEDRON<T>::
Cut_Simplex(ARRAY<VECTOR<T,3> >& X,const VECTOR<int,4>& indices,const VECTOR<VECTOR<T,3>,4>& X_nodes,const VECTOR<T,4>& phi_nodes,ARRAY<VECTOR<int,4> >& left_simplices,
    ARRAY<VECTOR<int,4> >& right_simplices)
{
    // left simplices are in the negative halfspace, right simplices in the positive halfspace
    int positive_count=0;
    for(int i=1;i<=4;i++) if(phi_nodes[i]>0) positive_count++;
    switch(positive_count){
      case 0: // we are in the negative halfspace
        Add_Points_As_Tetrahedron(X,left_simplices,indices[1],indices[2],indices[3],indices[4]);
        break;
      case 1: // tet in positive halfspace, three tets in negative
        for(int i=1;i<=4;i++)if(phi_nodes[i]>0){
            VECTOR<int,3> interface_locations;int index=i%4+1;
            VECTOR<int,3> other_indices;
            for(int j=1;j<=3;j++,index=index%4+1){
                other_indices[j]=indices[index];
                interface_locations[j]=X.Append(LINEAR_INTERPOLATION<T,VECTOR<T,3> >::Linear(X_nodes[i],X_nodes[index],LEVELSET_UTILITIES<T>::Theta(phi_nodes[i],phi_nodes[index])));}
            if(i%2 == 0){exchange(interface_locations[1],interface_locations[3]);exchange(other_indices[1],other_indices[3]);}
            Add_Points_As_Tetrahedron(X,right_simplices,indices[i],interface_locations[1],interface_locations[2],interface_locations[3]);
            // other three in left simplices
            // (i1,o1,o2,o3), (i2,i1,o2,o3), (i3,i1,i2,o3)
            Add_Points_As_Tetrahedron(X,left_simplices,interface_locations[1],other_indices[1],other_indices[2],other_indices[3]);
            Add_Points_As_Tetrahedron(X,left_simplices,interface_locations[2],interface_locations[1],other_indices[2],other_indices[3]);
            Add_Points_As_Tetrahedron(X,left_simplices,interface_locations[3],interface_locations[1],interface_locations[2],other_indices[3]);
            return;}
      case 2:{ // three tets in each halfspace
        positive_count=0;int negative_count=0;int negative_indices[2],positive_indices[2];
        for(int i=1;i<=4;i++){if(phi_nodes[i]<=0) negative_indices[negative_count++]=i;else positive_indices[positive_count++]=i;}
        if((negative_indices[1]-negative_indices[0])%2 == 1) exchange(positive_indices[0],positive_indices[1]);  // odd wrong, even right (odd=swap)
        int interface_locations[2][2];
        for(int j=0;j<=1;j++)for(int k=0;k<=1;k++){
            int n=negative_indices[j],p=positive_indices[k];
            interface_locations[j][k]=X.Append(LINEAR_INTERPOLATION<T,VECTOR<T,3> >::Linear(X_nodes[n],X_nodes[p],LEVELSET_UTILITIES<T>::Theta(phi_nodes[n],phi_nodes[p])));}
        // the tets are: (1-,1+2-,2+2-,1-2+) and (1-,1-1+,1+2-,1-2+), and (1-,1+2-,2-,2-2+)
        // the tets are: (1+,1-2+,2-2+,1+2-) and (1+,1+1-,1-2+,1+2-), and (1+,1-2+,2+,2+2-)
        int n0=indices[negative_indices[0]];int p0=indices[positive_indices[0]];
        Add_Points_As_Tetrahedron(X,left_simplices,n0,interface_locations[1][0],indices[negative_indices[1]],interface_locations[1][1]);
        Add_Points_As_Tetrahedron(X,left_simplices,n0,interface_locations[1][0],interface_locations[1][1],interface_locations[0][1]);
        Add_Points_As_Tetrahedron(X,left_simplices,n0,interface_locations[0][0],interface_locations[1][0],interface_locations[0][1]);
        Add_Points_As_Tetrahedron(X,right_simplices,p0,interface_locations[0][1],indices[positive_indices[1]],interface_locations[1][1]);
        Add_Points_As_Tetrahedron(X,right_simplices,p0,interface_locations[0][1],interface_locations[1][1],interface_locations[1][0]);
        Add_Points_As_Tetrahedron(X,right_simplices,p0,interface_locations[0][0],interface_locations[0][1],interface_locations[1][0]);}
        break;
      case 3: // tet in negative halfspace, three tets in positive
        for(int i=1;i<=4;i++)if(phi_nodes[i]<=0){
            VECTOR<int,3> interface_locations;int index=i%4+1;
            VECTOR<int,3> other_indices;
            for(int j=1;j<=3;j++,index=index%4+1){
                other_indices[j]=indices[index];
                interface_locations[j]=X.Append(LINEAR_INTERPOLATION<T,VECTOR<T,3> >::Linear(X_nodes[i],X_nodes[index],LEVELSET_UTILITIES<T>::Theta(phi_nodes[i],phi_nodes[index])));}
            if(i%2 == 0){exchange(interface_locations[1],interface_locations[3]);exchange(other_indices[1],other_indices[3]);}
            Add_Points_As_Tetrahedron(X,left_simplices,indices[i],interface_locations[1],interface_locations[2],interface_locations[3]);
            Add_Points_As_Tetrahedron(X,right_simplices,interface_locations[1],other_indices[1],other_indices[2],other_indices[3]);
            Add_Points_As_Tetrahedron(X,right_simplices,interface_locations[2],interface_locations[1],other_indices[2],other_indices[3]);
            Add_Points_As_Tetrahedron(X,right_simplices,interface_locations[3],interface_locations[1],interface_locations[2],other_indices[3]);
            return;}
      case 4:
        Add_Points_As_Tetrahedron(X,right_simplices,indices[1],indices[2],indices[3],indices[4]);
        break;}
}
//#####################################################################
template class TETRAHEDRON<float>;
#ifndef COMPILE_WITHOUT_DOUBLE_SUPPORT
template class TETRAHEDRON<double>;
#else
template double TETRAHEDRON<double>::Minimum_Altitude() const; // used by MATRIX<double,3>
#endif
