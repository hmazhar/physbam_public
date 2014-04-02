#ifndef COMPILE_WITHOUT_DYADIC_SUPPORT
//#####################################################################
// Copyright 2003-2005, Ron Fedkiw, Geoffrey Irving, Sergey Koltakov, Frank Losasso, Andrew Selle, Tamar Shinar.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class LINEAR_INTERPOLATION_OCTREE_HELPER
//#####################################################################
#ifndef __LINEAR_INTERPOLATION_OCTREE_HELPER__
#define __LINEAR_INTERPOLATION_OCTREE_HELPER__

#include <PhysBAM_Tools/Grids_Dyadic/MAP_OCTREE_MESH.h>
#include <PhysBAM_Tools/Grids_Dyadic_Interpolation/LINEAR_INTERPOLATION_DYADIC_HELPER.h>
#include <PhysBAM_Tools/Grids_Dyadic_Interpolation/LINEAR_INTERPOLATION_QUADTREE_HELPER.h>
#include <PhysBAM_Tools/Interpolation/LINEAR_INTERPOLATION.h>

namespace PhysBAM{

template<class T_input>
class LINEAR_INTERPOLATION_OCTREE_HELPER:public LINEAR_INTERPOLATION_DYADIC_HELPER<OCTREE_GRID<T_input> >
{
    typedef T_input T;typedef VECTOR<T,3> TV;
private:
    LINEAR_INTERPOLATION_OCTREE_HELPER(); // disallow construction
public:
    template<class T2> static T2 Interpolate_Nodes(const OCTREE_GRID<T>& grid,const ARRAY<T2>& u,const TV& X)
    {return Interpolate_Nodes(grid,grid.Leaf_Cell(X),u,X);}

    template<class T2> static T2 Interpolate_Nodes(const OCTREE_GRID<T>& grid,const OCTREE_CELL<T>* cell,const ARRAY<T2>& u,const TV& X)
    {assert(cell);TV center(cell->Center()),DX_over_two((T).5*cell->DX());TV minimum_corner(center-DX_over_two),maximum_corner(center+DX_over_two);
    return LINEAR_INTERPOLATION<T,T2>::Trilinear(u(cell->Node(0)),u(cell->Node(1)),u(cell->Node(2)),u(cell->Node(3)),u(cell->Node(4)),u(cell->Node(5)),u(cell->Node(6)),u(cell->Node(7)),
        minimum_corner,maximum_corner,X);}

    static TV Interpolate_Faces(const OCTREE_GRID<T>& grid,const OCTREE_CELL<T>* cell,const ARRAY<T>& u_face,const ARRAY<TV>& u_node,const TV& X)
    {assert(cell);TV DX_over_two((T).5*cell->DX());
    T x=(X.x-cell->Center().x)/DX_over_two.x,y=(X.y-cell->Center().y)/DX_over_two.y,z=(X.z-cell->Center().z)/DX_over_two.z; // normalize to be between -1 and 1
    return TV(Interpolate_X_Face_Transformed(grid,cell,u_face,u_node,x,y,z),Interpolate_Y_Face_Transformed(grid,cell,u_face,u_node,x,y,z),
                        Interpolate_Z_Face_Transformed(grid,cell,u_face,u_node,x,y,z));}

    static T Interpolate_Face(const OCTREE_GRID<T>& grid,const int axis,const OCTREE_CELL<T>* cell,const ARRAY<T>& u_face,const ARRAY<TV>& u_node,const TV& X)
    {assert(cell);TV DX_over_two((T).5*cell->DX());
    T x=(X.x-cell->Center().x)/DX_over_two.x,y=(X.y-cell->Center().y)/DX_over_two.y,z=(X.z-cell->Center().z)/DX_over_two.z; // normalize to be between -1 and 1
    switch(axis){
        case 0:return Interpolate_X_Face_Transformed(grid,cell,u_face,u_node,x,y,z);
        case 1:return Interpolate_Y_Face_Transformed(grid,cell,u_face,u_node,x,y,z);
        default:assert(axis==2);return Interpolate_Z_Face_Transformed(grid,cell,u_face,u_node,x,y,z);}}

private:
    static T Interpolate_X_Face_Transformed(const OCTREE_GRID<T>& grid,const OCTREE_CELL<T>* cell,const ARRAY<T>& u_face,const ARRAY<TV>& u_node,const T& x,const T& y,const T& z)
    {T t=(T).5*(x+1);
    if(y>=abs(z)){return (1-t)*LINEAR_INTERPOLATION_QUADTREE_HELPER<T>::Interpolate_Triangle_Helper(u_face(cell->Face(0)),u_node(cell->Node(2)).x,u_node(cell->Node(6)).x,y,z)+
                           (t)*LINEAR_INTERPOLATION_QUADTREE_HELPER<T>::Interpolate_Triangle_Helper(u_face(cell->Face(1)),u_node(cell->Node(3)).x,u_node(cell->Node(7)).x,y,z);}
    else if(-y>=abs(z)){return (1-t)*LINEAR_INTERPOLATION_QUADTREE_HELPER<T>::Interpolate_Triangle_Helper(u_face(cell->Face(0)),u_node(cell->Node(0)).x,u_node(cell->Node(4)).x,-y,z)+
                                 (t)*LINEAR_INTERPOLATION_QUADTREE_HELPER<T>::Interpolate_Triangle_Helper(u_face(cell->Face(1)),u_node(cell->Node(1)).x,u_node(cell->Node(5)).x,-y,z);}
    else if(z>=abs(y)){return (1-t)*LINEAR_INTERPOLATION_QUADTREE_HELPER<T>::Interpolate_Triangle_Helper(u_face(cell->Face(0)),u_node(cell->Node(4)).x,u_node(cell->Node(6)).x,z,y)+
                                (t)*LINEAR_INTERPOLATION_QUADTREE_HELPER<T>::Interpolate_Triangle_Helper(u_face(cell->Face(1)),u_node(cell->Node(5)).x,u_node(cell->Node(7)).x,z,y);}
    else{assert(-z>=abs(y));return (1-t)*LINEAR_INTERPOLATION_QUADTREE_HELPER<T>::Interpolate_Triangle_Helper(u_face(cell->Face(0)),u_node(cell->Node(0)).x,u_node(cell->Node(2)).x,-z,y)+
                                     (t)*LINEAR_INTERPOLATION_QUADTREE_HELPER<T>::Interpolate_Triangle_Helper(u_face(cell->Face(1)),u_node(cell->Node(1)).x,u_node(cell->Node(3)).x,-z,y);}}

    static T Interpolate_Y_Face_Transformed(const OCTREE_GRID<T>& grid,const OCTREE_CELL<T>* cell,const ARRAY<T>& u_face,const ARRAY<TV>& u_node,const T& x,const T& y,const T& z)
    {T t=(T).5*(y+1);
    if(x>=abs(z)){return (1-t)*LINEAR_INTERPOLATION_QUADTREE_HELPER<T>::Interpolate_Triangle_Helper(u_face(cell->Face(2)),u_node(cell->Node(1)).y,u_node(cell->Node(5)).y,x,z)+
                           (t)*LINEAR_INTERPOLATION_QUADTREE_HELPER<T>::Interpolate_Triangle_Helper(u_face(cell->Face(3)),u_node(cell->Node(3)).y,u_node(cell->Node(7)).y,x,z);}
    else if(-x>=abs(z)){return (1-t)*LINEAR_INTERPOLATION_QUADTREE_HELPER<T>::Interpolate_Triangle_Helper(u_face(cell->Face(2)),u_node(cell->Node(0)).y,u_node(cell->Node(4)).y,-x,z)+
                                 (t)*LINEAR_INTERPOLATION_QUADTREE_HELPER<T>::Interpolate_Triangle_Helper(u_face(cell->Face(3)),u_node(cell->Node(2)).y,u_node(cell->Node(6)).y,-x,z);}
    else if(z>=abs(x)){return (1-t)*LINEAR_INTERPOLATION_QUADTREE_HELPER<T>::Interpolate_Triangle_Helper(u_face(cell->Face(2)),u_node(cell->Node(4)).y,u_node(cell->Node(5)).y,z,x)+
                                (t)*LINEAR_INTERPOLATION_QUADTREE_HELPER<T>::Interpolate_Triangle_Helper(u_face(cell->Face(3)),u_node(cell->Node(6)).y,u_node(cell->Node(7)).y,z,x);}
    else{assert(-z>=abs(x));return (1-t)*LINEAR_INTERPOLATION_QUADTREE_HELPER<T>::Interpolate_Triangle_Helper(u_face(cell->Face(2)),u_node(cell->Node(0)).y,u_node(cell->Node(1)).y,-z,x)+
                                     (t)*LINEAR_INTERPOLATION_QUADTREE_HELPER<T>::Interpolate_Triangle_Helper(u_face(cell->Face(3)),u_node(cell->Node(2)).y,u_node(cell->Node(3)).y,-z,x);}}

    static T Interpolate_Z_Face_Transformed(const OCTREE_GRID<T>& grid,const OCTREE_CELL<T>* cell,const ARRAY<T>& u_face,const ARRAY<TV>& u_node,const T& x,const T& y,const T& z)
    {T t=(T).5*(z+1);
    if(x>=abs(y)){return (1-t)*LINEAR_INTERPOLATION_QUADTREE_HELPER<T>::Interpolate_Triangle_Helper(u_face(cell->Face(4)),u_node(cell->Node(1)).z,u_node(cell->Node(3)).z,x,y)+
                           (t)*LINEAR_INTERPOLATION_QUADTREE_HELPER<T>::Interpolate_Triangle_Helper(u_face(cell->Face(5)),u_node(cell->Node(5)).z,u_node(cell->Node(7)).z,x,y);}
    else if(-x>=abs(y)){return (1-t)*LINEAR_INTERPOLATION_QUADTREE_HELPER<T>::Interpolate_Triangle_Helper(u_face(cell->Face(4)),u_node(cell->Node(0)).z,u_node(cell->Node(2)).z,-x,y)+
                                 (t)*LINEAR_INTERPOLATION_QUADTREE_HELPER<T>::Interpolate_Triangle_Helper(u_face(cell->Face(5)),u_node(cell->Node(4)).z,u_node(cell->Node(6)).z,-x,y);}
    else if(y>=abs(x)){return (1-t)*LINEAR_INTERPOLATION_QUADTREE_HELPER<T>::Interpolate_Triangle_Helper(u_face(cell->Face(4)),u_node(cell->Node(2)).z,u_node(cell->Node(3)).z,y,x)+
                                (t)*LINEAR_INTERPOLATION_QUADTREE_HELPER<T>::Interpolate_Triangle_Helper(u_face(cell->Face(5)),u_node(cell->Node(6)).z,u_node(cell->Node(7)).z,y,x);}
    else{assert(-y>=abs(x));return (1-t)*LINEAR_INTERPOLATION_QUADTREE_HELPER<T>::Interpolate_Triangle_Helper(u_face(cell->Face(4)),u_node(cell->Node(0)).z,u_node(cell->Node(1)).z,-y,x)+
                                     (t)*LINEAR_INTERPOLATION_QUADTREE_HELPER<T>::Interpolate_Triangle_Helper(u_face(cell->Face(5)),u_node(cell->Node(4)).z,u_node(cell->Node(5)).z,-y,x);}}
public:

    template<class T2> static T2 Interpolate_Cells(const OCTREE_GRID<T>& grid,const OCTREE_CELL<T>* cell,const ARRAY<T2>& u_cell,const ARRAY<T2>& u_node,const TV& X)
    {assert(cell);TV DX_over_two((T).5*cell->DX());
    T x=(X.x-cell->Center().x)/DX_over_two.x,y=(X.y-cell->Center().y)/DX_over_two.y,z=(X.z-cell->Center().z)/DX_over_two.z; // normalize the coordinates to be between -1 and 1
    int node_index1,node_index2,node_index3,node_index4;T x_prime,y_prime,z_prime;
    if(abs(x)>abs(y)){ // rule y out
        if(abs(x)>abs(z)){ // rule out z
            if(x>0){node_index1=1;node_index2=3;node_index3=5;node_index4=7;x_prime=x;y_prime=y;z_prime=z;} // +x pyramid
            else{node_index1=0;node_index2=2;node_index3=4;node_index4=6;x_prime=-x;y_prime=y;z_prime=z;}} // -x pyramid
        else{ // rule out x
            if(z>0){node_index1=4;node_index2=5;node_index3=6;node_index4=7;x_prime=z;y_prime=x;z_prime=y;} // +z pyramid
            else{node_index1=0;node_index2=1;node_index3=2;node_index4=3;x_prime=-z;y_prime=x;z_prime=y;}}} // -z pyramid
    else{ // rule out x
        if(abs(y)>abs(z)){ // rule out z
            if(y>0){node_index1=2;node_index2=3;node_index3=6;node_index4=7;x_prime=y;y_prime=x;z_prime=z;} // +y pyramid
            else{node_index1=0;node_index2=1;node_index3=4;node_index4=5;x_prime=-y;y_prime=x;z_prime=z;}} // -y pyramid
        else{ // rule out y
            if(z>0){node_index1=4;node_index2=5;node_index3=6;node_index4=7;x_prime=z;y_prime=x;z_prime=y;} // +z pyramid
            else{node_index1=0;node_index2=1;node_index3=2;node_index4=3;x_prime=-z;y_prime=x;z_prime=y;}}} // -z pyramid
    const T2& center=u_cell(cell->Cell());
    if(x_prime<(T)1e-8) return center; // avoid division by zero
    const T2& node1=u_node(cell->Node(node_index1)),&node2=u_node(cell->Node(node_index2)),&node3=u_node(cell->Node(node_index3)),&node4=u_node(cell->Node(node_index4));
    T one_over_x_prime=1/x_prime;y_prime*=one_over_x_prime;z_prime*=one_over_x_prime;
    return (1-x_prime)*center+x_prime*(T).25*((1-z_prime)*((1-y_prime)*node1+(1+y_prime)*node2)+(1+z_prime)*((1-y_prime)*node3+(1+y_prime)*node4));}

//#####################################################################
};
}
#endif
#endif
