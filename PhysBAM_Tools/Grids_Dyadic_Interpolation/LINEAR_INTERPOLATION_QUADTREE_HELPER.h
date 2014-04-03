//#####################################################################
// Copyright 2004-2005, Ron Fedkiw, Geoffrey Irving, Frank Losasso, Andrew Selle, Tamar Shinar.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class LINEAR_INTERPOLATION_QUADTREE_HELPER
//#####################################################################
#ifndef COMPILE_WITHOUT_DYADIC_SUPPORT
#ifndef __LINEAR_INTERPOLATION_QUADTREE_HELPER__
#define __LINEAR_INTERPOLATION_QUADTREE_HELPER__

#include <PhysBAM_Tools/Grids_Dyadic/MAP_QUADTREE_MESH.h>
#include <PhysBAM_Tools/Grids_Dyadic_Interpolation/LINEAR_INTERPOLATION_DYADIC_HELPER.h>
#include <PhysBAM_Tools/Interpolation/LINEAR_INTERPOLATION.h>
namespace PhysBAM{

template<class T_input>
class LINEAR_INTERPOLATION_QUADTREE_HELPER:public LINEAR_INTERPOLATION_DYADIC_HELPER<QUADTREE_GRID<T_input> >
{
    typedef T_input T;typedef VECTOR<T,2> TV;
private:
    LINEAR_INTERPOLATION_QUADTREE_HELPER(); // disallow construction
public:
    template<class T2> static T2 Interpolate_Nodes(const QUADTREE_GRID<T>& grid,const ARRAY<T2>& u,const TV& X)
    {return Interpolate_Nodes(grid,grid.Leaf_Cell(X),u,X);}

    template<class T2> static T2 Interpolate_Nodes(const QUADTREE_GRID<T>& grid,const QUADTREE_CELL<T>* cell,const ARRAY<T2>& u,const TV& X)
    {assert(cell); // this means we are outside the domain (of everything... including the ghost cells)
    TV center(cell->Center()),DX_over_two((T).5*cell->DX());
    TV minimum_corner(center-DX_over_two),maximum_corner(center+DX_over_two);
    return LINEAR_INTERPOLATION<T,T2>::Bilinear(u(cell->Node(0)),u(cell->Node(1)),u(cell->Node(2)),u(cell->Node(3)),minimum_corner,maximum_corner,X);}

    static TV Interpolate_Faces(const QUADTREE_GRID<T>& grid,const QUADTREE_CELL<T>* cell,const ARRAY<T>& u_face,const ARRAY<TV>& u_node,const TV& X)
    {assert(cell);TV DX(cell->DX()),minimum_corner(cell->Center()-(T).5*DX);
    T x=(X.x-minimum_corner.x)/DX.x,y=(X.y-minimum_corner.y)/DX.y; // normalize the coordinates to be between 0 and 1
    return TV(Interpolate_X_Face_Transformed(grid,cell,u_face,u_node,x,y),Interpolate_Y_Face_Transformed(grid,cell,u_face,u_node,x,y));}

    static T Interpolate_Face(const QUADTREE_GRID<T>& grid,const int axis,const QUADTREE_CELL<T>* cell,const ARRAY<T>& u_face,const ARRAY<TV>& u_node,const TV& X)
    {assert(cell);TV DX(cell->DX()),minimum_corner(cell->Center()-(T).5*DX);
    T x=(X.x-minimum_corner.x)/DX.x,y=(X.y-minimum_corner.y)/DX.y; // normalize the coordinates to be between 0 and 1
    switch(axis){
      case 0:return Interpolate_X_Face_Transformed(grid,cell,u_face,u_node,x,y);
      default:assert(axis==1);return Interpolate_Y_Face_Transformed(grid,cell,u_face,u_node,x,y);}}

private:
    static T Interpolate_X_Face_Transformed(const QUADTREE_GRID<T>& grid,const QUADTREE_CELL<T>* cell,const ARRAY<T>& u_face,const ARRAY<TV>& u_node,const T& x,const T& y)
    {if(y>=(T).5) return LINEAR_INTERPOLATION<T,T>::Bilinear(u_face(cell->Face(0)),u_face(cell->Face(1)),u_node(cell->Node(2)).x,u_node(cell->Node(3)).x,TV(x,2*(y-(T).5)));
    else return LINEAR_INTERPOLATION<T,T>::Bilinear(u_node(cell->Node(0)).x,u_node(cell->Node(1)).x,u_face(cell->Face(0)),u_face(cell->Face(1)),TV(x,2*y));}

    static T Interpolate_Y_Face_Transformed(const QUADTREE_GRID<T>& grid,const QUADTREE_CELL<T>* cell,const ARRAY<T>& u_face,const ARRAY<TV>& u_node,const T& x,const T& y)
    {if(x>=(T).5) return LINEAR_INTERPOLATION<T,T>::Bilinear(u_face(cell->Face(2)),u_node(cell->Node(1)).y,u_face(cell->Face(3)),u_node(cell->Node(3)).y,TV(2*(x-(T).5),y));
    else return LINEAR_INTERPOLATION<T,T>::Bilinear(u_node(cell->Node(0)).y,u_face(cell->Face(2)),u_node(cell->Node(2)).y,u_face(cell->Face(3)),TV(2*x,y));}
public:

     // assumes that this is the +x triangle where the distance between node1 and node2 is 2, and center is at the origin
    template<class T2> static T2 Interpolate_Triangle_Helper(const T2& center,const T2& node1,const T2& node2,const T x,const T y)
    {return (1-x)*center+(T).5*((x-y)*node1+(x+y)*node2);}

    template<class T2> static T2 Interpolate_Cells(const QUADTREE_GRID<T>& grid,const QUADTREE_CELL<T>* cell,const ARRAY<T2>& u_cell,const ARRAY<T2>& u_node,const TV& X)
    {assert(cell);TV DX_over_two((T).5*cell->DX());
    T x=(X.x-cell->Center().x)/DX_over_two.x,y=(X.y-cell->Center().y)/DX_over_two.y; // normalize the coordinates to be between -1 and 1
    int node_index1,node_index2;T x_prime,y_prime;
    if(abs(x)>abs(y)){
        if(x>0){node_index1=1;node_index2=3;x_prime=x;y_prime=y;} // +x pyramid
        else{node_index1=0;node_index2=2;x_prime=-x;y_prime=y;}} // -x pyramid
    else{
        if(y>0){node_index1=2;node_index2=3;x_prime=y;y_prime=x;} // +y pyramid
        else{node_index1=0;node_index2=1;x_prime=-y;y_prime=x;}} // -y pyramid
    return Interpolate_Triangle_Helper(u_cell(cell->Cell()),u_node(cell->Node(node_index1)),u_node(cell->Node(node_index2)),x_prime,y_prime);} // do the calculation as if it was the +x pyramid

//#####################################################################
};
}
#endif
#endif
