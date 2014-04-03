//#####################################################################
// Copyright 2011, Mridul Aanjaneya, Jon Gretarsson.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class LINEAR_INTERPOLATION_BINTREE_HELPER
//#####################################################################
#if !COMPILE_WITHOUT_DYADIC_SUPPORT || COMPILE_WITH_BINTREE_SUPPORT
#ifndef __LINEAR_INTERPOLATION_BINTREE_HELPER__
#define __LINEAR_INTERPOLATION_BINTREE_HELPER__

#include <PhysBAM_Tools/Grids_Dyadic/MAP_BINTREE_MESH.h>
#include <PhysBAM_Tools/Grids_Dyadic_Interpolation/LINEAR_INTERPOLATION_DYADIC_HELPER.h>
#include <PhysBAM_Tools/Interpolation/LINEAR_INTERPOLATION.h>
namespace PhysBAM{

template<class T_input>
class LINEAR_INTERPOLATION_BINTREE_HELPER:public LINEAR_INTERPOLATION_DYADIC_HELPER<BINTREE_GRID<T_input> >
{
    typedef T_input T;typedef VECTOR<T,1> TV;
    typedef typename BINTREE_GRID<T_input>::INDEX INDEX;
private:
    LINEAR_INTERPOLATION_BINTREE_HELPER(); // disallow construction
public:
    template<class T2> static T2 Interpolate_Nodes(const BINTREE_GRID<T>& grid,const ARRAY<T2>& u,const TV& X)
    {return Interpolate_Nodes(grid,grid.Leaf_Cell(X),u,X);}

    template<class T2> static T2 Interpolate_Nodes(const BINTREE_GRID<T>& grid,const BINTREE_CELL<T>* cell,const ARRAY<T2>& u,const TV& X)
    {assert(cell); // this means we are outside the domain (of everything... including the ghost cells)
    TV center(cell->Center()),DX_over_two((T).5*cell->DX());
    TV minimum_corner(center-DX_over_two),maximum_corner(center+DX_over_two);
    return LINEAR_INTERPOLATION<T,T2>::Linear(minimum_corner.x,maximum_corner.x,u(cell->Node(0)),u(cell->Node(1)),X.x);}

    static TV Interpolate_Faces(const BINTREE_GRID<T>& grid,const BINTREE_CELL<T>* cell,const ARRAY<T>& u_face,const ARRAY<TV>& u_node,const TV& X)
    {return TV(Interpolate_Face(grid,1,cell,u_face,u_node,X));}

    static T Interpolate_Face(const BINTREE_GRID<T>& grid,const int axis,const BINTREE_CELL<T>* cell,const ARRAY<T>& u_face,const ARRAY<TV>& u_node,const TV& X)
    {assert(cell);TV center(cell->Center()),DX_over_two((T).5*cell->DX());
    TV minimum_corner(center-DX_over_two),maximum_corner(center+DX_over_two);
    return LINEAR_INTERPOLATION<T,T>::Linear(minimum_corner.x,maximum_corner.x,u_face(cell->Face(0)),u_face(cell->Face(1)),X.x);}

    template<class T2> static T2 Interpolate_Cells(const BINTREE_GRID<T>& grid,const BINTREE_CELL<T>* cell,const ARRAY<T2>& u_cell,const ARRAY<T2>& u_node,const TV& X)
    {assert(cell);TV DX_over_two((T).5*cell->DX());
     if((X.x-cell->Center().x)>0){
        BINTREE_CELL<T>* neighbor_cell=cell->Get_Neighbor(1,&grid,0,false);
        return LINEAR_INTERPOLATION<T,T2>::Linear(cell->Center().x,neighbor_cell->Center().x,u_cell(cell->Cell()),u_cell(neighbor_cell->Cell()),X.x);}
     else{
        BINTREE_CELL<T>* neighbor_cell=cell->Get_Neighbor(-1,&grid,0,false);
        return LINEAR_INTERPOLATION<T,T2>::Linear(neighbor_cell->Center().x,cell->Center().x,u_cell(neighbor_cell->Cell()),u_cell(cell->Cell()),X.x);}}

    template<class T2> static ARRAY<PAIR<INDEX,T> > Interpolate_Cell_Weights(const BINTREE_GRID<T>& grid,const ARRAY<T2>& u_cell,const TV& X)
    {
        ARRAY<PAIR<INDEX,T> > weights;
        BINTREE_CELL<T>* cell=grid.Leaf_Cell(X);
        if((X.x-cell->Center().x)>0){
            BINTREE_CELL<T>* neighbor_cell=cell->Get_Neighbor(1,&grid,0,false);
            T width=neighbor_cell->Center().x-cell->Center().x;
            T left_distance=(X.x-cell->Center().x)/width;
            weights.Append(PAIR<INDEX,T>(cell->Cell(),(T)1-left_distance));
            weights.Append(PAIR<INDEX,T>(neighbor_cell->Cell(),left_distance));}
        else{
            BINTREE_CELL<T>* neighbor_cell=cell->Get_Neighbor(-1,&grid,0,false);
            T width=cell->Center().x-neighbor_cell->Center().x;
            T left_distance=(X.x-neighbor_cell->Center().x)/width;
            weights.Append(PAIR<INDEX,T>(neighbor_cell->Cell(),(T)1-left_distance));
            weights.Append(PAIR<INDEX,T>(cell->Cell(),left_distance));}
        return weights;
    }
//#####################################################################
};
}
#endif
#endif
