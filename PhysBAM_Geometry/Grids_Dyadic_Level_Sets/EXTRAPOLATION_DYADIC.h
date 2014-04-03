//#####################################################################
// Copyright 2003-2005, Ron Fedkiw, Frank Losasso.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class EXTRAPOLATION_DYADIC  
//##################################################################### 
#ifndef COMPILE_WITHOUT_DYADIC_SUPPORT
#ifndef __EXTRAPOLATION_DYADIC__
#define __EXTRAPOLATION_DYADIC__

#include <PhysBAM_Tools/Grids_Dyadic/OCTREE_GRID.h>
#include <PhysBAM_Tools/Grids_Dyadic/QUADTREE_GRID.h>
#include <PhysBAM_Geometry/Grids_Dyadic_Level_Sets/LEVELSET_OCTREE.h>
#include <PhysBAM_Geometry/Grids_Dyadic_Level_Sets/LEVELSET_QUADTREE.h>
#include <PhysBAM_Geometry/Level_Sets/EXTRAPOLATION.h>
namespace PhysBAM{

template<class T_GRID,class T2>
class EXTRAPOLATION_DYADIC:public EXTRAPOLATION<T_GRID,T2>
{
    typedef typename T_GRID::VECTOR_T TV;typedef typename TV::SCALAR T;
    typedef typename T_GRID::CELL CELL;typedef typename LEVELSET_POLICY<T_GRID>::LEVELSET T_LEVELSET;
    typedef typename T_GRID::CELL_ITERATOR CELL_ITERATOR;typedef typename T_GRID::NODE_ITERATOR NODE_ITERATOR;typedef typename T_GRID::FACE_ITERATOR FACE_ITERATOR;
    typedef typename INTERPOLATION_POLICY<T_GRID>::LINEAR_INTERPOLATION_HELPER T_LINEAR_INTERPOLATION_DYADIC_HELPER;
    typedef typename BOUNDARY_POLICY<T_GRID>::BOUNDARY_SCALAR T_BOUNDARY_SCALAR;typedef typename REBIND<T_BOUNDARY_SCALAR,T2>::TYPE T_BOUNDARY_EXTRAPOLATE;
public:
    typedef EXTRAPOLATION<T_GRID,T2> BASE;
    using BASE::band_width;using BASE::boundary;using BASE::isobaric_fix_width;
    
    T_GRID& grid;
    bool currently_using_cells;
    bool currently_using_faces;
    ARRAY<CELL*>& cell_pointer_from_index;
    ARRAY<VECTOR<int,T_GRID::number_of_neighbors_per_cell> >& node_neighbors;
    ARRAY<VECTOR<CELL*,T_GRID::number_of_neighbors_per_cell> >& cell_neighbors;
    ARRAY<int>* seed_indices;
    ARRAY<bool>* seed_done;
    bool collision_aware_extrapolation;
    const ARRAY<VECTOR<bool,T_GRID::dimension> >* neighbors_visible; // used for collision aware extrapolation
    T_BOUNDARY_SCALAR* boundary_scalar;
    T_BOUNDARY_EXTRAPOLATE* boundary_extrapolate;
    T small_number;
private:
    T tolerance; // optimizations
    FACE_ITERATOR face_iterator_for_neighbors;
    T_BOUNDARY_SCALAR& boundary_scalar_default;
    T_BOUNDARY_EXTRAPOLATE& boundary_extrapolate_default;
public:
 
    EXTRAPOLATION_DYADIC(T_GRID& grid_input);
    ~EXTRAPOLATION_DYADIC();

    void Set_Band_Width(const T number_of_cells=(T)5)
    {band_width=number_of_cells*grid.Minimum_Edge_Length();}
    
    void Set_Small_Number(const T small_number_input=1e-8)
    {small_number=small_number_input;}
    
    void Set_Band_Width_To_Entire_Domain()
    {band_width=1e10;}
    
    void Set_Isobaric_Fix_Width(const int number_of_cells=0)
    {isobaric_fix_width=number_of_cells*grid.Minimum_Edge_Length();}

    void Set_Custom_Seed_Indices(ARRAY<int>* seed_indices_input=0)
    {seed_indices=seed_indices_input;}

    void Set_Custom_Seed_Done(ARRAY<bool>* seed_done_input=0)
    {seed_done=seed_done_input;}

    void Set_Collision_Aware_Extrapolation(const ARRAY<VECTOR<bool,T_GRID::dimension> >& neighbors_visible_input)
    {collision_aware_extrapolation=true;neighbors_visible=&neighbors_visible_input;}

private:
    int Neighbor(const int neighbor_number,const int current_index)
    {if(currently_using_cells){
        CELL* neighbor=cell_neighbors(current_index)(neighbor_number);
        if(!neighbor||neighbor->Depth_Of_This_Cell()<grid.maximum_depth) return 0;
        return neighbor->Cell();}
    else if(currently_using_faces){
        face_iterator_for_neighbors.Set_Index(current_index);
        return face_iterator_for_neighbors.Minimal_Neighbor(neighbor_number-1);}
    else return node_neighbors(current_index)(neighbor_number);}

    bool Neighbor_Visible(const int neighbor_number,const int current_index)
    {static const int opposite_face_when_used[]={0,1,0,2,0,3,0};
    if(neighbors_visible){
        if(currently_using_cells){assert(neighbors_visible->m==grid.number_of_nodes);
            if(current_index>cell_neighbors.m) return true;CELL* cell_neighbor=cell_neighbors(current_index)(neighbor_number);if(cell_neighbor==0) return true;// temporary cure
            if(opposite_face_when_used[neighbor_number]) return (*neighbors_visible)(cell_neighbors(current_index)(neighbor_number)->Cell())(opposite_face_when_used[neighbor_number]);
            else return (*neighbors_visible)(current_index)(opposite_face_when_used[neighbor_number-1]);}
        else if(currently_using_faces){assert(neighbors_visible->m==grid.number_of_faces);
            if(current_index>cell_neighbors.m) return true;CELL* cell_neighbor=cell_neighbors(current_index)(neighbor_number);if(cell_neighbor==0) return true;// temporary cure
            if(opposite_face_when_used[neighbor_number]) return (*neighbors_visible)(cell_neighbors(current_index)(neighbor_number)->Cell())(opposite_face_when_used[neighbor_number]);
            else return (*neighbors_visible)(current_index)(opposite_face_when_used[neighbor_number-1]);}}
    return true;} // always return true for nodes since we won't do nodes where there are objects

    int Number_Of_Computational_Nodes()
    {if(currently_using_cells)return grid.number_of_cells;
    else if(currently_using_faces)return grid.number_of_faces;
    else return grid.number_of_nodes;}
    
    bool Use_Computational_Node(const int node)
    {if(currently_using_cells)return cell_pointer_from_index(node)->Depth_Of_This_Cell()==grid.maximum_depth;
    else if(currently_using_faces)return true;
    else return true;}

    TV Computational_Node_Location(const int node)
    {if(currently_using_cells)return cell_pointer_from_index(node)->Center();
    else if(currently_using_faces)return grid.Face_Location(node);
    else return grid.Node_Location(node);}

//#####################################################################
public:
    void Extrapolate_Nodes(ARRAY<T> phi_nodes,ARRAY<T2>& u,const bool zero_outside_band=false,const T time=0);
    void Extrapolate_Faces(ARRAY<T> phi_faces,ARRAY<T>& u,const bool zero_outside_band=false,const T time=0,const bool extrapolate_nodes=false);
    void Extrapolate_Cells(ARRAY<T> phi,ARRAY<T2>& u,const bool zero_outside_band=false,const T time=0,const bool extrapolate_nodes=false);
private:
    template<class T3> void Extrapolate(ARRAY<T>& phi_ghost,ARRAY<T3>& u_ghost,const T time);
    void Initialize(const ARRAY<T>& phi,ARRAY<bool,VECTOR<int,1> >& done,ARRAY<bool>& close,ARRAY<int>& heap,int& heap_length);
    template<class T3> void Update_Close_Point(ARRAY<T3>& u,const ARRAY<T>& phi,const ARRAY<bool,VECTOR<int,1> >& done,const int i);
//#####################################################################
};   
}
#endif
#endif
