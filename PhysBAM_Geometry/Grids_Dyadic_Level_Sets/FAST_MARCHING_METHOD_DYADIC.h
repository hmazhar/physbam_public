//#####################################################################
// Copyright 2005, Frank Losasso, Andrew Selle.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class FAST_MARCHING_METHOD_DYADIC  
//#####################################################################
#ifndef COMPILE_WITHOUT_DYADIC_SUPPORT
#ifndef __FAST_MARCHING_METHOD_DYADIC__
#define __FAST_MARCHING_METHOD_DYADIC__

#include <PhysBAM_Geometry/Grids_Dyadic_Level_Sets/LEVELSET_DYADIC.h>
namespace PhysBAM{

template<class T_GRID>
class FAST_MARCHING_METHOD_DYADIC:public NONCOPYABLE
{
    typedef typename T_GRID::VECTOR_T TV;typedef typename TV::SCALAR T;
    typedef typename T_GRID::CELL CELL;typedef typename LEVELSET_POLICY<T_GRID>::LEVELSET T_LEVELSET;
    typedef typename T_GRID::CELL_ITERATOR CELL_ITERATOR;typedef typename T_GRID::NODE_ITERATOR NODE_ITERATOR;
public:
    T_GRID& grid;
    const LEVELSET_DYADIC<T_GRID>& levelset;
    bool currently_using_cells;
    ARRAY<CELL*>& cell_pointer_from_index;
    ARRAY<VECTOR<int,T_GRID::number_of_neighbors_per_node> >& node_neighbors;
    ARRAY<VECTOR<CELL*,T_GRID::number_of_neighbors_per_cell> >& cell_neighbors;

    const LEVELSET_DYADIC<T_GRID>* levelset_for_initialize;
    ARRAY<T>* phi_nodes_for_initialize;

    FAST_MARCHING_METHOD_DYADIC(const LEVELSET_DYADIC<T_GRID>& levelset);
    ~FAST_MARCHING_METHOD_DYADIC();

private:
    int Neighbor(const int neighbor_number,const int current_index)
    {if(currently_using_cells){
        CELL* neighbor=cell_neighbors(current_index)(neighbor_number);
        if(!neighbor||neighbor->Depth_Of_This_Cell()<grid.maximum_depth) return 0;
        return neighbor->Cell();}
    else return node_neighbors(current_index)(neighbor_number);}

    bool Neighbor_Visible(const int neighbor_number,const int current_index)
    {static const int opposite_face_when_used[]={0,1,0,2,0,3,0};
    if(levelset.collision_body_list){
        if(currently_using_cells){
            if(opposite_face_when_used[neighbor_number])
                return levelset.collision_body_list->cell_neighbors_visible(cell_neighbors(current_index)(neighbor_number)->Cell())(opposite_face_when_used[neighbor_number]);
            else return levelset.collision_body_list->cell_neighbors_visible(current_index)(opposite_face_when_used[neighbor_number-1]);}}
    return true;} // always return true for nodes since we won't do nodes where there are objects
    
    int Number_Of_Computational_Nodes()
    {if(currently_using_cells)return grid.number_of_cells;
    else return grid.number_of_nodes;}
    
    bool Use_Computational_Node(const int node)
    {if(currently_using_cells)return cell_pointer_from_index(node)->Depth_Of_This_Cell()==grid.maximum_depth;
    else return true;}

    TV Normal(const CELL* cell,const TV& X)
    {return levelset_for_initialize->Extended_Normal_From_Leaf_Cell(cell,X,phi_nodes_for_initialize);}

    T Phi(const CELL* cell,const TV& X)
    {return levelset_for_initialize->Phi_From_Close_Cell(cell,X,phi_nodes_for_initialize);}

    TV Computational_Node_Location(const int node)
    {if(currently_using_cells)return cell_pointer_from_index(node)->Center();
    else return grid.Node_Location(node);}

//#####################################################################
public:
    void Fast_Marching_Method_Cells(ARRAY<T>& phi_ghost,const T stopping_distance=0,const ARRAY<int>* seed_indices=0,const bool add_seed_indices_for_ghost_cells=false);
    void Fast_Marching_Method_Nodes(ARRAY<T>& phi_ghost_nodes,const T stopping_distance,const ARRAY<int>* seed_indices,const bool add_seed_indices_for_ghost_cells,ARRAY<T>* phi_ghost_cells=0);
private:
    void Fast_Marching_Method(ARRAY<T>& phi_ghost,const T stopping_distance=0,const ARRAY<int>* seed_indices=0,const bool add_seed_indices_for_ghost_cells=false);
    void Update_Or_Add_Neighbor(ARRAY<T>& phi_ghost,ARRAY<bool,VECTOR<int,1> >& done,ARRAY<int>& close_k,ARRAY<int>& heap,int& heap_length,int neighbor);
    void Initialize_Interface(ARRAY<T>& phi,ARRAY<bool,VECTOR<int,1> >& done,ARRAY<int>& close_k,ARRAY<int>& heap,int& heap_length,const ARRAY<int>* seed_indices=0,const bool add_seed_indices_for_ghost_cells=false);
    void Update_Close_Point(ARRAY<T>& phi,const ARRAY<bool,VECTOR<int,1> >& done,const int i);
    void Add_To_Initial(ARRAY<bool,VECTOR<int,1> >& done,ARRAY<int>& close_k,const int index);
//#####################################################################
};
}
#endif
#endif
