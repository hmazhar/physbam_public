//#####################################################################
// Copyright 2010, Jon Gretarsson.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#if !COMPILE_WITHOUT_DYADIC_SUPPORT || COMPILE_WITH_BINTREE_SUPPORT
#ifndef __MAP_BINTREE_MESH__
#define __MAP_BINTREE_MESH__

#include <PhysBAM_Tools/Grids_Dyadic/BINTREE_CELL.h>
#include <PhysBAM_Tools/Grids_Uniform/GRID.h>
#include <climits>
namespace PhysBAM{

template<class T>
class MAP_BINTREE_MESH
{
    typedef VECTOR<T,1> TV;
private:
    MAP_BINTREE_MESH()
    {}

public:
    // Used to access the face that is being iterated over
    //  Usage:  int face=face_index[axis (of the face)][deepest cell (0 or 1 depending on which of the two cells touching the face is finer)]
    static const int face_by_axis[2];
    static const int nodes_by_axis[][2][1];
    static const int opposite_face[];

    static void Get_Deepest_Cell(const int number_of_cells,const BINTREE_CELL<T>** cells,int& deepest_cell,int& maximum_depth)
    {deepest_cell=0;maximum_depth=0;
    for(int i=0;i<number_of_cells;i++) if(cells[i]!=0){int depth=cells[i]->Depth_Of_This_Cell();if(depth>maximum_depth){maximum_depth=depth;deepest_cell=i;}}
    assert(maximum_depth);}

    static int Get_Deepest_Cell_Of_Two(const BINTREE_CELL<T>** cells,const int cell1,const int cell2)
    {if(cells[cell1]){
        if(cells[cell2]){int depth1=Depth_Of_This_Cell(cells[cell1]),depth2=Depth_Of_This_Cell(cells[cell2]);if(depth1>=depth2)return cell1;else return cell2;}
        else return cell1;}
    else if(cells[cell2]) return cell2;
    else return -1;}

    static void Get_Deepest_And_Shallowest_Cell(const int number_of_cells,const BINTREE_CELL<T>** cells,int& deepest_cell,int& maximum_depth,int& shallowest_cell,int& minimum_depth)
    {deepest_cell=0;maximum_depth=0;shallowest_cell=0;minimum_depth=INT_MAX;
    for(int i=0;i<number_of_cells;i++) if(cells[i]!=0){
        int depth=cells[i]->Depth_Of_This_Cell();
        if(depth>maximum_depth){maximum_depth=depth;deepest_cell=i;}
        if(depth<minimum_depth){minimum_depth=depth;shallowest_cell=i;}}
    assert(maximum_depth);assert(minimum_depth!=INT_MAX);assert(minimum_depth<=maximum_depth);}

    static const BINTREE_CELL<T>* Child(const int child,const BINTREE_CELL<T>* cell)
    {if(cell == 0) return 0;else if(cell->children == 0) return cell;else return cell->Child(child);}

    static bool Has_Children(const BINTREE_CELL<T>* cell)
    {if(cell == 0) return false;else return cell->Has_Children();}

    static int Depth_Of_This_Cell(const BINTREE_CELL<T>* cell)
    {if(cell == 0) return -1;else return cell->Depth_Of_This_Cell();}

    static BINTREE_CELL<T>* Get_Cell(const ARRAY<BINTREE_CELL<T>*,VECTOR<int,1> >& cells,const int i)
    {if(cells.Valid_Index(i)) return(cells(i));else return 0;}

//#####################################################################
public:
    static void Map_Nodes(const GRID<TV>& grid,const ARRAY<BINTREE_CELL<T>*,VECTOR<int,1> >& cells,const int number_of_ghost_cells,void* pointer,
        void (*map_function)(void* pointer,const BINTREE_CELL<T>* c1,const BINTREE_CELL<T>* c2));
    static void Map_Nodes(const int m_start,const int m_end,const ARRAY<BINTREE_CELL<T>*,VECTOR<int,1> >& cells,void* pointer,
        void (*map_function)(void* pointer,const BINTREE_CELL<T>* c1,const BINTREE_CELL<T>* c2));

    static void Map_Faces(const GRID<TV>& grid,const ARRAY<BINTREE_CELL<T>*,VECTOR<int,1> >& cells,const int number_of_ghost_cells,void* pointer,
        void (*map_function)(void* pointer,const BINTREE_CELL<T>* c1,const BINTREE_CELL<T>* c2,int axis));
    static void Map_Faces(const int m_start,const int m_end,const ARRAY<BINTREE_CELL<T>*,VECTOR<int,1> >& cells,void* pointer,
        void (*map_function)(void* pointer,const BINTREE_CELL<T>* c1,const BINTREE_CELL<T>* c2,int axis));

    template<class T3> static void Map_Faces_For_Reduction(const GRID<TV>& grid,const ARRAY<BINTREE_CELL<T>*,VECTOR<int,1> >& cells,const int number_of_ghost_cells,void* pointer,ARRAY<T3>& face_values,
        T3 (*map_function)(void* pointer,const T3 face_value));

    static void Map_Cells(const GRID<TV>& grid,const ARRAY<BINTREE_CELL<T>*,VECTOR<int,1> >& cells,const int number_of_ghost_cells,void* pointer,
        void (*map_function)(void* pointer,const BINTREE_CELL<T>* c1));
    static void Map_Cells(const int m_start,const int m_end,const ARRAY<BINTREE_CELL<T>*,VECTOR<int,1> >& cells,void* pointer,
        void (*map_function)(void* pointer,const BINTREE_CELL<T>* c1));

    // only map the nodes touching both a normal cell and a ghost cell
    static void Map_Boundary_Nodes(const GRID<TV>& grid,const ARRAY<BINTREE_CELL<T>*,VECTOR<int,1> >& cells,void* pointer,
        void (*map_function)(void* pointer,const BINTREE_CELL<T>* c1,const BINTREE_CELL<T>* c2));
    static void Map_Left_Boundary_Nodes(const GRID<TV>& grid,const ARRAY<BINTREE_CELL<T>*,VECTOR<int,1> >& cells,void* pointer,
        void (*map_function)(void* pointer,const BINTREE_CELL<T>* c1,const BINTREE_CELL<T>* c2));
    static void Map_Right_Boundary_Nodes(const GRID<TV>& grid,const ARRAY<BINTREE_CELL<T>*,VECTOR<int,1> >& cells,void* pointer,
        void (*map_function)(void* pointer,const BINTREE_CELL<T>* c1,const BINTREE_CELL<T>* c2));
    static void Map_Ghost_Nodes(const GRID<TV>& grid,const ARRAY<BINTREE_CELL<T>*,VECTOR<int,1> >& cells,const int number_of_ghost_cells,void* pointer,
        void (*map_function)(void* pointer,const BINTREE_CELL<T>* c1,const BINTREE_CELL<T>* c2));
    static void Map_Left_Ghost_Nodes(const GRID<TV>& grid,const ARRAY<BINTREE_CELL<T>*,VECTOR<int,1> >& cells,const int number_of_ghost_cells,void* pointer,
        void (*map_function)(void* pointer,const BINTREE_CELL<T>* c1,const BINTREE_CELL<T>* c2));
    static void Map_Right_Ghost_Nodes(const GRID<TV>& grid,const ARRAY<BINTREE_CELL<T>*,VECTOR<int,1> >& cells,const int number_of_ghost_cells,void* pointer,
        void (*map_function)(void* pointer,const BINTREE_CELL<T>* c1,const BINTREE_CELL<T>* c2));
    static void Map_Exterior_Ghost_Cell_Nodes(const GRID<TV>& grid,const ARRAY<BINTREE_CELL<T>*,VECTOR<int,1> >& cells,const int number_of_ghost_cells,void* pointer,
        void (*map_function)(void* pointer,const BINTREE_CELL<T>* c1,const BINTREE_CELL<T>* c2));
    // only map the faces touching both a normal cell and a ghost cell
    static void Map_Boundary_Faces(const GRID<TV>& grid,const ARRAY<BINTREE_CELL<T>*,VECTOR<int,1> >& cells,void* pointer,
        void (*map_function)(void* pointer,const BINTREE_CELL<T>* c1,const BINTREE_CELL<T>* c2,int axis));
    static void Map_Left_Boundary_Faces(const GRID<TV>& grid,const ARRAY<BINTREE_CELL<T>*,VECTOR<int,1> >& cells,void* pointer,
        void (*map_function)(void* pointer,const BINTREE_CELL<T>* c1,const BINTREE_CELL<T>* c2,int axis));
    static void Map_Right_Boundary_Faces(const GRID<TV>& grid,const ARRAY<BINTREE_CELL<T>*,VECTOR<int,1> >& cells,void* pointer,
        void (*map_function)(void* pointer,const BINTREE_CELL<T>* c1,const BINTREE_CELL<T>* c2,int axis));
    static void Map_Domain_Side_Faces(const GRID<TV>& grid,const ARRAY<BINTREE_CELL<T>*,VECTOR<int,1> >& cells,const int side,const int number_of_ghost_cells,void* pointer,
        void (*map_function)(void* pointer,const BINTREE_CELL<T>* c1,const BINTREE_CELL<T>* c2,int axis));
    static void Map_Xmin_Faces(const GRID<TV>& grid,const ARRAY<BINTREE_CELL<T>*,VECTOR<int,1> >& cells,const int number_of_ghost_cells,void* pointer,
        void (*map_function)(void* pointer,const BINTREE_CELL<T>* c1,const BINTREE_CELL<T>* c2,int axis));
    static void Map_Xmax_Faces(const GRID<TV>& grid,const ARRAY<BINTREE_CELL<T>*,VECTOR<int,1> >& cells,const int number_of_ghost_cells,void* pointer,
        void (*map_function)(void* pointer,const BINTREE_CELL<T>* c1,const BINTREE_CELL<T>* c2,int axis));
    static void Map_Ghost_Faces(const GRID<TV>& grid,const ARRAY<BINTREE_CELL<T>*,VECTOR<int,1> >& cells,const int number_of_ghost_cells,void* pointer,
        void (*map_function)(void* pointer,const BINTREE_CELL<T>* c1,const BINTREE_CELL<T>* c2,int axis));
    static void Map_Exterior_Ghost_Cell_Faces(const GRID<TV>& grid,const ARRAY<BINTREE_CELL<T>*,VECTOR<int,1> >& cells,const int number_of_ghost_cells,void* pointer,
        void (*map_function)(void* pointer,const BINTREE_CELL<T>* c1,const BINTREE_CELL<T>* c2,int axis));
    // map the cells that are surrounding the area. (maps one coarse cell around the entire boundary)
    static void Map_Ghost_Cells(const GRID<TV>& grid,const ARRAY<BINTREE_CELL<T>*,VECTOR<int,1> >& cells,const int number_of_ghost_cells,void* pointer,
        void (*map_function)(void* pointer,const BINTREE_CELL<T>* c1));
    static void Map_Boundary_Cells(const GRID<TV>& grid,const ARRAY<BINTREE_CELL<T>*,VECTOR<int,1> >& cells,void* pointer,
        void (*map_function)(void* pointer,const BINTREE_CELL<T>* c1));
    static void Map_Left_Boundary_Cells(const GRID<TV>& grid,const ARRAY<BINTREE_CELL<T>*,VECTOR<int,1> >& cells,void* pointer,
        void (*map_function)(void* pointer,const BINTREE_CELL<T>* c1));
    static void Map_Right_Boundary_Cells(const GRID<TV>& grid,const ARRAY<BINTREE_CELL<T>*,VECTOR<int,1> >& cells,void* pointer,
        void (*map_function)(void* pointer,const BINTREE_CELL<T>* c1));
    // these are used mostly internally in this class, but some things outside may use them if they are careful about it
    static void Map_Node_Process_Cell(void* pointer,const BINTREE_CELL<T>* cell,
        void (*map_function)(void* pointer,const BINTREE_CELL<T>* c1,const BINTREE_CELL<T>* c2));
    static void Map_Node_Process_Node(void* pointer,const BINTREE_CELL<T>* c1,const BINTREE_CELL<T>* c2,
        void (*map_function)(void* pointer,const BINTREE_CELL<T>* c1,const BINTREE_CELL<T>* c2));
    static void Map_Face_Process_Cell(void* pointer,const BINTREE_CELL<T>* cell,
        void (*map_function)(void* pointer,const BINTREE_CELL<T>* c1,const BINTREE_CELL<T>* c2,int axis));
    template<class T3> static void Map_Face_Process_Cell_For_Reduction(void* pointer,const BINTREE_CELL<T>* cell,ARRAY<T3>& face_values,T3 (*map_function)(void* pointer,const T3 face_value));
    static void Map_Face_Process_Face(void* pointer,const BINTREE_CELL<T>* c1,const BINTREE_CELL<T>* c2,
        void (*map_function)(void* pointer,const BINTREE_CELL<T>* c1,const BINTREE_CELL<T>* c2,int axis));
    template<class T3> static T3 Map_Face_Process_Face_For_Reduction(void* pointer,const BINTREE_CELL<T>* c1,const BINTREE_CELL<T>* c2,ARRAY<T3>& face_values,
        T3 (*map_function)(void* pointer,const T3 face_value));
    static void Map_Cell_Process_Cell(void* pointer,const BINTREE_CELL<T>* cell,
        void (*map_function)(void* pointer,const BINTREE_CELL<T>* c1));
//#####################################################################
};
}
#endif
#endif
