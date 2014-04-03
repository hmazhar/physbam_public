//#####################################################################
// Copyright 2003-2008, Ron Fedkiw, Frank Losasso, Andrew Selle.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#ifndef COMPILE_WITHOUT_DYADIC_SUPPORT
#ifndef __MAP_OCTREE_MESH__
#define __MAP_OCTREE_MESH__

#include <PhysBAM_Tools/Grids_Dyadic/OCTREE_CELL.h>
#include <PhysBAM_Tools/Grids_Uniform/GRID.h>
#include <climits>
namespace PhysBAM{

template<class T>
class MAP_OCTREE_MESH
{
    typedef VECTOR<T,3> TV;
private:
    MAP_OCTREE_MESH()
    {}

public:
    // Face iterator lookup functions:

    // Used to access the 4 nodes that touch the face being iterated over
    //  Usage:  int node=nodes_by_axis[axis (of the face)][deepest cell (0 or 1 depending on which of the two cells touching the face is finer)][0...3 (which node is wanted)]
    static const int nodes_by_axis[3][2][4];

    // Used to access the face that is being iterated over
    //  Usage:  int face=face_index[axis (of the face)][deepest cell (0 or 1 depending on which of the two cells touching the face is finer)]
    static const int face_by_axis[3][2];

    // Edge iterator lookup functions:

    // Used to access the 2 nodes that touch the edge being iterated over
    //  Usage:  int node=nodes_by_axis[axis (of the edge)][deepest cell (0...3 depending on which of the two cells touching the edge is finer)][0...2 (which node is wanted)]
    static const int edge_nodes_by_axis[3][4][2];

    // Misc lookup function:

    // Used to get the opposite face from the index
    static const int opposite_face[6];

    static void Get_Deepest_Cell(const int number_of_cells,const OCTREE_CELL<T>** cells,int& deepest_cell,int& maximum_depth)
    {deepest_cell=0;maximum_depth=0;
    for(int i=0;i<number_of_cells;i++) if(cells[i]!=0){int depth=cells[i]->Depth_Of_This_Cell();if(depth>maximum_depth){maximum_depth=depth;deepest_cell=i;}}
    assert(maximum_depth);}

    static int Get_Deepest_Cell_Of_Two(const OCTREE_CELL<T>** cells,const int cell1,const int cell2)
    {if(cells[cell1]){
        if(cells[cell2]){int depth1=Depth_Of_This_Cell(cells[cell1]),depth2=Depth_Of_This_Cell(cells[cell2]);if(depth1>=depth2)return cell1;else return cell2;}
        else return cell1;}
    else if(cells[cell2]) return cell2;
    else return -1;}

    static void Get_Deepest_And_Shallowest_Cell(const int number_of_cells,const OCTREE_CELL<T>** cells,int& deepest_cell,int& maximum_depth,int& shallowest_cell,int& minimum_depth)
    {deepest_cell=0;maximum_depth=0;shallowest_cell=0;minimum_depth=INT_MAX;
    for(int i=0;i<number_of_cells;i++) if(cells[i]!=0){
        int depth=cells[i]->Depth_Of_This_Cell();
        if(depth>maximum_depth){maximum_depth=depth;deepest_cell=i;}
        if(depth<minimum_depth){minimum_depth=depth;shallowest_cell=i;}}
    assert(maximum_depth);assert(minimum_depth!=INT_MAX);assert(minimum_depth<=maximum_depth);}

    static const OCTREE_CELL<T>* Child(const int child,const OCTREE_CELL<T>* cell)
    {if(cell == 0) return 0;else if(cell->children == 0) return cell;else return cell->Child(child);}

    static bool Has_Children(const OCTREE_CELL<T>* cell)
    {if(cell == 0) return false;else return cell->Has_Children();}

    static int Depth_Of_This_Cell(const OCTREE_CELL<T>* cell)
    {if(cell == 0) return -1;else return cell->Depth_Of_This_Cell();}

    static OCTREE_CELL<T>* Get_Cell(const ARRAY<OCTREE_CELL<T>*,VECTOR<int,3> >& cells,const int i,const int j,const int ij)
    {if(cells.Valid_Index(i,j,ij)) return(cells(i,j,ij));else return 0;}

//#####################################################################
public:
    static void Map_Nodes(const GRID<TV>& grid,const ARRAY<OCTREE_CELL<T>*,VECTOR<int,3> >& cells,const int number_of_ghost_cells,void* pointer,
        void (*map_function)(void* pointer,const OCTREE_CELL<T>* c1,const OCTREE_CELL<T>* c2,const OCTREE_CELL<T>* c3,const OCTREE_CELL<T>* c4,const OCTREE_CELL<T>* c5,const OCTREE_CELL<T>* c6,const OCTREE_CELL<T>* c7,const OCTREE_CELL<T>* c8));
    static void Map_Nodes(int m_start,int m_end,int n_start,int n_end,int mn_start,int mn_end,const ARRAY<OCTREE_CELL<T>*,VECTOR<int,3> >& cells,void* pointer,
        void (*map_function)(void* pointer,const OCTREE_CELL<T>* c1,const OCTREE_CELL<T>* c2,const OCTREE_CELL<T>* c3,const OCTREE_CELL<T>* c4,const OCTREE_CELL<T>* c5,const OCTREE_CELL<T>* c6,const OCTREE_CELL<T>* c7,const OCTREE_CELL<T>* c8));

    static void Map_Edges(const GRID<TV>& grid,const ARRAY<OCTREE_CELL<T>*,VECTOR<int,3> >& cells,const int number_of_ghost_cells,void* pointer,
        void (*map_function)(void* pointer,const OCTREE_CELL<T>* c1,const OCTREE_CELL<T>* c2,const OCTREE_CELL<T>* c3,const OCTREE_CELL<T>* c4,int axis));
    static void Map_Edges(int m_start,int m_end,int n_start,int n_end,int mn_start,int mn_end,const ARRAY<OCTREE_CELL<T>*,VECTOR<int,3> >& cells,void* pointer,
        void (*map_function)(void* pointer,const OCTREE_CELL<T>* c1,const OCTREE_CELL<T>* c2,const OCTREE_CELL<T>* c3,const OCTREE_CELL<T>* c4,int axis));

    static void Map_Faces(const GRID<TV>& grid,const ARRAY<OCTREE_CELL<T>*,VECTOR<int,3> >& cells,const int number_of_ghost_cells,void* pointer,
        void (*map_function)(void* pointer,const OCTREE_CELL<T>* c1,const OCTREE_CELL<T>* c2,int axis));
    static void Map_Faces(int m_start,int m_end,int n_start,int n_end,int mn_start,int mn_end,const ARRAY<OCTREE_CELL<T>*,VECTOR<int,3> >& cells,void* pointer,
        void (*map_function)(void* pointer,const OCTREE_CELL<T>* c1,const OCTREE_CELL<T>* c2,int axis));

    template<class T3> static void Map_Faces_For_Reduction(const GRID<TV>& grid,const ARRAY<OCTREE_CELL<T>*,VECTOR<int,3> >& cells,const int number_of_ghost_cells,void* pointer,ARRAY<T3>& face_values,
        T3 (*map_function)(void* pointer,const T3 face_value_1,const T3 face_value_2,const T3 face_value_3,const T3 face_value_4));

    static void Map_Cells(const GRID<TV>& grid,const ARRAY<OCTREE_CELL<T>*,VECTOR<int,3> >& cells,const int number_of_ghost_cells,void* pointer,
        void (*map_function)(void* pointer,const OCTREE_CELL<T>* c1));
    static void Map_Cells(int m_start,int m_end,int n_start,int n_end,int mn_start,int mn_end,const ARRAY<OCTREE_CELL<T>*,VECTOR<int,3> >& cells,void* pointer,
        void (*map_function)(void* pointer,const OCTREE_CELL<T>* c1));

    // only map the nodes touching both a normal cell and a ghost cell
    static void Map_Boundary_Nodes(const GRID<TV>& grid,const ARRAY<OCTREE_CELL<T>*,VECTOR<int,3> >& cells,void* pointer,
        void (*map_function)(void* pointer,const OCTREE_CELL<T>* c1,const OCTREE_CELL<T>* c2,const OCTREE_CELL<T>* c3,const OCTREE_CELL<T>* c4,const OCTREE_CELL<T>* c5,const OCTREE_CELL<T>* c6,const OCTREE_CELL<T>* c7,const OCTREE_CELL<T>* c8));
    static void Map_Left_Boundary_Nodes(const GRID<TV>& grid,const ARRAY<OCTREE_CELL<T>*,VECTOR<int,3> >& cells,void* pointer,
        void (*map_function)(void* pointer,const OCTREE_CELL<T>* c1,const OCTREE_CELL<T>* c2,const OCTREE_CELL<T>* c3,const OCTREE_CELL<T>* c4,const OCTREE_CELL<T>* c5,const OCTREE_CELL<T>* c6,const OCTREE_CELL<T>* c7,const OCTREE_CELL<T>* c8));
    static void Map_Right_Boundary_Nodes(const GRID<TV>& grid,const ARRAY<OCTREE_CELL<T>*,VECTOR<int,3> >& cells,void* pointer,
        void (*map_function)(void* pointer,const OCTREE_CELL<T>* c1,const OCTREE_CELL<T>* c2,const OCTREE_CELL<T>* c3,const OCTREE_CELL<T>* c4,const OCTREE_CELL<T>* c5,const OCTREE_CELL<T>* c6,const OCTREE_CELL<T>* c7,const OCTREE_CELL<T>* c8));
    static void Map_Bottom_Boundary_Nodes(const GRID<TV>& grid,const ARRAY<OCTREE_CELL<T>*,VECTOR<int,3> >& cells,void* pointer,
        void (*map_function)(void* pointer,const OCTREE_CELL<T>* c1,const OCTREE_CELL<T>* c2,const OCTREE_CELL<T>* c3,const OCTREE_CELL<T>* c4,const OCTREE_CELL<T>* c5,const OCTREE_CELL<T>* c6,const OCTREE_CELL<T>* c7,const OCTREE_CELL<T>* c8));
    static void Map_Top_Boundary_Nodes(const GRID<TV>& grid,const ARRAY<OCTREE_CELL<T>*,VECTOR<int,3> >& cells,void* pointer,
        void (*map_function)(void* pointer,const OCTREE_CELL<T>* c1,const OCTREE_CELL<T>* c2,const OCTREE_CELL<T>* c3,const OCTREE_CELL<T>* c4,const OCTREE_CELL<T>* c5,const OCTREE_CELL<T>* c6,const OCTREE_CELL<T>* c7,const OCTREE_CELL<T>* c8));
    static void Map_Front_Boundary_Nodes(const GRID<TV>& grid,const ARRAY<OCTREE_CELL<T>*,VECTOR<int,3> >& cells,void* pointer,
        void (*map_function)(void* pointer,const OCTREE_CELL<T>* c1,const OCTREE_CELL<T>* c2,const OCTREE_CELL<T>* c3,const OCTREE_CELL<T>* c4,const OCTREE_CELL<T>* c5,const OCTREE_CELL<T>* c6,const OCTREE_CELL<T>* c7,const OCTREE_CELL<T>* c8));
    static void Map_Back_Boundary_Nodes(const GRID<TV>& grid,const ARRAY<OCTREE_CELL<T>*,VECTOR<int,3> >& cells,void* pointer,
        void (*map_function)(void* pointer,const OCTREE_CELL<T>* c1,const OCTREE_CELL<T>* c2,const OCTREE_CELL<T>* c3,const OCTREE_CELL<T>* c4,const OCTREE_CELL<T>* c5,const OCTREE_CELL<T>* c6,const OCTREE_CELL<T>* c7,const OCTREE_CELL<T>* c8));
    static void Map_Ghost_Nodes(const GRID<TV>& grid,const ARRAY<OCTREE_CELL<T>*,VECTOR<int,3> >& cells,const int number_of_ghost_cells,void* pointer,
        void (*map_function)(void* pointer,const OCTREE_CELL<T>* c1,const OCTREE_CELL<T>* c2,const OCTREE_CELL<T>* c3,const OCTREE_CELL<T>* c4,const OCTREE_CELL<T>* c5,const OCTREE_CELL<T>* c6,const OCTREE_CELL<T>* c7,const OCTREE_CELL<T>* c8));
    static void Map_Left_Ghost_Nodes(const GRID<TV>& grid,const ARRAY<OCTREE_CELL<T>*,VECTOR<int,3> >& cells,const int number_of_ghost_cells,void* pointer,
        void (*map_function)(void* pointer,const OCTREE_CELL<T>* c1,const OCTREE_CELL<T>* c2,const OCTREE_CELL<T>* c3,const OCTREE_CELL<T>* c4,const OCTREE_CELL<T>* c5,const OCTREE_CELL<T>* c6,const OCTREE_CELL<T>* c7,const OCTREE_CELL<T>* c8));
    static void Map_Right_Ghost_Nodes(const GRID<TV>& grid,const ARRAY<OCTREE_CELL<T>*,VECTOR<int,3> >& cells,const int number_of_ghost_cells,void* pointer,
        void (*map_function)(void* pointer,const OCTREE_CELL<T>* c1,const OCTREE_CELL<T>* c2,const OCTREE_CELL<T>* c3,const OCTREE_CELL<T>* c4,const OCTREE_CELL<T>* c5,const OCTREE_CELL<T>* c6,const OCTREE_CELL<T>* c7,const OCTREE_CELL<T>* c8));
    static void Map_Bottom_Ghost_Nodes(const GRID<TV>& grid,const ARRAY<OCTREE_CELL<T>*,VECTOR<int,3> >& cells,const int number_of_ghost_cells,void* pointer,
        void (*map_function)(void* pointer,const OCTREE_CELL<T>* c1,const OCTREE_CELL<T>* c2,const OCTREE_CELL<T>* c3,const OCTREE_CELL<T>* c4,const OCTREE_CELL<T>* c5,const OCTREE_CELL<T>* c6,const OCTREE_CELL<T>* c7,const OCTREE_CELL<T>* c8));
    static void Map_Top_Ghost_Nodes(const GRID<TV>& grid,const ARRAY<OCTREE_CELL<T>*,VECTOR<int,3> >& cells,const int number_of_ghost_cells,void* pointer,
        void (*map_function)(void* pointer,const OCTREE_CELL<T>* c1,const OCTREE_CELL<T>* c2,const OCTREE_CELL<T>* c3,const OCTREE_CELL<T>* c4,const OCTREE_CELL<T>* c5,const OCTREE_CELL<T>* c6,const OCTREE_CELL<T>* c7,const OCTREE_CELL<T>* c8));
    static void Map_Front_Ghost_Nodes(const GRID<TV>& grid,const ARRAY<OCTREE_CELL<T>*,VECTOR<int,3> >& cells,const int number_of_ghost_cells,void* pointer,
        void (*map_function)(void* pointer,const OCTREE_CELL<T>* c1,const OCTREE_CELL<T>* c2,const OCTREE_CELL<T>* c3,const OCTREE_CELL<T>* c4,const OCTREE_CELL<T>* c5,const OCTREE_CELL<T>* c6,const OCTREE_CELL<T>* c7,const OCTREE_CELL<T>* c8));
    static void Map_Back_Ghost_Nodes(const GRID<TV>& grid,const ARRAY<OCTREE_CELL<T>*,VECTOR<int,3> >& cells,const int number_of_ghost_cells,void* pointer,
        void (*map_function)(void* pointer,const OCTREE_CELL<T>* c1,const OCTREE_CELL<T>* c2,const OCTREE_CELL<T>* c3,const OCTREE_CELL<T>* c4,const OCTREE_CELL<T>* c5,const OCTREE_CELL<T>* c6,const OCTREE_CELL<T>* c7,const OCTREE_CELL<T>* c8));
    static void Map_Exterior_Ghost_Cell_Nodes(const GRID<TV>& grid,const ARRAY<OCTREE_CELL<T>*,VECTOR<int,3> >& cells,const int number_of_ghost_cells,void* pointer,
        void (*map_function)(void* pointer,const OCTREE_CELL<T>* c1,const OCTREE_CELL<T>* c2,const OCTREE_CELL<T>* c3,const OCTREE_CELL<T>* c4,const OCTREE_CELL<T>* c5,const OCTREE_CELL<T>* c6,const OCTREE_CELL<T>* c7,const OCTREE_CELL<T>* c8));
    // only map the faces touching both a normal cell and a ghost cell
    static void Map_Boundary_Faces(const GRID<TV>& grid,const ARRAY<OCTREE_CELL<T>*,VECTOR<int,3> >& cells,void* pointer,
        void (*map_function)(void* pointer,const OCTREE_CELL<T>* c1,const OCTREE_CELL<T>* c2,int axis));
    static void Map_Left_Boundary_Faces(const GRID<TV>& grid,const ARRAY<OCTREE_CELL<T>*,VECTOR<int,3> >& cells,void* pointer,
        void (*map_function)(void* pointer,const OCTREE_CELL<T>* c1,const OCTREE_CELL<T>* c2,int axis));
    static void Map_Right_Boundary_Faces(const GRID<TV>& grid,const ARRAY<OCTREE_CELL<T>*,VECTOR<int,3> >& cells,void* pointer,
        void (*map_function)(void* pointer,const OCTREE_CELL<T>* c1,const OCTREE_CELL<T>* c2,int axis));
    static void Map_Bottom_Boundary_Faces(const GRID<TV>& grid,const ARRAY<OCTREE_CELL<T>*,VECTOR<int,3> >& cells,void* pointer,
        void (*map_function)(void* pointer,const OCTREE_CELL<T>* c1,const OCTREE_CELL<T>* c2,int axis));
    static void Map_Top_Boundary_Faces(const GRID<TV>& grid,const ARRAY<OCTREE_CELL<T>*,VECTOR<int,3> >& cells,void* pointer,
        void (*map_function)(void* pointer,const OCTREE_CELL<T>* c1,const OCTREE_CELL<T>* c2,int axis));
    static void Map_Front_Boundary_Faces(const GRID<TV>& grid,const ARRAY<OCTREE_CELL<T>*,VECTOR<int,3> >& cells,void* pointer,
        void (*map_function)(void* pointer,const OCTREE_CELL<T>* c1,const OCTREE_CELL<T>* c2,int axis));
    static void Map_Back_Boundary_Faces(const GRID<TV>& grid,const ARRAY<OCTREE_CELL<T>*,VECTOR<int,3> >& cells,void* pointer,
        void (*map_function)(void* pointer,const OCTREE_CELL<T>* c1,const OCTREE_CELL<T>* c2,int axis));
    static void Map_Domain_Side_Faces(const GRID<TV>& grid,const ARRAY<OCTREE_CELL<T>*,VECTOR<int,3> >& cells,const int side,const int number_of_ghost_cells,void* pointer,
        void (*map_function)(void* pointer,const OCTREE_CELL<T>* c1,const OCTREE_CELL<T>* c2,int axis));
    static void Map_Xmin_Faces(const GRID<TV>& grid,const ARRAY<OCTREE_CELL<T>*,VECTOR<int,3> >& cells,const int number_of_ghost_cells,void* pointer,
        void (*map_function)(void* pointer,const OCTREE_CELL<T>* c1,const OCTREE_CELL<T>* c2,int axis));
    static void Map_Xmax_Faces(const GRID<TV>& grid,const ARRAY<OCTREE_CELL<T>*,VECTOR<int,3> >& cells,const int number_of_ghost_cells,void* pointer,
        void (*map_function)(void* pointer,const OCTREE_CELL<T>* c1,const OCTREE_CELL<T>* c2,int axis));
    static void Map_Ymin_Faces(const GRID<TV>& grid,const ARRAY<OCTREE_CELL<T>*,VECTOR<int,3> >& cells,const int number_of_ghost_cells,void* pointer,
        void (*map_function)(void* pointer,const OCTREE_CELL<T>* c1,const OCTREE_CELL<T>* c2,int axis));
    static void Map_Ymax_Faces(const GRID<TV>& grid,const ARRAY<OCTREE_CELL<T>*,VECTOR<int,3> >& cells,const int number_of_ghost_cells,void* pointer,
        void (*map_function)(void* pointer,const OCTREE_CELL<T>* c1,const OCTREE_CELL<T>* c2,int axis));
    static void Map_Zmin_Faces(const GRID<TV>& grid,const ARRAY<OCTREE_CELL<T>*,VECTOR<int,3> >& cells,const int number_of_ghost_cells,void* pointer,
        void (*map_function)(void* pointer,const OCTREE_CELL<T>* c1,const OCTREE_CELL<T>* c2,int axis));
    static void Map_Zmax_Faces(const GRID<TV>& grid,const ARRAY<OCTREE_CELL<T>*,VECTOR<int,3> >& cells,const int number_of_ghost_cells,void* pointer,
        void (*map_function)(void* pointer,const OCTREE_CELL<T>* c1,const OCTREE_CELL<T>* c2,int axis));
    static void Map_Ghost_Faces(const GRID<TV>& grid,const ARRAY<OCTREE_CELL<T>*,VECTOR<int,3> >& cells,const int number_of_ghost_cells,void* pointer,
        void (*map_function)(void* pointer,const OCTREE_CELL<T>* c1,const OCTREE_CELL<T>* c2,int axis));
    static void Map_Exterior_Ghost_Cell_Faces(const GRID<TV>& grid,const ARRAY<OCTREE_CELL<T>*,VECTOR<int,3> >& cells,const int number_of_ghost_cells,void* pointer,
        void (*map_function)(void* pointer,const OCTREE_CELL<T>* c1,const OCTREE_CELL<T>* c2,int axis));
    // map the cells that are surrounding the area. (maps one coarse cell around the entire boundary)
    static void Map_Ghost_Cells(const GRID<TV>& grid,const ARRAY<OCTREE_CELL<T>*,VECTOR<int,3> >& cells,const int number_of_ghost_cells,void* pointer,
        void (*map_function)(void* pointer,const OCTREE_CELL<T>* c1));
    static void Map_Boundary_Cells(const GRID<TV>& grid,const ARRAY<OCTREE_CELL<T>*,VECTOR<int,3> >& cells,void* pointer,
        void (*map_function)(void* pointer,const OCTREE_CELL<T>* c1));
    static void Map_Left_Boundary_Cells(const GRID<TV>& grid,const ARRAY<OCTREE_CELL<T>*,VECTOR<int,3> >& cells,void* pointer,
        void (*map_function)(void* pointer,const OCTREE_CELL<T>* c1));
    static void Map_Right_Boundary_Cells(const GRID<TV>& grid,const ARRAY<OCTREE_CELL<T>*,VECTOR<int,3> >& cells,void* pointer,
        void (*map_function)(void* pointer,const OCTREE_CELL<T>* c1));
    static void Map_Bottom_Boundary_Cells(const GRID<TV>& grid,const ARRAY<OCTREE_CELL<T>*,VECTOR<int,3> >& cells,void* pointer,
        void (*map_function)(void* pointer,const OCTREE_CELL<T>* c1));
    static void Map_Top_Boundary_Cells(const GRID<TV>& grid,const ARRAY<OCTREE_CELL<T>*,VECTOR<int,3> >& cells,void* pointer,
        void (*map_function)(void* pointer,const OCTREE_CELL<T>* c1));
    static void Map_Front_Boundary_Cells(const GRID<TV>& grid,const ARRAY<OCTREE_CELL<T>*,VECTOR<int,3> >& cells,void* pointer,
        void (*map_function)(void* pointer,const OCTREE_CELL<T>* c1));
    static void Map_Back_Boundary_Cells(const GRID<TV>& grid,const ARRAY<OCTREE_CELL<T>*,VECTOR<int,3> >& cells,void* pointer,
        void (*map_function)(void* pointer,const OCTREE_CELL<T>* c1));

    // these are used mostly internally in this class, but some things outside may use them if they are careful about it
    static void Map_Node_Process_Cell(void* pointer,const OCTREE_CELL<T>* cell,
        void (*map_function)(void* pointer,const OCTREE_CELL<T>* c1,const OCTREE_CELL<T>* c2,const OCTREE_CELL<T>* c3,const OCTREE_CELL<T>* c4,const OCTREE_CELL<T>* c5,const OCTREE_CELL<T>* c6,const OCTREE_CELL<T>* c7,const OCTREE_CELL<T>* c8));
    static void Map_Node_Process_Face(void* pointer,const OCTREE_CELL<T>* c1,const OCTREE_CELL<T>* c2,const int orient,
        void (*map_function)(void* pointer,const OCTREE_CELL<T>* c1,const OCTREE_CELL<T>* c2,const OCTREE_CELL<T>* c3,const OCTREE_CELL<T>* c4,const OCTREE_CELL<T>* c5,const OCTREE_CELL<T>* c6,const OCTREE_CELL<T>* c7,const OCTREE_CELL<T>* c8));
    static void Map_Node_Process_Edge(void* pointer,const OCTREE_CELL<T>* c1,const OCTREE_CELL<T>* c2,const OCTREE_CELL<T>* c3,const OCTREE_CELL<T>* c4,const int orient,
        void (*map_function)(void* pointer,const OCTREE_CELL<T>* c1,const OCTREE_CELL<T>* c2,const OCTREE_CELL<T>* c3,const OCTREE_CELL<T>* c4,const OCTREE_CELL<T>* c5,const OCTREE_CELL<T>* c6,const OCTREE_CELL<T>* c7,const OCTREE_CELL<T>* c8));
    static void Map_Node_Process_Node(void* pointer,const OCTREE_CELL<T>* c1,const OCTREE_CELL<T>* c2,const OCTREE_CELL<T>* c3,const OCTREE_CELL<T>* c4,const OCTREE_CELL<T>* c5,const OCTREE_CELL<T>* c6,const OCTREE_CELL<T>* c7,const OCTREE_CELL<T>* c8,
        void (*map_function)(void* pointer,const OCTREE_CELL<T>* c1,const OCTREE_CELL<T>* c2,const OCTREE_CELL<T>* c3,const OCTREE_CELL<T>* c4,const OCTREE_CELL<T>* c5,const OCTREE_CELL<T>* c6,const OCTREE_CELL<T>* c7,const OCTREE_CELL<T>* c8));
    static void Map_Edge_Process_Cell(void* pointer,const OCTREE_CELL<T>* cell,
        void (*map_function)(void* pointer,const OCTREE_CELL<T>* c1,const OCTREE_CELL<T>* c2,const OCTREE_CELL<T>* c3,const OCTREE_CELL<T>* c4,int axis));
    static void Map_Edge_Process_Face(void* pointer,const OCTREE_CELL<T>* c1,const OCTREE_CELL<T>* c2,const int orient,
        void (*map_function)(void* pointer,const OCTREE_CELL<T>* c1,const OCTREE_CELL<T>* c2,const OCTREE_CELL<T>* c3,const OCTREE_CELL<T>* c4,int axis));
    static void Map_Edge_Process_Edge(void* pointer,const OCTREE_CELL<T>* c1,const OCTREE_CELL<T>* c2,const OCTREE_CELL<T>* c3,const OCTREE_CELL<T>* c4,const int orient,
        void (*map_function)(void* pointer,const OCTREE_CELL<T>* c1,const OCTREE_CELL<T>* c2,const OCTREE_CELL<T>* c3,const OCTREE_CELL<T>* c4,int axis));
    static void Map_Face_Process_Cell(void* pointer,const OCTREE_CELL<T>* cell,
        void (*map_function)(void* pointer,const OCTREE_CELL<T>* c1,const OCTREE_CELL<T>* c2,int axis));
    static void Map_Face_Process_Face(void* pointer,const OCTREE_CELL<T>* c1,const OCTREE_CELL<T>* c2,const int orient,
        void (*map_function)(void* pointer,const OCTREE_CELL<T>* c1,const OCTREE_CELL<T>* c2,int axis));
    template<class T3> static void Map_Face_Process_Cell_For_Reduction(void* pointer,const OCTREE_CELL<T>* cell,ARRAY<T3>& face_values,
        T3 (*map_function)(void* pointer,const T3 face_value_1,const T3 face_value_2,const T3 face_value_3,const T3 face_value_4));
    template<class T3> static T3 Map_Face_Process_Face_For_Reduction(void* pointer,const OCTREE_CELL<T>* c1,const OCTREE_CELL<T>* c2,const int orient,ARRAY<T3>& face_values,
        T3 (*map_function)(void* pointer,const T3 face_value_1,const T3 face_value_2,const T3 face_value_3,const T3 face_value_4));
    static void Map_Cell_Process_Cell(void* pointer,const OCTREE_CELL<T>* cell,
        void (*map_function)(void* pointer,const OCTREE_CELL<T>* c1));
//#####################################################################
};
}
#endif
#endif
