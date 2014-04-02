//#####################################################################
// Copyright 2002-2008, Zhaosheng Bao, Robert Bridson, Ronald Fedkiw, Eran Guendelman, Geoffrey Irving, Sergey Koltakov, Nipun Kwatra, Frank Losasso, Andrew Selle, Tamar Shinar, Eftychios Sifakis, Jerry Talton, Joseph Teran.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class GRID
//#####################################################################
// A number of functions (e.g. the Clamp functions) assume the grid indexing starts at 1.  This way we can use truncation rather than floor because floor is really slow.
//#####################################################################
#ifndef __GRID__
#define __GRID__

#include <PhysBAM_Tools/Arrays/ARRAYS_FORWARD.h>
#include <PhysBAM_Tools/Grids_Uniform/FACE_INDEX.h>
#include <PhysBAM_Tools/Grids_Uniform/UNIFORM_GRID_FORWARD.h>
#include <PhysBAM_Tools/Grids_Uniform_Arrays/ARRAYS_ND.h>
#include <PhysBAM_Tools/Math_Tools/pow.h>
#include <PhysBAM_Tools/Math_Tools/RANGE.h>
#include <PhysBAM_Tools/Vectors/VECTOR_3D.h>
namespace PhysBAM{

template<class T> class UNIFORM_GRID_ITERATOR;
template<class T> class UNIFORM_GRID_ITERATOR_NODE;
template<class T> class UNIFORM_GRID_ITERATOR_CELL;
template<class T> class UNIFORM_GRID_ITERATOR_FACE;
template<class T_GRID> class BLOCK_UNIFORM;

template<class TV>
class GRID
{
    typedef typename TV::SCALAR T;typedef VECTOR<int,TV::dimension> TV_INT;
public:
    enum REGION {WHOLE_REGION,GHOST_REGION,BOUNDARY_REGION,INTERIOR_REGION,BOUNDARY_INTERIOR_REGION}; // for iterators
    static const int dimension=TV::dimension;
    static const int number_of_cells_per_block=1<<dimension;
    static const int number_of_faces_per_block=dimension==1?3:dimension==2?12:36;
    static const int number_of_incident_faces_per_block=2*dimension-2;
    static const int number_of_nodes_per_face=1<<(dimension-(dimension!=0));
    static const int number_of_nodes_per_cell=1<<dimension;
    static const int number_of_cells_per_node=1<<dimension;
    static const int number_of_neighbors_per_node=2*dimension;
    static const int number_of_neighbors_per_cell=2*dimension;
    static const int number_of_faces_per_cell=2*dimension;
    static const int number_of_one_ring_neighbors_per_cell=dimension==1?2:dimension==2?8:26;
    typedef GRID<TV> UNIFORM_GRID;
    typedef UNIFORM_GRID_ITERATOR<TV> BASE_ITERATOR;
    typedef UNIFORM_GRID_ITERATOR_NODE<TV> NODE_ITERATOR;
    typedef UNIFORM_GRID_ITERATOR_CELL<TV> CELL_ITERATOR;
    typedef UNIFORM_GRID_ITERATOR_FACE<TV> FACE_ITERATOR;

    typedef T SCALAR;
    typedef TV VECTOR_T;
    typedef TV_INT VECTOR_INT;
    typedef TV_INT INDEX; typedef FACE_INDEX<TV::m> INDEX_FACE;
    typedef UNIFORM_TAG<TV> GRID_TAG;
    typedef BLOCK_UNIFORM<GRID<TV> > BLOCK;

    TV_INT counts;
    RANGE<TV> domain;
    TV dX;
    TV one_over_dX;
    TV_INT numbers_of_cells; // not saved to file
    T min_dX;
    T MAC_offset; // 0 for a regular grid and .5 for a MAC grid

    GRID()
    {
        Initialize(TV_INT(),RANGE<TV>::Unit_Box());
    }

    GRID(const int m_input,const int n_input,const int mn_input,const T xmin_input,const T xmax_input,const T ymin_input,const T ymax_input,const T zmin_input,const T zmax_input,
        const bool MAC_grid=false)
    {
        STATIC_ASSERT(dimension==3);
        Initialize(m_input,n_input,mn_input,xmin_input,xmax_input,ymin_input,ymax_input,zmin_input,zmax_input,MAC_grid);
    }
    
    GRID(const int m_input,const int n_input,const T xmin_input,const T xmax_input,const T ymin_input,const T ymax_input,const bool MAC_grid=false)
    {
        STATIC_ASSERT(dimension==2);
        Initialize(m_input,n_input,xmin_input,xmax_input,ymin_input,ymax_input,MAC_grid);
    }

    GRID(const int m_input,const T xmin_input,const T xmax_input,const bool MAC_grid=false)
    {
        STATIC_ASSERT(dimension==1);
        Initialize(m_input,xmin_input,xmax_input,MAC_grid);
    }

    GRID(const int m_input,const int n_input,const int mn_input,const RANGE<TV>& box,const bool MAC_grid=false)
    {
        STATIC_ASSERT(dimension==3);
        Initialize(m_input,n_input,mn_input,box,MAC_grid);
    }

    GRID(const int m_input,const int n_input,const RANGE<TV>& box,const bool MAC_grid=false)
    {
        STATIC_ASSERT(dimension==2);
        Initialize(m_input,n_input,box,MAC_grid);
    }

    GRID(const int m_input,const RANGE<TV>& box,const bool MAC_grid=false)
    {
        STATIC_ASSERT(dimension==1);
        Initialize(m_input,box,MAC_grid);
    }

    GRID(const T dx,const RANGE<TV>& box,const bool MAC_grid=false)
    {
        Initialize(dx,box,MAC_grid);
    }

    GRID(const TV_INT& counts,const RANGE<TV>& box,const bool MAC_grid=false)
    {
        Initialize(counts,box,MAC_grid);
    }

    GRID(const GRID<TV>& grid_input)
    {
        Initialize(grid_input.counts,grid_input.domain,grid_input.Is_MAC_Grid());
    }

    template<class T2>
    GRID(const GRID<T2>& grid_input)
    {
        Initialize(grid_input.counts,RANGE<TV>(grid_input.domain),grid_input.Is_MAC_Grid());
    }

    GRID<TV>& operator=(const GRID<TV>& grid_input)
    {if(this!=&grid_input) Initialize(grid_input.counts,grid_input.domain,grid_input.Is_MAC_Grid());return *this;}

    bool operator==(const GRID<TV>& grid) const
    {return grid.counts==counts && grid.domain==domain && grid.MAC_offset==MAC_offset;}

    inline bool operator!=(const GRID<TV>& grid) const
    {return !operator==(grid);}

    void Set_To_Double_Resolution_Grid(const GRID<TV>& grid)
    {Initialize(2*grid.numbers_of_cells+!grid.Is_MAC_Grid(),grid.domain,grid.Is_MAC_Grid());}

    void Initialize(const int m_input,const int n_input,const int mn_input,const T xmin_input,const T xmax_input,const T ymin_input,const T ymax_input,const T zmin_input,
        const T zmax_input,const bool MAC_grid=false)
    {STATIC_ASSERT(dimension==3);Initialize(TV_INT(m_input,n_input,mn_input),RANGE<TV>(xmin_input,xmax_input,ymin_input,ymax_input,zmin_input,zmax_input),MAC_grid);}

    void Initialize(const int m_input,const int n_input,const T xmin_input,const T xmax_input,const T ymin_input,const T ymax_input,const bool MAC_grid=false)
    {STATIC_ASSERT(dimension==2);Initialize(TV_INT(m_input,n_input),RANGE<TV>(xmin_input,xmax_input,ymin_input,ymax_input),MAC_grid);}

    void Initialize(const int m_input,const T xmin_input,const T xmax_input,const bool MAC_grid=false)
    {STATIC_ASSERT(dimension==1);Initialize(TV_INT(m_input),RANGE<TV>(xmin_input,xmax_input),MAC_grid);}

    void Initialize(const int m_input,const int n_input,const int mn_input,const RANGE<TV>& box,const bool MAC_grid=false)
    {STATIC_ASSERT(dimension==3);Initialize(TV_INT(m_input,n_input,mn_input),box,MAC_grid);}

    void Initialize(const int m_input,const int n_input,const RANGE<TV>& box,const bool MAC_grid=false)
    {STATIC_ASSERT(dimension==2);Initialize(TV_INT(m_input,n_input),box,MAC_grid);}

    void Initialize(const int m_input,const RANGE<TV>& box,const bool MAC_grid=false)
    {STATIC_ASSERT(dimension==1);Initialize(TV_INT(m_input),box,MAC_grid);}

    void Initialize(const T dx,const RANGE<TV>& box,const bool MAC_grid=false)
    {Initialize(TV_INT(box.Edge_Lengths()/dx),box,MAC_grid);}

    bool Is_MAC_Grid() const
    {return MAC_offset==.5;}

    bool Is_Isotropic(const T tolerance=1e-4) const
    {return (dX-dX.x).Max_Abs()<=tolerance;}

    const TV& DX() const
    {return dX;}

    const TV& One_Over_DX() const
    {return one_over_dX;}

    const TV_INT& Numbers_Of_Cells() const
    {return numbers_of_cells;}

    TV_INT Numbers_Of_Nodes() const
    {return numbers_of_cells+1;}

    T Minimum_Edge_Length() const
    {return min_dX;}

    T Maximum_Edge_Length() const
    {return dX.Max();}

    T One_Over_Cell_Size() const
    {return one_over_dX.Product();}

    T Cell_Size() const
    {return dX.Product();}

    TV Face_Sizes() const
    {return Cell_Size()*one_over_dX;}

    T Face_Size(const int axis) const
    {return Cell_Size()*one_over_dX(axis);}

    TV Xmin() const
    {return domain.min_corner;}

    TV Xmax() const
    {return domain.max_corner;}

    T Axis_X(const int i,const int axis) const // axis component of position at the i=1 to i=counts(axis) grid points
    {return domain.min_corner(axis)+(i-1+MAC_offset)*dX(axis);}

    T Axis_X_plus_half(const int i,const int axis) const
    {return domain.min_corner(axis)+(i-(T).5+MAC_offset)*dX(axis);}

    T Axis_X_minus_half(const int i,const int axis) const
    {return domain.min_corner(axis)+(i-(T)1.5+MAC_offset)*dX(axis);}

    TV X_plus_half(const TV_INT& index) const
    {return domain.min_corner+(TV(index)-((T).5+MAC_offset))*dX;}

    TV X_minus_half(const TV_INT& index) const
    {return domain.min_corner+(TV(index)-((T)1.5+MAC_offset))*dX;}

    TV X(const TV_INT& index) const
    {return domain.min_corner+(TV(index)+(MAC_offset-1))*dX;}

    TV X(const int i,const int j,const int ij) const
    {STATIC_ASSERT(dimension==3);return X(TV_INT(i,j,ij));}

    TV X(const int i,const int j) const
    {STATIC_ASSERT(dimension==2);return X(TV_INT(i,j));}

    TV X(const int i) const
    {STATIC_ASSERT(dimension==1);return X(TV_INT(i));}

    TV Node(const int i,const int j,const int ij) const
    {STATIC_ASSERT(dimension==3);return Node(TV_INT(i,j,ij));}

    TV Node(const int i,const int j) const
    {STATIC_ASSERT(dimension==2);return Node(TV_INT(i,j));}

    TV Node(const int i) const
    {STATIC_ASSERT(dimension==1);return Node(TV_INT(i));}

    TV Node(const TV_INT& index) const
    {return domain.min_corner+TV(index-1)*dX;}

    TV Center(const int i,const int j,const int ij) const
    {STATIC_ASSERT(dimension==3);return Center(TV_INT(i,j,ij));}

    TV Center(const int i,const int j) const
    {STATIC_ASSERT(dimension==2);return Center(TV_INT(i,j));}

    TV Center(const int i) const
    {STATIC_ASSERT(dimension==1);return Center(TV_INT(i));}

    TV Center(const TV_INT& index) const
    {return domain.min_corner+(TV(index)-(T).5)*dX;}

    TV Axis_X_Face(const INDEX_FACE& face) const
    {return Axis_X_Face(face.index,face.axis);}

    TV Axis_X_Face(const TV_INT& index,const int axis) const
    {TV adjusted=TV(index)-(T).5;adjusted(axis)-=(T).5;return domain.min_corner+adjusted*dX;}

    TV X_Face(const int i,const int j,const int ij) const
    {STATIC_ASSERT(dimension==3);return X_Face(TV_INT(i,j,ij));}

    TV X_Face(const int i,const int j) const
    {STATIC_ASSERT(dimension==2);return X_Face(TV_INT(i,j));}

    TV X_Face(const int i) const
    {STATIC_ASSERT(dimension==1);return X_Face(TV_INT(i));}

    TV X_Face(const TV_INT& index) const
    {return Axis_X_Face(index,1);}

    TV Y_Face(const int i,const int j,const int ij) const
    {STATIC_ASSERT(dimension==3);return Y_Face(TV_INT(i,j,ij));}

    TV Y_Face(const int i,const int j) const
    {STATIC_ASSERT(dimension==2);return Y_Face(TV_INT(i,j));}

    TV Y_Face(const TV_INT& index) const
    {return Axis_X_Face(index,2);}

    TV Z_Face(const int i,const int j,const int ij) const
    {STATIC_ASSERT(dimension==3);return Z_Face(TV_INT(i,j,ij));}

    TV Z_Face(const TV_INT& index) const
    {return Axis_X_Face(index,3);}

    TV Face(const int axis,const TV_INT& index) const
    {TV shifted_index(TV(index)-(T).5);shifted_index(axis)-=(T).5;return domain.min_corner+shifted_index*dX;}

    RANGE<TV> Face_Domain(const int axis,const TV_INT& index,const T thickness_over_two=0) const
    {TV dimensions=(T).5*dX;dimensions[axis]=thickness_over_two;RANGE<TV> domain(Face(axis,index));domain.Change_Size(dimensions);return domain;}

    static const VECTOR<TV_INT,8>& Binary_Counts(const VECTOR<int,3>&)
    {static const VECTOR<TV_INT,8> array(TV_INT(0,0,0),TV_INT(1,0,0),TV_INT(0,1,0),TV_INT(1,1,0),TV_INT(0,0,1),TV_INT(1,0,1),TV_INT(0,1,1),TV_INT(1,1,1));return array;}

    static const VECTOR<TV_INT,4>& Binary_Counts(const VECTOR<int,2>&)
    {static const VECTOR<TV_INT,4> array(TV_INT(0,0),TV_INT(1,0),TV_INT(0,1),TV_INT(1,1));return array;}

    static const VECTOR<TV_INT,2>& Binary_Counts(const VECTOR<int,1>&)
    {static const VECTOR<TV_INT,2> array(TV_INT(0),TV_INT(1));return array;}

    static const VECTOR<TV_INT,1>& Binary_Counts(const VECTOR<int,0>&)
    {static const VECTOR<TV_INT,1> array;return array;}

    static TV_INT Node_Cell_Index(const TV_INT& node,const int cell)
    {return node+Binary_Counts(node)(cell)-1;}

    static TV_INT Face_Node_Index(const int axis,const VECTOR<int,3>& face_index,const int node)
    {static const TV_INT corner_from_face_offset[3][4]={
        {TV_INT(0,0,0),TV_INT(0,1,0),TV_INT(0,0,1),TV_INT(0,1,1)},
        {TV_INT(0,0,0),TV_INT(1,0,0),TV_INT(0,0,1),TV_INT(1,0,1)},
        {TV_INT(0,0,0),TV_INT(1,0,0),TV_INT(0,1,0),TV_INT(1,1,0)}};
    assert(1<=axis&&axis<=3&&1<=node&&node<=4);return face_index+corner_from_face_offset[axis-1][node-1];}

    static TV_INT Face_Node_Index(const int axis,const VECTOR<int,2>& face_index,const int node)
    {assert(1<=node&&node<=2);TV_INT index=face_index;index[3-axis]+=node-1;return index;}

    static TV_INT Face_Node_Index(const int axis,const VECTOR<int,1>& face_index,const int node)
    {assert(axis==1&&node==1);return face_index;}

    static TV_INT Node_Face_Index(const int axis,const VECTOR<int,3>& node_index,const int face)
    {static const TV_INT face_from_node_offset[3][4]={
        {TV_INT(0,-1,-1),TV_INT(0,0,-1),TV_INT(0,-1,0),TV_INT(0,0,0)},
        {TV_INT(-1,0,-1),TV_INT(-1,0,0),TV_INT(0,0,-1),TV_INT(0,0,0)},
        {TV_INT(-1,-1,0),TV_INT(0,-1,0),TV_INT(-1,0,0),TV_INT(0,0,0)}};
    assert(1<=face&&face<=4&&1<=axis&&axis<=3);return node_index+face_from_node_offset[axis-1][face-1];}

    static TV_INT Node_Face_Index(const int axis,const VECTOR<int,2>& node_index,const int face)
    {assert(1<=face&&face<=2);TV_INT index=node_index;index[3-axis]+=face-2;return index;}

    static TV_INT Node_Face_Index(const int axis,const VECTOR<int,1>& node_index,const int face)
    {assert(axis==1&&face==1);return node_index;}

    void Face_Corner_To_Opposite_Corner_Vectors(const int axis,VECTOR<T,3> vectors[4])
    {static const TV multipliers[3][3]={
        {TV(0,-1,-1),TV(0,1,-1),TV(0,-1,1)},
        {TV(-1,0,-1),TV(1,0,-1),TV(-1,0,1)},
        {TV(-1,-1,0),TV(1,-1,0),TV(-1,1,0)}};
    assert(1<=axis&&axis<=3);vectors[3]=dX;vectors[3][axis]=0;for(int i=0;i<3;i++) vectors[i]=vectors[3]*multipliers[axis-1][i];}

    void Face_Corner_To_Opposite_Corner_Vectors(const int axis,VECTOR<T,2> vectors[2])
    {vectors[1]=dX;vectors[1][axis]=0;vectors[0]=(T)-1*vectors[1];}

    T Face_Corner_To_Opposite_Corner_Length(const int axis) const
    {TV face_dimensions=dX;face_dimensions[axis]=0;return face_dimensions.Magnitude();}

    TV_INT Index(const TV& location) const // returns the left, bottom and front indices on a regular grid
    {return TV_INT(floor((location-domain.min_corner)*one_over_dX+(1-MAC_offset)));} // note that "floor" is expensive

    void Cell(const TV& location,TV_INT& index,const int number_of_ghost_cells) const // returns the left, bottom and front
    {int number_of_ghost_cells_plus_one=number_of_ghost_cells+1; // Add before casting to avoid negatives rounding up to zero
    index=TV_INT((location-domain.min_corner)*one_over_dX+number_of_ghost_cells_plus_one)-number_of_ghost_cells;}

    TV_INT Cell(const TV& location,const int number_of_ghost_cells) const // returns the left, bottom and front
    {TV_INT index;Cell(location,index,number_of_ghost_cells);return index;}

    RANGE<TV> Cell_Domain(const TV_INT& index) const
    {TV corner=domain.min_corner+TV(index)*dX;return RANGE<TV>(corner-dX,corner);}

    static void Cells_Touching_Face(const int axis,const TV_INT& face_index,TV_INT& cell1,TV_INT& cell2)
    {cell2=face_index;cell1=face_index;cell1[axis]-=1;}

    TV Clamp(const TV& location) const
    {return domain.Clamp(location);}

    TV Clamp_Component(const TV& location,const int component) const
    {TV loc(location);loc(component)=clamp(loc(component),domain.min_corner(component),domain.max_corner(component));return loc;}

    TV Clamp(const TV& location,int number_of_ghost_cells) const // clamps to the grid (with ghost cells)
    {TV extra=(T)number_of_ghost_cells*dX;return clamp(location,domain.min_corner-extra,domain.max_corner+extra);}

    TV_INT Clamp_Min(const TV_INT& index) const
    {return clamp_min(index,TV_INT::All_Ones_Vector());}

    TV_INT Clamp_Max(const TV_INT& index) const
    {return clamp_max(index,counts);}

    void Clamp(TV_INT& index) const
    {index=clamp(index,TV_INT::All_Ones_Vector(),counts);}
    
    void Clamp(TV_INT& index,const int number_of_ghost_cells) const
    {index=clamp(index,TV_INT()-(number_of_ghost_cells-1),counts+number_of_ghost_cells);}

    TV_INT Clamped_Index(const TV& location) const
    {return clamp(1+TV_INT(((location-domain.min_corner)*one_over_dX-MAC_offset)),1+TV_INT(),counts);}
    
    TV_INT Clamped_Face_Index(const TV& location) const
    {return clamp(1+TV_INT(((location-domain.min_corner)*one_over_dX-MAC_offset)),1+TV_INT(),counts+1);}

    TV_INT Clamped_Index_End_Minus_One(const TV& location) const
    {return clamp(1+TV_INT(((location-domain.min_corner)*one_over_dX-MAC_offset)),1+TV_INT(),counts-1);}
        
    TV_INT Clamped_Index_End_Minus_One(const TV& location,const int number_of_ghost_cells) const
    {return clamp(1+TV_INT((location-domain.min_corner)*one_over_dX-MAC_offset+number_of_ghost_cells)-number_of_ghost_cells,TV_INT()-(number_of_ghost_cells-1),counts+(number_of_ghost_cells-1));}

    TV_INT Clamp_To_Cell(const TV& location) const
    {return clamp(1+TV_INT((location-domain.min_corner)*one_over_dX),TV_INT::All_Ones_Vector(),numbers_of_cells);}

    TV_INT Clamp_To_Cell(const TV& location,const int number_of_ghost_cells) const
    {return clamp(1+TV_INT((location-domain.min_corner)*one_over_dX+number_of_ghost_cells)-number_of_ghost_cells,TV_INT::All_Ones_Vector()-number_of_ghost_cells,numbers_of_cells+number_of_ghost_cells);}

    RANGE<TV_INT> Clamp_To_Cell(const RANGE<TV>& box,const int number_of_ghost_cells) const
    {return RANGE<TV_INT>(Clamp_To_Cell(box.Minimum_Corner(),number_of_ghost_cells),Clamp_To_Cell(box.Maximum_Corner(),number_of_ghost_cells));}

    TV_INT Block_Index(const TV& X,const int number_of_ghost_cells) const // index of node at center of block
    {assert(Is_MAC_Grid());return Clamped_Index_End_Minus_One(X,number_of_ghost_cells)+TV_INT::All_Ones_Vector();}

    TV_INT Closest_Node(const TV& location) const
    {return clamp(TV_INT((location-domain.min_corner)*one_over_dX+(T)1.5),TV_INT::All_Ones_Vector(),counts);}

    const RANGE<TV>& Domain() const
    {return domain;}

    RANGE<TV> Ghost_Domain(const int number_of_ghost_cells) const
    {TV expand(dX*number_of_ghost_cells);return RANGE<TV>(domain.min_corner-expand,domain.max_corner+expand);}

    const TV_INT& Counts() const
    {return counts;}

    RANGE<TV_INT> Domain_Indices(const int ghost_cells=0) const
    {return RANGE<TV_INT>(TV_INT::All_Ones_Vector()-ghost_cells,counts+ghost_cells);}

    RANGE<TV_INT> Node_Indices(const int ghost_cells=0) const
    {return RANGE<TV_INT>(TV_INT()+1,Numbers_Of_Nodes()).Thickened(ghost_cells);}

    RANGE<TV_INT> Cell_Indices(const int ghost_cells=0) const
    {return RANGE<TV_INT>(TV_INT()+1,numbers_of_cells).Thickened(ghost_cells);}

    RANGE<TV_INT> Block_Indices(const int ghost_cells=0) const
    {return Node_Indices(ghost_cells);}

    bool Inside_Domain(const TV_INT cell_index,const int ghost_cells=0) const
    {return Domain_Indices(ghost_cells).Lazy_Inside(cell_index);}

    VECTOR<RANGE<TV_INT>,TV::dimension> Face_Indices(const int ghost_cells=0) const
    {VECTOR<RANGE<TV_INT>,TV::dimension> v;for(int i=1;i<=TV::dimension;i++) v(i)=Get_Axis_X_Face_Grid(i).Node_Indices(ghost_cells);return v;}

    bool Outside(const TV& location) const
    {return domain.Lazy_Outside(location);}

    GRID<TV> Get_MAC_Grid() const
    {return GRID<TV>(numbers_of_cells,domain,true);}

    GRID<TV> Get_Regular_Grid() const
    {return GRID<TV>(Numbers_Of_Nodes(),domain,false);}

    GRID<TV> Get_Axis_X_Face_Grid(const int axis) const
    {TV_INT numbers=numbers_of_cells;numbers(axis)++;TV offset((T).5*dX);offset(axis)=0;
    return GRID<TV>(TV_INT::Componentwise_Max(TV_INT(),numbers),RANGE<TV>(domain.min_corner+offset,domain.max_corner-offset));}

    GRID<TV> Get_X_Face_Grid() const
    {return Get_Axis_X_Face_Grid(1);}

    GRID<TV> Get_Y_Face_Grid() const
    {return Get_Axis_X_Face_Grid(2);}

    GRID<TV> Get_Z_Face_Grid() const
    {return Get_Axis_X_Face_Grid(3);}

    GRID<TV> Get_Face_Grid(const int axis) const
    {switch(axis){case 1:return Get_X_Face_Grid();case 2:return Get_Y_Face_Grid();default:assert(axis==3);return Get_Z_Face_Grid();}}

    GRID<TV> Get_Regular_Grid_At_MAC_Positions() const
    {assert(Is_MAC_Grid());TV expansion=(T).5*dX;return GRID<TV>(counts,RANGE<TV>(domain.min_corner+expansion,domain.max_corner-expansion));}

    GRID<TV> Get_MAC_Grid_At_Regular_Positions() const
    {assert(!Is_MAC_Grid());TV expansion=(T).5*dX;return GRID<TV>(counts,RANGE<TV>(domain.min_corner-expansion,domain.max_corner+expansion),true);}

    GRID<VECTOR<T,TV::dimension-1> > Remove_Dimension(int dimension) const
    {return GRID<VECTOR<T,TV::dimension-1> >(counts.Remove_Index(dimension),domain.Remove_Dimension(dimension),Is_MAC_Grid());}

    GRID<VECTOR<T,TV::dimension-1> > Get_Horizontal_Grid() const
    {return Remove_Dimension(2);}

    GRID<VECTOR<T,1> > Get_1D_Grid(const int axis) const
    {return GRID<VECTOR<T,1> >(counts[axis],domain.Minimum_Corner()[axis],domain.Maximum_Corner()[axis],Is_MAC_Grid());}

    static VECTOR<int,3> Node_Neighbor(const VECTOR<int,3>& index,const int i)
    {static const VECTOR<int,3> neighbor_offset[6]={VECTOR<int,3>(-1,0,0),VECTOR<int,3>(1,0,0),VECTOR<int,3>(0,-1,0),VECTOR<int,3>(0,1,0),VECTOR<int,3>(0,0,-1),VECTOR<int,3>(0,0,1)};
    assert(1<=i&&i<=6);return index+neighbor_offset[i-1];}

    static VECTOR<int,2> Node_Neighbor(const VECTOR<int,2>& index,const int i)
    {static const VECTOR<int,2> neighbor_offset[4]={VECTOR<int,2>(-1,0),VECTOR<int,2>(1,0),VECTOR<int,2>(0,-1),VECTOR<int,2>(0,1)};
    assert(1<=i&&i<=4);return index+neighbor_offset[i-1];}

    static VECTOR<int,1> Node_Neighbor(const VECTOR<int,1>& index,const int i) // i=1 to 2
    {static const VECTOR<int,1> neighbor_offset[2]={VECTOR<int,1>(-1),VECTOR<int,1>(1)};
    assert(1<=i&&i<=2);return index+neighbor_offset[i-1];}

    static VECTOR<int,0> Node_Neighbor(const VECTOR<int,0>& index,const int i) // i=1
    {assert(1==i);return index;}

    void Cells_Neighboring_Node(const VECTOR<int,3>& node,VECTOR<int,3> cells[8]) const
    {cells[0]=VECTOR<int,3>(node.x-1,node.y-1,node.z-1);cells[1]=VECTOR<int,3>(node.x,node.y-1,node.z-1);cells[2]=VECTOR<int,3>(node.x-1,node.y,node.z-1);
    cells[3]=VECTOR<int,3>(node.x,node.y,node.z-1);cells[4]=VECTOR<int,3>(node.x-1,node.y-1,node.z);cells[5]=VECTOR<int,3>(node.x,node.y-1,node.z);
    cells[6]=VECTOR<int,3>(node.x-1,node.y,node.z);cells[7]=VECTOR<int,3>(node.x,node.y,node.z);}

    void Cells_Neighboring_Node(const VECTOR<int,2>& node,VECTOR<int,2> cells[4]) const
    {cells[0]=TV_INT(node.x-1,node.y-1);cells[1]=TV_INT(node.x,node.y-1);cells[2]=TV_INT(node.x-1,node.y);cells[3]=TV_INT(node.x,node.y);}

    void Cells_Neighboring_Node(const VECTOR<int,1>& node,VECTOR<int,1> cells[2]) const
    {cells[0]=TV_INT(node.x-1);cells[1]=TV_INT(node.x);}

    void Nodes_In_Cell_From_Minimum_Corner_Node(const VECTOR<int,3>& minimum_corner_node,VECTOR<int,3> nodes[8]) const
    {int i=minimum_corner_node.x,j=minimum_corner_node.y,ij=minimum_corner_node.z;
    nodes[0]=minimum_corner_node;nodes[1]=VECTOR<int,3>(i+1,j,ij);nodes[2]=VECTOR<int,3>(i,j+1,ij);nodes[3]=VECTOR<int,3>(i+1,j+1,ij);
    nodes[4]=VECTOR<int,3>(i,j,ij+1);nodes[5]=VECTOR<int,3>(i+1,j,ij+1);nodes[6]=VECTOR<int,3>(i,j+1,ij+1);nodes[7]=VECTOR<int,3>(i+1,j+1,ij+1);}

    void Nodes_In_Cell_From_Minimum_Corner_Node(const VECTOR<int,2>& minimum_corner_node,VECTOR<int,2> nodes[4]) const
    {nodes[0]=minimum_corner_node;nodes[1]=VECTOR<int,2>(minimum_corner_node.x+1,minimum_corner_node.y);nodes[2]=VECTOR<int,2>(minimum_corner_node.x,minimum_corner_node.y+1);
    nodes[3]=VECTOR<int,2>(minimum_corner_node.x+1,minimum_corner_node.y+1);}

    void Nodes_In_Cell_From_Minimum_Corner_Node(const VECTOR<int,1>& minimum_corner_node,VECTOR<int,1> nodes[2]) const
    {nodes[0]=minimum_corner_node;nodes[1]=VECTOR<int,1>(minimum_corner_node.x+1);}

    void Node_Locations_In_Cell_From_Minimum_Corner_Node(const TV_INT& minimum_corner_node,ARRAY<TV>& node_locations) const
    {assert(node_locations.m==number_of_nodes_per_cell);
    TV_INT nodes[number_of_nodes_per_cell];Nodes_In_Cell_From_Minimum_Corner_Node(minimum_corner_node,nodes);
    for(int i=1;i<=number_of_nodes_per_cell;i++) node_locations(i)=Node(nodes[i-1]);}

    static TV_INT First_Face_Index_In_Cell(const int axis,const TV_INT& cell_index)
    {return cell_index;}

    static TV_INT Second_Face_Index_In_Cell(const int axis,const TV_INT& cell_index)
    {TV_INT face_index=cell_index;face_index[axis]++;return face_index;}

    static TV_INT One_Ring_Neighbor(const VECTOR<int,3>& index,const int i)
    {static const TV_INT neighbor_offset[26]={TV_INT(-1,-1,-1),TV_INT(0,-1,-1),TV_INT(1,-1,-1),TV_INT(-1,0,-1),TV_INT(0,0,-1),
    TV_INT(1,0,-1),TV_INT(-1,1,-1),TV_INT(0,1,-1),TV_INT(1,1,-1),TV_INT(-1,-1,0),TV_INT(0,-1,0),TV_INT(1,-1,0),
    TV_INT(-1,0,0),TV_INT(1,0,0),TV_INT(-1,1,0),TV_INT(0,1,0),TV_INT(1,1,0),TV_INT(-1,-1,1),TV_INT(0,-1,1),
    TV_INT(1,-1,1),TV_INT(-1,0,1),TV_INT(0,0,1),TV_INT(1,0,1),TV_INT(-1,1,1),TV_INT(0,1,1), TV_INT(1,1,1)};
    assert(1<=i&&i<=26);return index+neighbor_offset[i-1];}

    static TV_INT One_Ring_Neighbor(const VECTOR<int,2>& index,const int i)
    {static const TV_INT neighbor_offset[8]={TV_INT(-1,-1),TV_INT(0,-1),TV_INT(1,-1),TV_INT(-1,0),TV_INT(1,0),TV_INT(-1,1),
        TV_INT(0,1),TV_INT(1,1)};
    assert(1<=i&&i<=8);return index+neighbor_offset[i-1];}

    static TV_INT One_Ring_Neighbor(const VECTOR<int,1>& index,const int i)
    {return Node_Neighbor(index,i);}

    static void Neighboring_Faces(VECTOR<INDEX_FACE,number_of_faces_per_cell>& n,const TV_INT& index)
    {for(int a=1;a<=TV::m;a++){INDEX_FACE fi(a,index);n(a*2-1)=fi;fi.index(a)++;n(a*2)=fi;}}

    template<class T2> void Put_Ghost(const T2& constant,ARRAYS_ND_BASE<VECTOR<T2,3> >& array,const int ghost_cells) const
    {for(int j=1-ghost_cells;j<=counts.y+ghost_cells;j++) for(int ij=1-ghost_cells;ij<=counts.z+ghost_cells;ij++) for(int s=1;s<=ghost_cells;s++) array(1-s,j,ij)=array(counts.x+s,j,ij)=constant;
    for(int i=1;i<=counts.x;i++) for(int ij=1-ghost_cells;ij<=counts.z+ghost_cells;ij++) for(int s=1;s<=ghost_cells;s++) array(i,1-s,ij)=array(i,counts.y+s,ij)=constant;
    for(int i=1;i<=counts.x;i++) for(int j=1;j<=counts.y;j++) for(int s=1;s<=ghost_cells;s++) array(i,j,1-s)=array(i,j,counts.z+s)=constant;}

    template<class T2> void Put_Ghost(const T2& constant,ARRAYS_ND_BASE<VECTOR<T2,2> >& array,const int ghost_cells) const
    {for(int j=1-ghost_cells;j<=counts.y+ghost_cells;j++) for(int s=1;s<=ghost_cells;s++) array(1-s,j)=array(counts.x+s,j)=constant;
    for(int i=1;i<=counts.x;i++) for(int s=1;s<=ghost_cells;s++) array(i,1-s)=array(i,counts.y+s)=constant;}

    template<class T2> void Put_Ghost(const T2& constant,ARRAYS_ND_BASE<VECTOR<T2,1> >& array,const int ghost_cells) const
    {for(int s=1;s<=ghost_cells;s++) array(1-s)=array(counts.x+s)=constant;}
    
//#####################################################################
    void Initialize(const TV_INT& counts_input,const RANGE<TV>& box,const bool MAC_grid=false);
    static GRID<TV> Create_Grid_Given_Cell_Size(const RANGE<TV>& domain,const T cell_size,const bool mac_grid,const int boundary_nodes=0);
    static GRID<TV> Create_Even_Sized_Grid_Given_Cell_Size(const RANGE<TV>& domain,const T cell_size,const bool mac_grid,const int boundary_nodes=0);
//#####################################################################
};
}
#endif
