//#####################################################################
// Copyright 2005-2010, Geoffrey Irving, Michael Lentine, Frank Losasso, Andrew Selle, Tamar Shinar.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#ifndef __INTERPOLATION_UNIFORM__
#define __INTERPOLATION_UNIFORM__

#include <PhysBAM_Tools/Grids_Uniform/BLOCK_UNIFORM.h>
#include <PhysBAM_Tools/Grids_Uniform_Arrays/GRID_ARRAYS_POLICY_UNIFORM.h>
#include <PhysBAM_Tools/Grids_Uniform_Interpolation/FACE_LOOKUP_UNIFORM.h>
#include <PhysBAM_Tools/Grids_Uniform_Interpolation/INTERPOLATION_UNIFORM_FORWARD.h>
#include <PhysBAM_Tools/Log/DEBUG_UTILITIES.h>
#include <PhysBAM_Tools/Utilities/NONCOPYABLE.h>
#include <PhysBAM_Tools/Utilities/PHYSBAM_OVERRIDE.h>
namespace PhysBAM{

template<class T_GRID,class T2,class T_FACE_LOOKUP> // T_FACE_LOOKUP=FACE_LOOKUP_UNIFORM<T_GRID>
class INTERPOLATION_UNIFORM:public NONCOPYABLE
{
    typedef typename T_GRID::VECTOR_T TV;typedef typename T_GRID::SCALAR T;
    STATIC_ASSERT((IS_SAME<typename T_GRID::GRID_TAG,UNIFORM_TAG<TV> >::value));
    typedef typename GRID_ARRAYS_POLICY<T_GRID>::ARRAYS_SCALAR T_ARRAYS_SCALAR;
    //typedef typename GRID_ARRAYS_POLICY<T_GRID>::ARRAYS_SCALAR::template REBIND<T2>::TYPE T_ARRAYS_T2;
    typedef typename T_GRID::VECTOR_INT TV_INT;
    typedef ARRAYS_ND_BASE<VECTOR<T2,TV::dimension> > T_ARRAYS_T2;
    typedef typename T_GRID::CELL_ITERATOR CELL_ITERATOR;typedef typename T_GRID::FACE_ITERATOR FACE_ITERATOR;
    typedef typename GRID_ARRAYS_POLICY<T_GRID>::FACE_ARRAYS T_FACE_ARRAYS_SCALAR;
public:
    template<class T3> struct REBIND{typedef INTERPOLATION_UNIFORM<T_GRID,T3,T_FACE_LOOKUP> TYPE;};

    virtual ~INTERPOLATION_UNIFORM()
    {}

    static TV_INT Clamped_Index(const T_GRID& grid,const T_ARRAYS_T2& u,const TV& location,int start_offset=0,int end_offset=0) 
    {TV_INT index;RANGE<TV_INT> u_domain=u.Domain_Indices();RANGE<TV> grid_domain=grid.Domain(); 
    TV_INT M_start=u_domain.Minimum_Corner(),M_end=u_domain.Maximum_Corner();TV min_X=grid_domain.Minimum_Corner();TV one_over_DX=grid.One_Over_DX(); 
    for(int axis=1;axis<=T_GRID::dimension;axis++){ 
        index[axis]=min(M_end[axis]-end_offset,M_start[axis]+start_offset+max(0,(int)((location[axis]-min_X[axis])*one_over_DX[axis]-M_start[axis]-start_offset+1-grid.MAC_offset)));} 
    return index;}

    static TV_INT Clamped_Index_End_Minus_One(const T_GRID& grid,const T_ARRAYS_T2& u,const TV& location)
    {return Clamped_Index(grid,u,location,0,1);}

    static TV_INT Clamped_Index_Interior(const T_GRID& grid,const T_ARRAYS_T2& u,const TV& location)
    {return Clamped_Index(grid,u,location,1,1);}

    static TV_INT Clamped_Index_Interior_End_Minus_One(const T_GRID& grid,const T_ARRAYS_T2& u,const TV& location)
    {return Clamped_Index(grid,u,location,1,2);}

    void Populate_New_Array(const T_GRID& grid,const T_ARRAYS_T2& u,const T_GRID& grid_new,T_ARRAYS_T2& u_new)
    {for(CELL_ITERATOR iterator(grid_new,u_new.Domain_Indices());iterator.Valid();iterator.Next()){ // CELL ITERATOR works for nodal
        u_new(iterator.Cell_Index())=Clamped_To_Array(grid,u,grid_new.X(iterator.Cell_Index()));}}

    void Populate_New_Array(const T_GRID& grid,const T_FACE_ARRAYS_SCALAR& u,const T_GRID& grid_new,T_FACE_ARRAYS_SCALAR& u_new)
    {FACE_LOOKUP_UNIFORM<T_GRID> lookup(u);
    for(FACE_ITERATOR iterator(grid_new);iterator.Valid();iterator.Next()){int axis=iterator.Axis();
        u_new.Component(axis)(iterator.Face_Index())=Clamped_To_Array_Face_Component(axis,grid,lookup,iterator.Location());}}

    T2 Clamped_To_Array_Cell(const T_GRID& grid,const T_ARRAYS_T2& u,const TV& X) const
    {return Clamped_To_Array(grid,u,X);}

    T2 Clamped_To_Array_Node(const T_GRID& grid,const T_ARRAYS_T2& u,const TV& X) const
    {return Clamped_To_Array(grid,u,X);}

    TV Clamped_To_Array_Face(const T_GRID& grid,const typename T_FACE_LOOKUP::LOOKUP& u,const TV& X) const
    {assert(grid.Is_MAC_Grid());return From_Block_Face(grid,BLOCK_UNIFORM<T_GRID>(grid,X,u.Number_Of_Ghost_Cells()),u,X);}

    T Clamped_To_Array_Face_Component(const int axis,const T_GRID& grid,const typename T_FACE_LOOKUP::LOOKUP& u,const TV& X) const
    {assert(grid.Is_MAC_Grid());return From_Block_Face_Component(axis,grid,BLOCK_UNIFORM<T_GRID>(grid,X,u.Number_Of_Ghost_Cells()),u,X);}
    
    ARRAY<PAIR<FACE_INDEX<TV::dimension>,T> > Clamped_To_Array_Face_Component_Weights(const int axis,const T_GRID& grid,const typename T_FACE_LOOKUP::LOOKUP& u,const TV& X) const
    {assert(grid.Is_MAC_Grid());return From_Block_Face_Component_Weights(axis,grid,BLOCK_UNIFORM<T_GRID>(grid,X,u.Number_Of_Ghost_Cells()),u,X);}

    VECTOR<T2,2> Extrema_Clamped_To_Array_Cell(const T_GRID& grid,const T_ARRAYS_T2& u_min,const T_ARRAYS_T2& u_max,const TV& X) const
    {return Extrema_Clamped_To_Array(grid,u_min,u_max,X);}

    VECTOR<T2,2> Extrema_Clamped_To_Array_Node(const T_GRID& grid,const T_ARRAYS_T2& u_min,const T_ARRAYS_T2& u_max,const TV& X) const
    {return Extrema_Clamped_To_Array(grid,u_min,u_max,X);}

    VECTOR<TV,2> Extrema_Clamped_To_Array_Face(const T_GRID& grid,const typename T_FACE_LOOKUP::LOOKUP& u_min,
        const typename T_FACE_LOOKUP::LOOKUP& u_max,const TV& X) const
    {assert(grid.Is_MAC_Grid());return Extrema_From_Block_Face(grid,BLOCK_UNIFORM<T_GRID>(grid,X,u_min.Number_Of_Ghost_Cells()),u_min,u_max,X);}

    VECTOR<T,2> Extrema_Clamped_To_Array_Face_Component(const int axis,const T_GRID& grid,const typename T_FACE_LOOKUP::LOOKUP& u_min,
        const typename T_FACE_LOOKUP::LOOKUP& u_max,const TV& X) const
    {assert(grid.Is_MAC_Grid());return Extrema_From_Block_Face_Component(axis,grid,BLOCK_UNIFORM<T_GRID>(grid,X,u_min.Number_Of_Ghost_Cells()),u_min,u_max,X);}

//#####################################################################
    virtual T2 Clamped_To_Array_No_Extrema(const T_GRID& grid,const T_ARRAYS_T2& u,const TV& X) const {PHYSBAM_FUNCTION_IS_NOT_DEFINED();}
    virtual ARRAY<PAIR<TV_INT,T> > Clamped_To_Array_Weights(const T_GRID& grid,const T_ARRAYS_T2& u,const TV& X) const {PHYSBAM_FUNCTION_IS_NOT_DEFINED();}
    virtual T2 Clamped_To_Array(const T_GRID& grid,const T_ARRAYS_T2& u,const TV& X) const {PHYSBAM_FUNCTION_IS_NOT_DEFINED();}
    virtual T2 From_Base_Node(const T_GRID& grid,const T_ARRAYS_T2& u,const TV& X,const TV_INT& index) const {PHYSBAM_FUNCTION_IS_NOT_DEFINED();}
    virtual TV From_Block_Face(const T_GRID& grid,const BLOCK_UNIFORM<T_GRID>& block,const typename T_FACE_LOOKUP::LOOKUP& u,const TV& X) const {PHYSBAM_FUNCTION_IS_NOT_DEFINED();}
    virtual T From_Block_Face_Component(const int axis,const T_GRID& grid,const BLOCK_UNIFORM<T_GRID>& block,const typename T_FACE_LOOKUP::LOOKUP& u,const TV& X) const
    {PHYSBAM_FUNCTION_IS_NOT_DEFINED();}
    virtual ARRAY<PAIR<FACE_INDEX<TV::dimension>,T> > From_Block_Face_Component_Weights(const int axis,const T_GRID& grid,const BLOCK_UNIFORM<T_GRID>& block,const typename T_FACE_LOOKUP::LOOKUP& u,const TV& X) const
    {PHYSBAM_FUNCTION_IS_NOT_DEFINED();}
    virtual VECTOR<T2,2> Extrema_Clamped_To_Array(const T_GRID& grid,const T_ARRAYS_T2& u_min,const T_ARRAYS_T2& u_max,const TV& X) const 
    {PHYSBAM_FUNCTION_IS_NOT_DEFINED();}
    virtual VECTOR<T2,2> Extrema_From_Base_Node(const T_GRID& grid,const T_ARRAYS_T2& u_min,const T_ARRAYS_T2& u_max,const TV& X,const TV_INT& index) const
    {PHYSBAM_FUNCTION_IS_NOT_DEFINED();}
    virtual VECTOR<TV,2> Extrema_From_Block_Face(const T_GRID& grid,const BLOCK_UNIFORM<T_GRID>& block,const typename T_FACE_LOOKUP::LOOKUP& u_min,
        const typename T_FACE_LOOKUP::LOOKUP& u_max,const TV& X) const
    {PHYSBAM_FUNCTION_IS_NOT_DEFINED();}
    virtual VECTOR<T,2> Extrema_From_Block_Face_Component(const int axis,const T_GRID& grid,const BLOCK_UNIFORM<T_GRID>& block,
        const typename T_FACE_LOOKUP::LOOKUP& u_min,const typename T_FACE_LOOKUP::LOOKUP& u_max,const TV& X) const
    {PHYSBAM_FUNCTION_IS_NOT_DEFINED();}
//#####################################################################
};
}
#endif
