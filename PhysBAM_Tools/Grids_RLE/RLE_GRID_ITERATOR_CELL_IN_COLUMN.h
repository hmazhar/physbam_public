//#####################################################################
// Copyright 2005, Eran Guendelman, Geoffrey Irving, Jerry Talton.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class RLE_GRID_ITERATOR_CELL_IN_COLUMN
//#####################################################################
#ifndef COMPILE_WITHOUT_RLE_SUPPORT
#ifndef __RLE_GRID_ITERATOR_CELL_IN_COLUMN__
#define __RLE_GRID_ITERATOR_CELL_IN_COLUMN__

#include <PhysBAM_Tools/Arrays/ARRAY.h>
#include <PhysBAM_Tools/Interpolation/LINEAR_INTERPOLATION.h>
#include <PhysBAM_Tools/Math_Tools/RANGE.h>
namespace PhysBAM{

template<class T_GRID>
class RLE_GRID_ITERATOR_CELL_IN_COLUMN
{
    typedef typename T_GRID::SCALAR T;
    typedef typename T_GRID::RUN T_RUN;typedef typename T_GRID::BOX_HORIZONTAL_INT T_BOX_HORIZONTAL_INT;
public:
    const T_GRID& grid;
    const bool bottom_sentinels_negative,bottom_sentinels_positive;
    const short top_sentinels;
    int j_end,length,j,dj;
    const ARRAY<T_RUN>* column;
    const T_RUN *run_end,*run;

protected:
    RLE_GRID_ITERATOR_CELL_IN_COLUMN(const T_GRID& grid_input,const bool top_sentinels_input)
        :grid(grid_input),bottom_sentinels_negative(false),bottom_sentinels_positive(false),top_sentinels(top_sentinels_input)
    {}

    RLE_GRID_ITERATOR_CELL_IN_COLUMN(const T_GRID& grid_input,const short bottom_sentinels,const short top_sentinels_input)
        :grid(grid_input),bottom_sentinels_negative(bottom_sentinels<0),bottom_sentinels_positive(bottom_sentinels>0),top_sentinels(top_sentinels_input)
    {}

    void Initialize_Column()
    {assert(-1<=top_sentinels && top_sentinels<=1);
    run_end=&(*column)(column->m-1+top_sentinels);run=&(*column)(1+bottom_sentinels_positive);
    if(run->is_long){j_end=run->jmin;length=(run+1)->jmin-run->jmin;j=run->jmin-bottom_sentinels_negative;dj=0;}
    else{j_end=(run+1)->jmin-1;length=1;j=run->jmin;dj=0;}}

    void Initialize_Column(const int j_input,const T_RUN* run_input) // for internal use in block iterator
    {j=j_input;column=0;run_end=0;run=run_input;j_end=0;
    if(run->is_long){j=run->jmin;length=(run+1)->jmin-run->jmin;dj=0;}
    else{length=1;dj=j-run->jmin;}}
public:

    static T_BOX_HORIZONTAL_INT Sentinels()
    {return T_BOX_HORIZONTAL_INT::Zero_Box();}

    operator bool() const
    {return run<=run_end;}

    RLE_GRID_ITERATOR_CELL_IN_COLUMN<T_GRID>& operator++()
    {assert(*this);dj++;
    if(++j>j_end){
        if(++run>run_end) return *this; // finished
        if(run->is_long){j_end=run->jmin;length=run<&(*column)(column->m)?(run+1)->jmin-run->jmin:0;} // start a new run
        else{j_end=(run+1)->jmin-1;length=1;}
        j=run->jmin;dj=0;}
    return *this;}

    int Cell() const
    {return run->cell+dj;}

    int First_Face_Index(const int axis) const
    {return run->faces[2*axis-2]+dj;}

    int Second_Face_Index(const int axis) const
    {return run->faces[2*axis-1]+dj;}

    int Face(const int face) const
    {assert(1<=face&&face<=T_GRID::number_of_faces_per_cell);return run->faces[face-1]+dj;}

    int Face(const int axis,const int side) const
    {assert(1<=axis&&axis<=T_GRID::dimension && 0<=side&&side<=1);return run->faces[2*axis-2+side]+dj;}

    int Face_X(const int side) const
    {assert(0<=side&&side<=1);return run->faces[side]+dj;}

    int Face_Z(const int side) const
    {assert(T_GRID::dimension==3 && 0<=side&&side<=1);return run->faces[4+side]+dj;}

    int Face_Y() const
    {return run->cell+dj;}

    int jmax() const
    {return j+length;}

    bool Long() const
    {return run->is_long;}

    bool Short() const
    {return !Long();}

    bool First_In_Column() const
    {return !dj && run==&(*column)(1);}

    bool Last_In_Column() const
    {return run==run_end;}

    int Run_Index() const
    {return run-&(*column)(1);}

    T y() const
    {assert(Short());return grid.uniform_grid.Axis_X_plus_half(j,2);}

    template<class TV> TV Cell_Value(const ARRAY<TV>& u,const int j_goal) const
    {if(Short()){assert(j_goal==j);return u(Cell());}else return Cell_Value_Long(u,j_goal);}

    template<class TV> TV Cell_Value_Long(const ARRAY<TV>& u,const int j_goal) const
    {assert(Long() && grid.long_run_cells==2 && j<=j_goal && j_goal<jmax());int c=run->cell;
    return LINEAR_INTERPOLATION<T,TV>::Linear(u(c),u(c+1),(T)(j_goal-j)/length);}

    static RANGE<VECTOR<int,1> > Indices_In_Slice(const T_GRID& grid,const int i)
    {return grid.Cells_In_Slice(i);}

//#####################################################################
};
}
#endif
#endif
