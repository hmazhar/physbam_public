//#####################################################################
// Copyright 2005, Geoffrey Irving.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class RLE_RUN
//##################################################################### 
#ifndef COMPILE_WITHOUT_RLE_SUPPORT
#ifndef __RLE_RUN__
#define __RLE_RUN__

#include <PhysBAM_Tools/Arrays/ARRAYS_FORWARD.h>
namespace PhysBAM{

class RLE_RUN
{
public:
    bool is_long;
    short jmin;

    RLE_RUN()
        :is_long(false),jmin(0)
    {}

    RLE_RUN(const bool is_long_input,const short jmin_input)
        :is_long(is_long_input),jmin(jmin_input)
    {}

    RLE_RUN(const RLE_RUN& run_input)
        :is_long(run_input.is_long),jmin(run_input.jmin)
    {}

    bool operator==(const RLE_RUN& run) const
    {return is_long==run.is_long && jmin==run.jmin;}

    bool operator!=(const RLE_RUN& run) const
    {return !(*this==run);}

//#####################################################################
    template<class T_ARRAY_RUN> static void Column_Union(const T_ARRAY_RUN& column1,const T_ARRAY_RUN& column2,ARRAY<RLE_RUN>& column_union,const int minimum_long_run_length);
    static void Dilate_Vertical(const ARRAY<RLE_RUN>& column,ARRAY<RLE_RUN>& dilated_column,const int extra_cells);
    static void Dilate_Horizontal(const ARRAY<RLE_RUN>* columns,ARRAY<RLE_RUN>* new_columns,const int m_start,const int m_end,const int stride,const int extra_cells,
        const int minimum_long_run_length);
    template<class T_RUN> static void Split_At_Ground(const ARRAY<RLE_RUN>& column,ARRAY<T_RUN>& new_column,const int ground_j,const int minimum_long_run_length);
//#####################################################################
};
}
#endif
#endif
