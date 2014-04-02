//#####################################################################
// Copyright 2009, Jon Gretarsson.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Namespace SIGNED_DISTANCE
//##################################################################### 
#include <PhysBAM_Tools/Grids_Uniform/GRID.h>
#include <PhysBAM_Tools/Grids_Uniform_Arrays/ARRAYS_ND.h>
#include <PhysBAM_Tools/Log/LOG.h>
#include <PhysBAM_Tools/Math_Tools/constants.h>
#include <PhysBAM_Geometry/Basic_Geometry/BOX.h>
#include <PhysBAM_Geometry/Grids_Uniform_Computations/SEGMENTED_CURVE_2D_SIGNED_DISTANCE.h>
#include <PhysBAM_Geometry/Topology_Based_Geometry/SEGMENTED_CURVE_2D.h>
namespace PhysBAM{

namespace SIGNED_DISTANCE{
//#####################################################################
// Function Calculate
//#####################################################################
template<class T> void Calculate(SEGMENTED_CURVE_2D<T>& curve,const GRID<VECTOR<T,2> >& grid,ARRAY<T,VECTOR<int,2> >& phi,bool print_progress)
{
    bool bounding_box_defined=curve.bounding_box!=0;if(!bounding_box_defined) curve.Update_Bounding_Box();

    phi.Resize(1,grid.counts.x,1,grid.counts.y);

    T epsilon=(T)1e-8*grid.min_dX;
    int total_cells=grid.counts.x*grid.counts.y,cells_done=0,progress=-1;
    for(int i=1;i<=grid.counts.x;i++) for(int j=1;j<=grid.counts.y;j++){
        phi(i,j)=curve.Calculate_Signed_Distance(grid.X(i,j),epsilon);
        if(print_progress){
            cells_done++;int new_progress=(int)((T)100*cells_done/total_cells);
            if(new_progress > progress){
#ifndef COMPILE_WITHOUT_READ_WRITE_SUPPORT
                LOG::cout<<new_progress<<"% "<<std::flush;
#endif
                progress=new_progress;}}}
#ifndef COMPILE_WITHOUT_READ_WRITE_SUPPORT
    if(print_progress) LOG::cout<<std::endl;
#endif

    if(!bounding_box_defined){delete curve.bounding_box;curve.bounding_box=0;}
}
//####################################################################
template void Calculate(SEGMENTED_CURVE_2D<float>&,const GRID<VECTOR<float,2> >&,ARRAY<float,VECTOR<int,2> >&,bool);
#ifndef COMPILE_WITHOUT_DOUBLE_SUPPORT
template void Calculate(SEGMENTED_CURVE_2D<double>&,const GRID<VECTOR<double,2> >&,ARRAY<double,VECTOR<int,2> >&,bool);
#endif
};
};
