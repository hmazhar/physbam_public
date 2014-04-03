//#####################################################################
// Copyright 2004-2005, Andrew Selle.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#include <PhysBAM_Tools/Arrays/ARRAY.h>
#include <PhysBAM_Tools/Math_Tools/RANGE.h>
#include <PhysBAM_Rendering/PhysBAM_OpenGL/OpenGL/OPENGL_TWO_PHASE_VELOCITY_MAGNITUDE_2D.h>
using namespace PhysBAM;
using namespace std;
//#####################################################################
// Function OPENGL_TWO_PHASE_VELOCITY_MAGNITUDE_2D
//#####################################################################
template<class T> OPENGL_TWO_PHASE_VELOCITY_MAGNITUDE_2D<T>::
OPENGL_TWO_PHASE_VELOCITY_MAGNITUDE_2D(GRID<TV>& grid,ARRAY<VECTOR<T,2> ,VECTOR<int,2> >& V_minus,ARRAY<VECTOR<T,2> ,VECTOR<int,2> >& V_plus,LEVELSET_2D<GRID<TV> >& levelset)
    :height_scale(0),grid(grid),V_minus(V_minus),V_plus(V_plus),levelset(levelset),minus(*(new ARRAY<VECTOR<T,3> >),*(new ARRAY<VECTOR<T,3> >)),
      plus(*(new ARRAY<VECTOR<T,3> >),*(new ARRAY<VECTOR<T,3> >))
{
}
//#####################################################################
// Function ~OPENGL_TWO_PHASE_VELOCITY_MAGNITUDE_2D
//#####################################################################
template<class T> OPENGL_TWO_PHASE_VELOCITY_MAGNITUDE_2D<T>::
~OPENGL_TWO_PHASE_VELOCITY_MAGNITUDE_2D()
{
    delete &(minus.vector_field);
    delete &(minus.vector_locations);
    delete &(plus.vector_field);
    delete &(plus.vector_locations);
}
//#####################################################################
// Function Display
//#####################################################################
template<class T> void OPENGL_TWO_PHASE_VELOCITY_MAGNITUDE_2D<T>::
Display(const int in_color) const
{
    minus.Display(in_color);plus.Display(in_color);
}
//#####################################################################
// Function Update
//#####################################################################
template<class T> void OPENGL_TWO_PHASE_VELOCITY_MAGNITUDE_2D<T>::
Update()
{
    minus.vector_field.Resize(0);minus.vector_locations.Resize(0);
    plus.vector_field.Resize(0);plus.vector_locations.Resize(0);

    int minus_count=0,plus_count=0;

    for(int i=1;i<=grid.counts.x;i++)for(int j=1;j<=grid.counts.y;j++) if(levelset.phi(i,j)<=0)minus_count++;else plus_count++;

    int minus_index=1,plus_index=1;
    minus.vector_field.Resize(minus_count);minus.vector_locations.Resize(minus_count);
    plus.vector_field.Resize(plus_count);plus.vector_locations.Resize(plus_count);
    
    levelset.Compute_Normals();
    for (int i=V_minus.domain.min_corner.x;i<=V_minus.domain.max_corner.x;i++) for(int j=V_minus.domain.min_corner.y;j<=V_minus.domain.max_corner.y;j++){
            VECTOR<T,2> velocity_2d;
            if(levelset.phi(i,j)<=0){
                velocity_2d=V_minus(i,j);
                minus.vector_field(minus_index)=VECTOR<T,3>(velocity_2d.x,velocity_2d.y,0);
                VECTOR<T,2> X=grid.X(VECTOR<int,2>(i,j));
                minus.vector_locations(minus_index)=VECTOR<T,3>(X.x,X.y,-height_scale*abs(VECTOR<T,2>::Dot_Product((*levelset.normals)(i,j),velocity_2d)));
                minus_index++;}
            else{
                velocity_2d=V_plus(i,j);
                plus.vector_field(plus_index)=VECTOR<T,3>(velocity_2d.x,velocity_2d.y,0);
                VECTOR<T,2> X=grid.X(VECTOR<int,2>(i,j));
                plus.vector_locations(plus_index)=VECTOR<T,3>(X.x,X.y,-height_scale*abs(VECTOR<T,2>::Dot_Product((*levelset.normals)(i,j),velocity_2d)));
                plus_index++;}}
}
//#####################################################################
// Function Bounding_Box
//#####################################################################
template<class T> RANGE<VECTOR<float,3> > OPENGL_TWO_PHASE_VELOCITY_MAGNITUDE_2D<T>::
Bounding_Box() const
{
    return RANGE<VECTOR<float,3> >(grid.domain.min_corner.x,grid.domain.max_corner.x,grid.domain.min_corner.y,grid.domain.max_corner.y,0,0);
}
//#####################################################################
// Function Scale_Vector_Size
//#####################################################################
template<class T> void OPENGL_TWO_PHASE_VELOCITY_MAGNITUDE_2D<T>::
Scale_Vector_Size(const T scale)
{
    plus.size*=scale;minus.size*=scale;
}
//#####################################################################
// Function Scale_Height
//#####################################################################
template<class T> void OPENGL_TWO_PHASE_VELOCITY_MAGNITUDE_2D<T>::
Scale_Height(const T scale)
{
    height_scale=scale;
}
//#####################################################################
template class OPENGL_TWO_PHASE_VELOCITY_MAGNITUDE_2D<float>;
#ifndef COMPILE_WITHOUT_DOUBLE_SUPPORT
template class OPENGL_TWO_PHASE_VELOCITY_MAGNITUDE_2D<double>;
#endif
