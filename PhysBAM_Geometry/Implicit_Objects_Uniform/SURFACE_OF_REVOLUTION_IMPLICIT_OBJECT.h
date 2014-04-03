//#####################################################################
// Copyright 2007, Geoffrey Irving, Jerry Talton.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#ifndef COMPILE_WITHOUT_READ_WRITE_SUPPORT
#ifndef __SURFACE_OF_REVOLUTION_IMPLICIT_OBJECT__
#define __SURFACE_OF_REVOLUTION_IMPLICIT_OBJECT__

#include <PhysBAM_Tools/Grids_Uniform/UNIFORM_GRID_ITERATOR_NODE.h>
#include <PhysBAM_Tools/Images/IMAGE.h>
#include <PhysBAM_Geometry/Grids_Uniform_Level_Sets/LEVELSET_2D.h>
#include <PhysBAM_Geometry/Implicit_Objects_Uniform/LEVELSET_IMPLICIT_OBJECT.h>
namespace PhysBAM{

template<class T_input>
class SURFACE_OF_REVOLUTION_IMPLICIT_OBJECT:public IMPLICIT_OBJECT<VECTOR<T_input,3> >
{
    typedef T_input T;typedef VECTOR<T,3> TV;
    using IMPLICIT_OBJECT<TV>::box;
public:
    GRID<VECTOR<T,2> > slice_grid;
    ARRAY<T,VECTOR<int,2> > slice_phi;
    LEVELSET_2D<GRID<VECTOR<T,2> > > slice_levelset;
    T width,height,tolerance;

    SURFACE_OF_REVOLUTION_IMPLICIT_OBJECT(std::string& filename,const T width_input,const T height_input,const T tolerance_input=1e-6)
        :slice_levelset(slice_grid,slice_phi),width(width_input),height(height_input),tolerance(tolerance_input)
    {
        ARRAY<TV,VECTOR<int,2> > slice_image;IMAGE<T>::Read(filename,slice_image);
        slice_grid=GRID<VECTOR<T,2> >(slice_image.counts.x,slice_image.counts.y,0,width,0,height);slice_phi.Resize(slice_grid.Domain_Indices());
        for(typename GRID<VECTOR<T,2> >::NODE_ITERATOR iterator(slice_grid);iterator.Valid();iterator.Next()){VECTOR<int,2> node=iterator.Node_Index();
            slice_phi(node)=2*slice_image(node).Average()-1;}
        slice_levelset.Fast_Marching_Method();slice_levelset.Compute_Normals();
        Update_Box();
    }

    void Update_Box() PHYSBAM_OVERRIDE
    {box=RANGE<TV>(TV(-width,0,-width),TV(width,height,width));}

    T Integration_Step(const T phi) const PHYSBAM_OVERRIDE
    {return max(phi,tolerance);}

private:
    static VECTOR<T,2> Slice_Location(const TV& X)
    {return VECTOR<T,2>(sqrt(sqr(X.x)+sqr(X.z)),X.y);}
public:

    T operator()(const TV& X) const PHYSBAM_OVERRIDE
    {return slice_levelset.Phi(Slice_Location(X));}

    TV Normal(const TV& X,const int aggregate=-1) const PHYSBAM_OVERRIDE
    {assert((aggregate>=1 && aggregate<=6) || aggregate==-1);if(aggregate!=-1) return box.Normal(aggregate);
    VECTOR<T,2> slice_X=Slice_Location(X),horizontal_direction=VECTOR<T,2>(X.x,X.z).Normalized(); // cse should remove duplicate sqrt
    VECTOR<T,2> slice_normal=slice_levelset.Normal(slice_X);
    VECTOR<T,2> horizontal_normal=slice_normal.x*horizontal_direction;
    return TV(horizontal_normal.x,slice_normal.y,horizontal_normal.y);}

//#####################################################################
};
}
#endif
#endif
