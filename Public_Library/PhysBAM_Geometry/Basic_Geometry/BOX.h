//#####################################################################
// Copyright 2002-2007, Robert Bridson, Ronald Fedkiw, Eran Guendelman, Geoffrey Irving, Sergey Koltakov, Nipun Kwatra, Neil Molino, Igor Neverov, Duc Nguyen, Avi Robinson-Mosher,
//     Craig Schroeder, Andrew Selle, Tamar Shinar, Jonathan Su, Jerry Talton, Joseph Teran.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class BOX
//#####################################################################
#ifndef __BOX__
#define __BOX__

#include <PhysBAM_Tools/Log/DEBUG_UTILITIES.h>
#include <PhysBAM_Tools/Math_Tools/RANGE.h>
#include <PhysBAM_Tools/Utilities/TYPE_UTILITIES.h>
#include <PhysBAM_Geometry/Basic_Geometry/BASIC_GEOMETRY_FORWARD.h>
namespace PhysBAM{

template<class TV> struct IS_SCALAR_BLOCK<BOX<TV> >:public IS_SCALAR_BLOCK<TV>{};
template<class TV,class RW> struct IS_BINARY_IO_SAFE<BOX<TV>,RW> {static const bool value=false;}; // required since memory format differs from disk format

template<class TV>
class BOX:public RANGE<TV>
{
    typedef typename TV::SCALAR T;
    STATIC_ASSERT(IS_FLOAT_OR_DOUBLE<T>::value); // if you want a box of ints, use RANGE
public:
    template<class T2> struct REBIND{typedef BOX<T2> TYPE;};
    typedef T SCALAR;
    typedef TV VECTOR_T;
    enum WORKAROUND {d=TV::dimension};

    typedef RANGE<TV> BASE;
    using BASE::min_corner;using BASE::max_corner;using BASE::Edge_Lengths;using BASE::Center;using BASE::Size;using BASE::Intersection;
    using BASE::Lazy_Outside;using BASE::Bounding_Box;using BASE::Inside;

    BOX()
    {}

    BOX(const T xmin,const T xmax)
        :RANGE<TV>(xmin,xmax)
    {}

    BOX(const T xmin,const T xmax,const T ymin,const T ymax)
        :RANGE<TV>(xmin,xmax,ymin,ymax)
    {}

    BOX(const T xmin,const T xmax,const T ymin,const T ymax,const T zmin,const T zmax)
        :RANGE<TV>(xmin,xmax,ymin,ymax,zmin,zmax)
    {}

    BOX(const TV& minimum_corner,const TV& maximum_corner)
        :RANGE<TV>(minimum_corner,maximum_corner)
    {}

    BOX(const RANGE<TV>& box)
        :RANGE<TV>(box)
    {}

    template<class T2> explicit BOX(const RANGE<T2>& box)
        :RANGE<TV>(box)
    {}

    BOX(const TV& point)
        :RANGE<TV>(point)
    {}

    BOX(const BOX<TV>& box,const FRAME<VECTOR<T,1> >& frame) // allow 1d boxes to be used as oriented boxes
        :RANGE<TV>(frame*box.min_corner,frame*box.max_corner)
    {}

    BOX<TV> Axis_Aligned_Bounding_Box() const
    {return *this;}

    const RANGE<TV>& Bounding_Box() const // for templatization purposes
    {return *this;}

    VECTOR<T,TV::dimension-1> Principal_Curvatures(const TV& X) const
    {return VECTOR<T,TV::dimension-1>();}

//#####################################################################
    TV Normal(const int aggregate) const;
    TV Normal(const TV& X) const;
    TV Surface(const TV& X) const;
    T Signed_Distance(const TV& X) const;
    static std::string Name();
//#####################################################################
};
}
#endif
