//#####################################################################
// Copyright 2002-2006, Ronald Fedkiw, Geoffrey Irving, Frank Losasso, Andrew Selle, Tamar Shinar.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class ADVECTION
//#####################################################################
//
// Solves Z_t + u Z_x + v Z_y + w Z_z = 0 for one Euler step. Conservative schemes solve Z_t +(uZ)_x + (vZ)_y + (wZ)_z = 0.
// Z is updated using u, v, w, & Z_ghost. Input Z as (1,m) by (1,n) by (1,mn), and Z_ghost as (-2,m+3) by (-2,n+3) by (-2,mn+3).
// Input u, v, & w as (1,m) by (1,n) by (1,mn) for nonconservative schemes. Input u, v, & w as (-2,m+3) by (-2,n+3) by (-2,mn+3) for conservative schemes.
//
//#####################################################################
#ifndef __ADVECTION__
#define __ADVECTION__

#include <PhysBAM_Tools/Advection/ADVECTION_FORWARD.h>
#include <PhysBAM_Tools/Grids_Dyadic_Boundaries/BOUNDARY_POLICY_DYADIC.h>
#include <PhysBAM_Tools/Grids_Uniform_Boundaries/BOUNDARY_POLICY_UNIFORM.h>
#include <PhysBAM_Tools/Log/DEBUG_UTILITIES.h>
#include <PhysBAM_Tools/Utilities/NONCOPYABLE.h>
#include <PhysBAM_Tools/Utilities/PHYSBAM_OVERRIDE.h>
namespace PhysBAM{

template<class T_GRID> struct BOUNDARY_POLICY;
template<class T_GRID> struct GRID_ARRAYS_POLICY;

template<class T_GRID,class T2,class T_FACE_LOOKUP> // T_FACE_LOOKUP=typename T_GRID::FACE_LOOKUP
class ADVECTION:public NONCOPYABLE
{
    typedef typename T_GRID::VECTOR_T TV;typedef typename TV::SCALAR T;
    typedef typename GRID_ARRAYS_POLICY<T_GRID>::ARRAYS_SCALAR T_ARRAYS_SCALAR;typedef typename GRID_ARRAYS_POLICY<T_GRID>::FACE_ARRAYS T_FACE_ARRAYS_SCALAR;
    typedef typename REBIND<T_ARRAYS_SCALAR,T2>::TYPE T_ARRAYS_T2;typedef typename REBIND<T_ARRAYS_SCALAR,TV>::TYPE T_ARRAYS_VECTOR;
    typedef typename BOUNDARY_POLICY<T_GRID>::BOUNDARY_SCALAR T_BOUNDARY;typedef typename REBIND<T_BOUNDARY,T2>::TYPE T_BOUNDARY_T2;
public:

    ADVECTION()
    {}

    virtual ~ADVECTION()
    {}

    void Update_Advection_Equation_Cell(const T_GRID& grid,T_ARRAYS_T2& Z,const T_ARRAYS_T2& Z_ghost,
        const T_FACE_ARRAYS_SCALAR& face_velocities,T_BOUNDARY_T2& boundary,const T dt,const T time,
        const T_ARRAYS_T2* Z_min_ghost=0,const T_ARRAYS_T2* Z_max_ghost=0,T_ARRAYS_T2* Z_min=0,T_ARRAYS_T2* Z_max=0)
    {Update_Advection_Equation_Cell_Lookup(grid,Z,Z_ghost,T_FACE_LOOKUP(face_velocities),boundary,dt,time,Z_min_ghost,Z_max_ghost,Z_min,Z_max);}

    void Update_Advection_Equation_Face(const T_GRID& grid,T_FACE_ARRAYS_SCALAR& Z,const T_FACE_ARRAYS_SCALAR& Z_ghost,
        const T_FACE_ARRAYS_SCALAR& face_velocities,T_BOUNDARY& boundary,const T dt,const T time,
        const T_FACE_ARRAYS_SCALAR* Z_min_ghost=0,const T_FACE_ARRAYS_SCALAR* Z_max_ghost=0,T_FACE_ARRAYS_SCALAR* Z_min=0,T_FACE_ARRAYS_SCALAR* Z_max=0)
    {if(Z_min_ghost && Z_max_ghost){
        T_FACE_LOOKUP Z_min_lookup(*Z_min_ghost),Z_max_lookup(*Z_max_ghost);
        Update_Advection_Equation_Face_Lookup(grid,Z,T_FACE_LOOKUP(Z_ghost),T_FACE_LOOKUP(face_velocities),boundary,dt,time,&Z_min_lookup,&Z_max_lookup,Z_min,Z_max);}
    else Update_Advection_Equation_Face_Lookup(grid,Z,T_FACE_LOOKUP(Z_ghost),T_FACE_LOOKUP(face_velocities),boundary,dt,time,0,0,0,0);}

//#####################################################################
    virtual void Update_Advection_Equation_Node(const T_GRID& grid,T_ARRAYS_T2& Z,const T_ARRAYS_T2& Z_ghost,
        const T_ARRAYS_VECTOR& V,T_BOUNDARY_T2& boundary,const T dt,const T time,
        const T_ARRAYS_T2* Z_min_ghost=0,const T_ARRAYS_T2* Z_max_ghost=0,T_ARRAYS_T2* Z_min=0,T_ARRAYS_T2* Z_max=0)
    {PHYSBAM_FUNCTION_IS_NOT_DEFINED();}
    virtual void Update_Advection_Equation_Cell(const T_GRID& grid,T_ARRAYS_T2& Z,const T_ARRAYS_T2& Z_ghost,
        const T_ARRAYS_VECTOR& V,T_BOUNDARY_T2& boundary,const T dt,const T time,
        const T_ARRAYS_T2* Z_min_ghost=0,const T_ARRAYS_T2* Z_max_ghost=0,T_ARRAYS_T2* Z_min=0,T_ARRAYS_T2* Z_max=0)
    {PHYSBAM_FUNCTION_IS_NOT_DEFINED();}
    virtual void Update_Advection_Equation_Cell_Lookup(const T_GRID& grid,T_ARRAYS_T2& Z,const T_ARRAYS_T2& Z_ghost,
        const T_FACE_LOOKUP& face_velocities,T_BOUNDARY_T2& boundary,const T dt,const T time,
        const T_ARRAYS_T2* Z_min_ghost,const T_ARRAYS_T2* Z_max_ghost,T_ARRAYS_T2* Z_min,T_ARRAYS_T2* Z_max)
    {PHYSBAM_FUNCTION_IS_NOT_DEFINED();}
    virtual void Update_Advection_Equation_Face_Lookup(const T_GRID& grid,T_FACE_ARRAYS_SCALAR& Z,const T_FACE_LOOKUP& Z_ghost,
        const T_FACE_LOOKUP& face_velocities,T_BOUNDARY& boundary,const T dt,const T time,
        const T_FACE_LOOKUP* Z_min_ghost,const T_FACE_LOOKUP* Z_max_ghost,T_FACE_ARRAYS_SCALAR* Z_min,T_FACE_ARRAYS_SCALAR* Z_max)
    {PHYSBAM_FUNCTION_IS_NOT_DEFINED();}
//#####################################################################
};
}
#endif
