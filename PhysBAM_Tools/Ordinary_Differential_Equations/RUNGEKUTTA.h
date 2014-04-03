//#####################################################################
// Copyright 2002-2010, Ronald Fedkiw, Geoffrey Irving, Michael Lentine, Frank Losasso, Duc Nguyen, Eran Guendelman.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class RUNGEKUTTA
//#####################################################################
#ifndef __RUNGEKUTTA__
#define __RUNGEKUTTA__

#include <PhysBAM_Tools/Arrays/ARRAY.h>
#include <PhysBAM_Tools/Arrays/Scalar_View.h>
#include <PhysBAM_Tools/Log/DEBUG_UTILITIES.h>
#include <PhysBAM_Tools/Utilities/NONCOPYABLE.h>
#include <PhysBAM_Tools/Vectors/SCALAR_POLICY.h>
namespace PhysBAM{

//#####################################################################
// Class RUNGEKUTTA_BOUNDARY_CONDITION_HELPER
//#####################################################################
template<class TV> class RUNGEKUTTA;

template<class T>
class RUNGEKUTTA_BOUNDARY_CONDITION_HELPER_BASE
{
public:
    virtual ~RUNGEKUTTA_BOUNDARY_CONDITION_HELPER_BASE() {};
    virtual void operator()(ARRAY_VIEW<T> u_view, const T time)=0;
};

template<class TV, class T_GRID, class T_BOUNDARY>
class RUNGEKUTTA_BOUNDARY_CONDITION_HELPER:public RUNGEKUTTA_BOUNDARY_CONDITION_HELPER_BASE<typename TV::SCALAR>
{
    typedef typename TV::SCALAR T;
public:
    RUNGEKUTTA<TV>& instance;
    T_GRID& grid;
    T_BOUNDARY& boundary;

    RUNGEKUTTA_BOUNDARY_CONDITION_HELPER(RUNGEKUTTA<TV>& _instance, T_GRID& _grid, T_BOUNDARY& _boundary):
        instance(_instance),grid(_grid),boundary(_boundary) {};
    ~RUNGEKUTTA_BOUNDARY_CONDITION_HELPER() {};
    void operator()(ARRAY_VIEW<T> u_view, const T time)
    {
        instance.Enforce_Boundary_Condition(grid,boundary,u_view,time);
    }
};
    
//#####################################################################
// Class RUNGEKUTTA_CORE
//#####################################################################
template<class T>
class RUNGEKUTTA_CORE:public NONCOPYABLE
{
public:
    int order; // 1, 2 or 3
    T dt; // size of time step
    int substep; // current substep
    ARRAY_VIEW<T> u; // variable being advanced in time
    ARRAY<T> u_copy; // a copy of the variable
protected:
    T time; // current time
    int pseudo_time; // real time or pseudo time - for boundary conditions
    RUNGEKUTTA_BOUNDARY_CONDITION_HELPER_BASE<T>* enforce_boundary_condition;
public:

    template<class TV>
    RUNGEKUTTA_CORE(TV& u_input)
        :u(Scalar_View(u_input)),u_copy(u),enforce_boundary_condition(0)
    {
        Set_Order();Set_Time();Real_Time();
    }
    
    ~RUNGEKUTTA_CORE()
    {
        if(enforce_boundary_condition)
        {
            delete enforce_boundary_condition;
            enforce_boundary_condition=0;
        }
    }

    void Set_Order(const int runge_kutta_order=3)
    {order=runge_kutta_order;}

    void Set_Time(const T starting_time=0)
    {time=starting_time;}

    void Real_Time()
    {pseudo_time=0;}

    void Pseudo_Time()
    {pseudo_time=1;}

//#####################################################################
    void Start(const T dt_input);
    T Main();
//#####################################################################
};

//#####################################################################
// Class RUNGEKUTTA
//#####################################################################
template<class TV>
class RUNGEKUTTA:public RUNGEKUTTA_CORE<typename SCALAR_POLICY<TV>::TYPE>
{
    typedef typename SCALAR_POLICY<TV>::TYPE T;
private:
    TV& u;
public:

    RUNGEKUTTA(TV& u)
        :RUNGEKUTTA_CORE<T>(u),u(u)
    {
        STATIC_ASSERT(SCALAR_VIEW_IS_VECTOR_SPACE<TV>::value);
    }

    static RUNGEKUTTA* Create(TV& u,const int order,const T dt,const T time)
    {RUNGEKUTTA* runge_kutta=new RUNGEKUTTA(u);runge_kutta->Set_Order(order);runge_kutta->Set_Time(time);
    runge_kutta->dt=dt;runge_kutta->substep=0;return runge_kutta;}

    template<class T_GRID,class T_BOUNDARY>
    void Set_Grid_And_Boundary_Condition(T_GRID& grid,T_BOUNDARY& boundary)
    {
        if(this->enforce_boundary_condition)
            delete this->enforce_boundary_condition;
        this->enforce_boundary_condition=new RUNGEKUTTA_BOUNDARY_CONDITION_HELPER<TV,T_GRID,T_BOUNDARY>(*this,grid,boundary);
    }

// private:
    template<class T_GRID,class T_BOUNDARY> void
    Enforce_Boundary_Condition(const T_GRID& grid,T_BOUNDARY& boundary,ARRAY_VIEW<T> u_view,const T time)
    {ARRAY_VIEW<T> real_u_view(Scalar_View(u));PHYSBAM_ASSERT(real_u_view.Size()==u_view.Size() && real_u_view.Get_Array_Pointer()==u_view.Get_Array_Pointer());
    boundary.Apply_Boundary_Condition(grid,u,time);}
};
template<class T,int d>
class RUNGEKUTTA<ARRAY<T,FACE_INDEX<d> > >:public RUNGEKUTTA_CORE<T>
{
private:
    ARRAY<T,FACE_INDEX<d> >& u;
public:

    RUNGEKUTTA(ARRAY<T,FACE_INDEX<d> >& u)
        :RUNGEKUTTA_CORE<T>(u),u(u)
    {}

    static RUNGEKUTTA* Create(ARRAY<T,FACE_INDEX<d> >& u,const int order,const T dt,const T time)
    {RUNGEKUTTA* runge_kutta=new RUNGEKUTTA(u);runge_kutta->Set_Order(order);runge_kutta->Set_Time(time);
    runge_kutta->dt=dt;runge_kutta->substep=0;return runge_kutta;}

    template<class T_GRID,class T_BOUNDARY>
    void Set_Grid_And_Boundary_Condition(T_GRID& grid,T_BOUNDARY& boundary)
    {
        if(this->enforce_boundary_condition)
            delete this->enforce_boundary_condition;
        this->enforce_boundary_condition=new RUNGEKUTTA_BOUNDARY_CONDITION_HELPER<ARRAY<T,FACE_INDEX<d> >,T_GRID,T_BOUNDARY>(*this,grid,boundary);
    }

// private:
    template<class T_GRID,class T_BOUNDARY> void
    Enforce_Boundary_Condition(const T_GRID& grid,T_BOUNDARY& boundary,ARRAY_VIEW<T> u_view,const T time)
    {ARRAY_VIEW<T> real_u_view(Scalar_View(u));PHYSBAM_ASSERT(real_u_view.Size()==u_view.Size() && real_u_view.Get_Array_Pointer()==u_view.Get_Array_Pointer());
    boundary.Apply_Boundary_Condition_Face(grid,u,time);}
};
}
#endif
