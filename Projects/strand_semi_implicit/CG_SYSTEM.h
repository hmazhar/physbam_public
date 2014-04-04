// Copyright (c) 2011, Eftychios Sifakis.
// Distributed under the FreeBSD license (see license.txt)

#pragma once

#include <PhysBAM_Geometry/Solids_Geometry/DEFORMABLE_GEOMETRY_COLLECTION.h>
#include <PhysBAM_Tools/Krylov_Solvers/KRYLOV_SYSTEM_BASE.h>

#include "SIMULATION_LAYOUT.h"
#include "CG_VECTOR.h"

namespace PhysBAM{

template<class T>
class CG_SYSTEM:public KRYLOV_SYSTEM_BASE<T>
{
    typedef KRYLOV_SYSTEM_BASE<T> BASE;
    typedef KRYLOV_VECTOR_BASE<T> VECTOR_BASE;

    SIMULATION_LAYOUT<T>& layout;
    const T time;
    const T dt;

public:
    CG_SYSTEM(SIMULATION_LAYOUT<T>& layout_input,const T time_input,const T dt_input)
        :BASE(false,false),layout(layout_input),time(time_input),dt(dt_input) {}

    void Multiply(const VECTOR_BASE& v,VECTOR_BASE& result) const
    {
        const ARRAY<VECTOR<T,3> >& v_array=CG_VECTOR<T>::Array(v);
        ARRAY<VECTOR<T,3> >& result_array=CG_VECTOR<T>::Array(result);

        result_array.Fill(VECTOR<T,3>());
        layout.Add_Damping_Forces(layout.particles.X,v_array,result_array);
        for(int p=1;p<=v_array.m;p++)
            result_array(p)=layout.mass(p)*v_array(p)-dt*result_array(p);
    }

    double Inner_Product(const VECTOR_BASE& x,const VECTOR_BASE& y) const
    {
        const ARRAY<VECTOR<T,3> >& x_array=CG_VECTOR<T>::Array(x);
        const ARRAY<VECTOR<T,3> >& y_array=CG_VECTOR<T>::Array(y);

        double result=0.;
        for(int i=1;i<=x_array.m;i++)
            result+=VECTOR<T,3>::Dot_Product(x_array(i),y_array(i));
        return result;
    }

    T Convergence_Norm(const VECTOR_BASE& x) const
    {
        const ARRAY<VECTOR<T,3> >& x_array=CG_VECTOR<T>::Array(x);

        T result=0.;
        for(int i=1;i<=x_array.m;i++)
            result=std::max(result,x_array(i).Magnitude());
        return result;
    }

    void Project(VECTOR_BASE& x) const
    {
        ARRAY<VECTOR<T,3> >& x_array=CG_VECTOR<T>::Array(x);

        layout.Clear_Values_Of_Kinematic_Particles(x_array);
    }

    void Set_Boundary_Conditions(VECTOR_BASE& v) const
    {
        ARRAY<VECTOR<T,3> >& v_array=CG_VECTOR<T>::Array(v);

        layout.Set_Kinematic_Velocities(time+dt,v_array);
    }
    
    void Project_Nullspace(VECTOR_BASE& x) const {} // Just a stub (for solids)
};

}


