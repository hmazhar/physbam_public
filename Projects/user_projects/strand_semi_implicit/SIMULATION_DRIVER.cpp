// Copyright (c) 2011, Eftychios Sifakis.
// Distributed under the FreeBSD license (see license.txt)

#include <PhysBAM_Tools/Krylov_Solvers/CONJUGATE_GRADIENT.h>

#include "SIMULATION_DRIVER.h"
#include "CG_VECTOR.h"
#include "CG_SYSTEM.h"
using namespace PhysBAM;

template<class T>
void SIMULATION_DRIVER<T>::Simulate_Frame(const int frame)
{
    T frame_end_time=layout.frame_time*(T)frame,dt_max,dt;

    for(;time<frame_end_time;time+=dt){
        dt_max=layout.Maximum_Dt();
        dt=std::min(dt_max,(T)1.001*(frame_end_time-time));
        Simulate_Time_Step(time,dt);}
}

template<class T>
void SIMULATION_DRIVER<T>::Simulate_Time_Step(const T time,const T dt)
{
    layout.Set_Kinematic_Positions(time,layout.particles.X);
    layout.Set_Kinematic_Velocities(time,layout.particles.V);

    // Construct right-hand-side
    ARRAY<TV> rhs(layout.n);
    layout.Add_Elastic_Forces(layout.particles.X,rhs);
    layout.Add_External_Forces(rhs);
    for(int p=1;p<=layout.n;p++)
        rhs(p)=layout.mass(p)*layout.particles.V(p)+dt*rhs(p);

    // Use previous velocities as initial guess for next velocities
    ARRAY<TV> V_next(layout.particles.V);

    // Temporary vectors, required by Conjugate Gradients
    ARRAY<TV> temp_q(layout.n),temp_s(layout.n),temp_r(layout.n),
        temp_k(layout.n),temp_z(layout.n);

    // Encapsulate all vectors in CG-mandated format
    CG_VECTOR<T> cg_x(V_next),cg_b(rhs),
        cg_q(temp_q),cg_s(temp_s),cg_r(temp_r),cg_k(temp_k),cg_z(temp_z);

    // Generate CG-formatted system object
    CG_SYSTEM<T> cg_system(layout,time,dt);

    // Generate Conjugate Gradients solver object
    CONJUGATE_GRADIENT<T> cg;
    cg.print_residuals=true;
    cg.print_diagnostics=true;

    // Solve linear system using CG
    cg.Solve(cg_system,
        cg_x,cg_b,cg_q,cg_s,cg_r,cg_k,cg_z,
        1e-6,0,100);

    // Trapezoidal rule for positions
    ARRAY<TV> dX(layout.n);
    dX =(dt/2)*layout.particles.V;
    dX+=(dt/2)*V_next;
    layout.Clear_Values_Of_Kinematic_Particles(dX);

    layout.particles.X+=dX;    // Update particle positions and velocities
    layout.particles.V=V_next;

    layout.Set_Kinematic_Positions(time+dt,layout.particles.X);
    layout.Set_Kinematic_Velocities(time+dt,layout.particles.V);
}

template class SIMULATION_DRIVER<float>;

