// Copyright (c) 2011, Taylor Patterson, Eftychios Sifakis.
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

    const int number_of_particles=layout.particles.array_collection->Size();
	Vp=layout.particles.V;
#if 1
        for(int p=1;p<=number_of_particles;p++)
            layout.particles.X(p)+=dt*layout.particles.V(p);
#else
        for(int p=1;p<=number_of_particles;p++)
            layout.particles.V(p)=TV();
#endif

	for(int j=1;j<=10;j++)
	    Simulate_Time_Step(time,dt);}
}

template<class T>
void SIMULATION_DRIVER<T>::Simulate_Time_Step(const T time,const T dt)
{
    layout.Set_Kinematic_Positions(time+dt,layout.particles.X);
    layout.Set_Kinematic_Velocities(time+dt,layout.particles.V);

    // Construct right-hand-side
    const int number_of_particles=layout.particles.array_collection->Size();
    ARRAY<TV> rhs(number_of_particles);
    layout.Add_Elastic_Forces(layout.particles.X,rhs);
    layout.Add_Damping_Forces(layout.particles.X,layout.particles.V,rhs);
    layout.Add_External_Forces(rhs);

    for(int p=1;p<=number_of_particles;p++)
        rhs(p)+=(layout.mass(p)*(Vp(p)-layout.particles.V(p)))/dt;

    // Use (0,0,0) as initial guess for change in positions
    ARRAY<TV> delta_X(number_of_particles);
    // not necessary: delta_X.Fill(TV(0.,0.,0.));

    // Temporary vectors, required by Conjugate Gradients
    ARRAY<TV> temp_q(number_of_particles),temp_s(number_of_particles),temp_r(number_of_particles),
        temp_k(number_of_particles),temp_z(number_of_particles);

    // Encapsulate all vectors in CG-mandated format
    CG_VECTOR<T> cg_x(delta_X),cg_b(rhs),
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
        1e-3,0,100);

    // Remove the solved change in position for kinematic particles
    layout.Clear_Values_Of_Kinematic_Particles(delta_X);

    // Update particle positions and velocities
    layout.particles.X+=delta_X;
    layout.particles.V+=delta_X/dt;

//    layout.Set_Kinematic_Positions(time+dt,layout.particles.X);
//    layout.Set_Kinematic_Velocities(time+dt,layout.particles.V);
}

template class SIMULATION_DRIVER<float>;

