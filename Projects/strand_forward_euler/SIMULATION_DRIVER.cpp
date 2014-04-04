// Copyright (c) 2011, Eftychios Sifakis.
// Distributed under the FreeBSD license (see license.txt)

#include <PhysBAM_Geometry/Solids_Geometry/DEFORMABLE_GEOMETRY_COLLECTION.h>

#include "SIMULATION_DRIVER.h"
#include "SIMULATION_LAYOUT.h"
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

    ARRAY<TV> force(layout.n),dX(layout.n),dV(layout.n);
    layout.Add_Elastic_Forces(layout.particles.X,force);
    layout.Add_Damping_Forces(layout.particles.X,layout.particles.V,force);
    layout.Add_External_Forces(force);

    // Apply the Forward Euler method
    dX=dt*layout.particles.V;               // Compute position change
    for(int p=1;p<=layout.n;p++)
        dV(p)=(dt/layout.mass(p))*force(p); // Compute velocity change

    layout.Clear_Values_Of_Kinematic_Particles(dX);
    layout.Clear_Values_Of_Kinematic_Particles(dV);

    layout.particles.X+=dX;                 // Update particle positions and velocities
    layout.particles.V+=dV;

    layout.Set_Kinematic_Positions(time+dt,layout.particles.X);
    layout.Set_Kinematic_Velocities(time+dt,layout.particles.V);
}

template class SIMULATION_DRIVER<float>;
