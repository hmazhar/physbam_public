// Copyright (c) 2011, Eftychios Sifakis.
// Distributed under the FreeBSD license (see license.txt)

#include <PhysBAM_Geometry/Geometry_Particles/REGISTER_GEOMETRY_READ_WRITE.h>
#include <PhysBAM_Geometry/Solids_Geometry/DEFORMABLE_GEOMETRY_COLLECTION.h>
#include <PhysBAM_Geometry/Topology_Based_Geometry/SEGMENTED_CURVE.h>
#include <PhysBAM_Geometry/Topology_Based_Geometry/FREE_PARTICLES.h>
#include "SIMULATION_LAYOUT.h"
using namespace PhysBAM;

template<class T>
SIMULATION_LAYOUT<T>::SIMULATION_LAYOUT(const STREAM_TYPE stream_type_input)
    :stream_type(stream_type_input),
     n(11),youngs_modulus(100.),damping_coefficient(10.),
    wire_mass(1.),wire_restlength(1.),
    collection(particles),
    number_of_frames(100),frame_time(.05),CFL_number(0.03)
{
    Initialize_Geometry_Particle();Initialize_Read_Write_Structures();
}

template<class T>
void SIMULATION_LAYOUT<T>::Initialize()
{
    particles.Store_Velocity();

    wire_curve=SEGMENTED_CURVE<TV>::Create(particles);
    wire_curve->mesh.Initialize_Straight_Mesh(n);particles.array_collection->Add_Elements(n);
    for(int p=1;p<=n;p++) particles.X(p)=TV(0,(T)(1-p)/(T)(n-1),.5);

    wire_particles=FREE_PARTICLES<TV>::Create(particles);
    for(int p=1;p<=n;p++) wire_particles->nodes.Append(p);

    collection.Add_Structure(wire_curve);collection.Add_Structure(wire_particles);

    mass.Resize(n);mass.Fill(wire_mass/(T)n);
    restlength.Resize(n-1);restlength.Fill(wire_restlength/(T)(n-1));
}

template<class T>
void SIMULATION_LAYOUT<T>::Add_Elastic_Forces(ARRAY_VIEW<const TV> X,ARRAY_VIEW<TV> force) const
{
    for(int s=1;s<=wire_curve->mesh.elements.m;s++){
        int p1,p2;wire_curve->mesh.elements(s).Get(p1,p2);
        TV X1=X(p1),X2=X(p2);
        TV normal=(X1-X2).Normalized();
        T length=(X1-X2).Magnitude();
        TV f=-normal*youngs_modulus*(length/restlength(s)-1.);
        force(p1)+=f;force(p2)-=f;
    }
}

template<class T>
void SIMULATION_LAYOUT<T>::Add_Damping_Forces(ARRAY_VIEW<const TV> X,ARRAY_VIEW<const TV> V,
    ARRAY_VIEW<TV> force) const
{
    for(int s=1;s<=wire_curve->mesh.elements.m;s++){
        int p1,p2;wire_curve->mesh.elements(s).Get(p1,p2);
        TV X1=X(p1),X2=X(p2);
        TV V1=V(p1),V2=V(p2);
        TV normal=(X1-X2).Normalized();
        T vrel=TV::Dot_Product(normal,V1-V2);
        TV f=-damping_coefficient*vrel*normal;
        force(p1)+=f;force(p2)-=f;
    }
}

template<class T>
void SIMULATION_LAYOUT<T>::Add_External_Forces(ARRAY_VIEW<TV> force) const
{
    for(int p=1;p<=n;p++)
        force(p)-=TV::Axis_Vector(2)*mass(p)*9.81;
}

template<class T>
T SIMULATION_LAYOUT<T>::Maximum_Dt() const
{
    T maximum_dt=FLT_MAX;
    for(int s=1;s<=wire_curve->mesh.elements.m;s++)
        maximum_dt=std::min(maximum_dt,2*damping_coefficient*restlength(s)/youngs_modulus);

    return CFL_number*maximum_dt;
}

template<class T>
void SIMULATION_LAYOUT<T>::Write_Output(const int frame) const
{
    FILE_UTILITIES::Create_Directory("output/"+STRING_UTILITIES::Value_To_String(frame));
    collection.Write(stream_type,"output",frame,0,true);
}

template<class T>
void SIMULATION_LAYOUT<T>::Set_Kinematic_Positions(const T time,ARRAY_VIEW<TV> X) const
{
    // Set upper endpoint position
    T angular_velocity=two_pi/(frame_time*(T)number_of_frames);
    X(1)=TV(.5*sin(time*angular_velocity),0,.5*cos(time*angular_velocity));
}

template<class T>
void SIMULATION_LAYOUT<T>::Set_Kinematic_Velocities(const T time,ARRAY_VIEW<TV> V) const
{
    // Set upper endpoint position
    T angular_velocity=two_pi/(frame_time*(T)number_of_frames);
    V(1)=TV(.5*angular_velocity*cos(time*angular_velocity),0,
        -.5*angular_velocity*sin(time*angular_velocity));
}

template<class T>
void SIMULATION_LAYOUT<T>::Clear_Values_Of_Kinematic_Particles(ARRAY_VIEW<TV> array) const
{
    // Only the first particle is constrained
    array(1)=TV();
}

template class SIMULATION_LAYOUT<float>;
