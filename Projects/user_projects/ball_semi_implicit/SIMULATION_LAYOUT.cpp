// Copyright (c) 2011, Eftychios Sifakis.
// Distributed under the FreeBSD license (see license.txt)

#include <PhysBAM_Geometry/Geometry_Particles/REGISTER_GEOMETRY_READ_WRITE.h>
#include <PhysBAM_Geometry/Read_Write/Geometry/READ_WRITE_TETRAHEDRALIZED_VOLUME.h>
#include <PhysBAM_Geometry/Solids_Geometry/DEFORMABLE_GEOMETRY_COLLECTION.h>
#include <PhysBAM_Geometry/Topology_Based_Geometry/FREE_PARTICLES.h>
#include <PhysBAM_Geometry/Topology_Based_Geometry/TETRAHEDRALIZED_VOLUME.h>
#include <PhysBAM_Tools/Grids_Uniform/GRID.h>
#include "SIMULATION_LAYOUT.h"

using namespace PhysBAM;

template<class T>
SIMULATION_LAYOUT<T>::SIMULATION_LAYOUT(const STREAM_TYPE stream_type_input)
    :stream_type(stream_type_input),
    youngs_modulus(1e5),poissons_ratio(.4),rayleigh_coefficient(2e-2),density(1e3),
    collection(particles),number_of_frames(200),frame_time(.05)
    // ,CFL_number(0.03)
{
    lambda=youngs_modulus*poissons_ratio/((1+poissons_ratio)*(1-2*poissons_ratio));
    mu=youngs_modulus/(2*(1+poissons_ratio));
    alpha=rayleigh_coefficient*lambda;
    beta=rayleigh_coefficient*mu;

    Initialize_Geometry_Particle();Initialize_Read_Write_Structures();
}

template<class T>
void SIMULATION_LAYOUT<T>::Initialize()
{
    particles.Store_Velocity();

    // Read a tetrahedralized sphere model from a file

    ball_volume=TETRAHEDRALIZED_VOLUME<T>::Create(particles);
    FILE_UTILITIES::Read_From_File(stream_type,"sphere.tet",*ball_volume);
    LOG::cout<<"Initializing sphere geometry ..."<<std::endl;
    LOG::cout<<"  Particles="<<ball_volume->particles.array_collection->Size()<<std::endl;
    LOG::cout<<"  Tetrahedrons="<<ball_volume->mesh.elements.m<<std::endl;

    // Mark particles which are kinematically animated

    const int number_of_particles=particles.array_collection->Size();

    constrained_particles=FREE_PARTICLES<TV>::Create(particles);
    for(int p=1;p<=number_of_particles;p++)
        if(particles.X(p).y>.9)
            constrained_particles->nodes.Append(p);

    // Visualize both constraint particles, and tetrahedalized sphere

    collection.Add_Structure(ball_volume);collection.Add_Structure(constrained_particles);

    // Assign mass (uniformly)

    T ball_mass=density*(4./3.)*pi;
    mass.Resize(number_of_particles);mass.Fill(ball_mass/(T)number_of_particles);

    // Save reference positions of all particles (used by Set_Kinematic_Positions/Velocities)
    X_reference=particles.X;

    // Precompute Dm_inverse

    Dm_inverse.Resize(ball_volume->mesh.elements.m);
    tet_rest_volume.Resize(ball_volume->mesh.elements.m);
    for(int tet=1;tet<=ball_volume->mesh.elements.m;tet++){
        int p0,p1,p2,p3;ball_volume->mesh.elements(tet).Get(p0,p1,p2,p3);
        const TV& X0=ball_volume->particles.X(p0);
        const TV& X1=ball_volume->particles.X(p1);
        const TV& X2=ball_volume->particles.X(p2);
        const TV& X3=ball_volume->particles.X(p3);
        MATRIX<T,3> Dm(X1-X0,X2-X0,X3-X0);
        Dm_inverse(tet)=Dm.Inverse();
        tet_rest_volume(tet)=fabs(Dm.Determinant())/6.;
    }
}

template<class T>
void SIMULATION_LAYOUT<T>::Add_Elastic_Forces(ARRAY_VIEW<const TV> X,ARRAY_VIEW<TV> force) const
{
    for(int tet=1;tet<=ball_volume->mesh.elements.m;tet++){
        int p0,p1,p2,p3;ball_volume->mesh.elements(tet).Get(p0,p1,p2,p3);
        const TV& x0=X(p0);
        const TV& x1=X(p1);
        const TV& x2=X(p2);
        const TV& x3=X(p3);
        MATRIX<T,3> Ds(x1-x0,x2-x0,x3-x0);
        MATRIX<T,3> F=Ds*Dm_inverse(tet);
        MATRIX<T,3> epsilon=(F+F.Transposed())*.5-MATRIX<T,3>::Identity_Matrix();
        MATRIX<T,3> P=epsilon*2.*mu+MATRIX<T,3>::Identity_Matrix()*lambda*epsilon.Trace();
        MATRIX<T,3> G=-P*Dm_inverse(tet).Transposed()*tet_rest_volume(tet);
        force(p1)+=G.Column(1);
        force(p2)+=G.Column(2);
        force(p3)+=G.Column(3);
        force(p0)+=(-G.Column(1)-G.Column(2)-G.Column(3));
    }
}

template<class T>
void SIMULATION_LAYOUT<T>::Add_Damping_Forces(ARRAY_VIEW<const TV> X,ARRAY_VIEW<const TV> V,
    ARRAY_VIEW<TV> force) const
{
    for(int tet=1;tet<=ball_volume->mesh.elements.m;tet++){
        int p0,p1,p2,p3;ball_volume->mesh.elements(tet).Get(p0,p1,p2,p3);
        const TV& v0=V(p0);
        const TV& v1=V(p1);
        const TV& v2=V(p2);
        const TV& v3=V(p3);
        MATRIX<T,3> Ds_dot(v1-v0,v2-v0,v3-v0);
        MATRIX<T,3> F_dot=Ds_dot*Dm_inverse(tet);
        MATRIX<T,3> epsilon_dot=(F_dot+F_dot.Transposed())*.5;
        MATRIX<T,3> Pd=epsilon_dot*2.*alpha+MATRIX<T,3>::Identity_Matrix()*beta*epsilon_dot.Trace();
        MATRIX<T,3> Gd=-Pd*Dm_inverse(tet).Transposed()*tet_rest_volume(tet);
        force(p1)+=Gd.Column(1);
        force(p2)+=Gd.Column(2);
        force(p3)+=Gd.Column(3);
        force(p0)+=(-Gd.Column(1)-Gd.Column(2)-Gd.Column(3));
    }
}

template<class T>
void SIMULATION_LAYOUT<T>::Add_External_Forces(ARRAY_VIEW<TV> force) const
{
    const int number_of_particles=particles.array_collection->Size();
    for(int p=1;p<=number_of_particles;p++)
        force(p)-=TV::Axis_Vector(2)*mass(p)*9.81;
}

template<class T>
T SIMULATION_LAYOUT<T>::Maximum_Dt() const
{
    // This example does not use a CFL number (we instead specify at dt, directly)
    return 1e-3;
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
    T start_time=1.;
    T end_time=5.;
    T period=2.;
    TV dX;
    if (time>=start_time && time<=end_time)
        dX=TV(1-cos(2.*pi*(time-start_time)/period),0,0);

    for(int i=1;i<=constrained_particles->nodes.m;i++){
        int p=constrained_particles->nodes(i);
        X(p)=X_reference(p)+dX;
    }
}

template<class T>
void SIMULATION_LAYOUT<T>::Set_Kinematic_Velocities(const T time,ARRAY_VIEW<TV> V) const
{
    T start_time=1.;
    T end_time=5.;
    T period=2.;
    TV dV;
    if (time>=start_time && time<=end_time)
        dV=TV(2*pi*sin(2.*pi*(time-start_time)/period)/period,0,0);

    for(int i=1;i<=constrained_particles->nodes.m;i++){
        int p=constrained_particles->nodes(i);
        V(p)=dV;
    }
}

template<class T>
void SIMULATION_LAYOUT<T>::Clear_Values_Of_Kinematic_Particles(ARRAY_VIEW<TV> array) const
{
    for(int i=1;i<=constrained_particles->nodes.m;i++){
        int p=constrained_particles->nodes(i);
        array(p)=TV();
    }
}

template class SIMULATION_LAYOUT<float>;
