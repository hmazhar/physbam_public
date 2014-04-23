// Copyright (c) 2011, Eftychios Sifakis.
// Distributed under the FreeBSD license (see license.txt)

#include <PhysBAM_Geometry/Geometry_Particles/REGISTER_GEOMETRY_READ_WRITE.h>
#include <PhysBAM_Geometry/Read_Write/Geometry/READ_WRITE_TETRAHEDRALIZED_VOLUME.h>
#include <PhysBAM_Geometry/Solids_Geometry/DEFORMABLE_GEOMETRY_COLLECTION.h>
#include <PhysBAM_Geometry/Topology_Based_Geometry/FREE_PARTICLES.h>
#include <PhysBAM_Geometry/Topology_Based_Geometry/TETRAHEDRALIZED_VOLUME.h>
#include <PhysBAM_Geometry/Topology_Based_Geometry/TRIANGULATED_SURFACE.h>
#include <PhysBAM_Tools/Grids_Uniform/GRID.h>
#include "SIMULATION_LAYOUT.h"

using namespace PhysBAM;

template<class T>
SIMULATION_LAYOUT<T>::SIMULATION_LAYOUT(const STREAM_TYPE stream_type_input)
    :stream_type(stream_type_input),
    youngs_modulus(1e5),collision_coefficient(1e5),poissons_ratio(.4),rayleigh_coefficient(1e-2),density(1e3),ground_height(-3.),ground_thickness(.01),
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

    // Generate a ground place

    ground_surface=TRIANGULATED_SURFACE<T>::Create(particles);
    int p1=particles.array_collection->Add_Element();particles.X(p1)=TV(-10,ground_height,-10);
    int p2=particles.array_collection->Add_Element();particles.X(p2)=TV(-10,ground_height, 10);
    int p3=particles.array_collection->Add_Element();particles.X(p3)=TV( 10,ground_height, 10);
    int p4=particles.array_collection->Add_Element();particles.X(p4)=TV( 10,ground_height,-10);
    ground_surface->mesh.elements.Append(VECTOR<int,3>(p1,p2,p3));
    ground_surface->mesh.elements.Append(VECTOR<int,3>(p1,p3,p4));

    // Update nodes of both structures (i.e. notify their meshes that extra particles, those of the surface, have been generated)

    ball_volume->Update_Number_Nodes();
    ground_surface->Update_Number_Nodes();

    const int number_of_particles=particles.array_collection->Size();

    // Visualize both ground plane, and tetrahedalized sphere

    collection.Add_Structure(ball_volume);collection.Add_Structure(ground_surface);

    // Assign mass (uniformly)

    T ball_mass=density*(4./3.)*pi;
    mass.Resize(number_of_particles);mass.Fill(ball_mass/(T)number_of_particles);

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

    for(int p=1;p<=particles.array_collection->Size();p++){
        T depth=X(p).y-ground_height-ground_thickness;
        if(depth<0)
            force(p)-=TV(0,depth*collision_coefficient,0);
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

    for(int p=1;p<=particles.array_collection->Size();p++){
        T depth=X(p).y-ground_height-ground_thickness;
        if(depth<0)
            force(p)-=TV(0,V(p).y*2.*sqrt(mass(p)*collision_coefficient),0);
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
    // Essentially, we hold the 4 vertices of the ground plane fixed!

    const int p1=particles.array_collection->Size()-3,p2=p1+1,p3=p2+1,p4=p3+1;
    X(p1)=TV(-10,ground_height,-10);
    X(p2)=TV(-10,ground_height, 10);
    X(p3)=TV( 10,ground_height, 10);
    X(p4)=TV( 10,ground_height,-10);
}

template<class T>
void SIMULATION_LAYOUT<T>::Set_Kinematic_Velocities(const T time,ARRAY_VIEW<TV> V) const
{
    // Essentially, we hold the 4 vertices of the ground plane fixed!

    const int p1=particles.array_collection->Size()-3,p2=p1+1,p3=p2+1,p4=p3+1;
    V(p1)=V(p2)=V(p3)=V(p4)=TV();
}

template<class T>
void SIMULATION_LAYOUT<T>::Clear_Values_Of_Kinematic_Particles(ARRAY_VIEW<TV> array) const
{
    // Essentially, we hold the 4 vertices of the ground plane fixed!

    const int p1=particles.array_collection->Size()-3,p2=p1+1,p3=p2+1,p4=p3+1;
    array(p1)=array(p2)=array(p3)=array(p4)=TV();
}

template class SIMULATION_LAYOUT<float>;
