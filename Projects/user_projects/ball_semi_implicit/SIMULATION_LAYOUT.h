// Copyright (c) 2011, Eftychios Sifakis.
// Distributed under the FreeBSD license (see license.txt)

#pragma once

#include <PhysBAM_Tools/Matrices/MATRIX_3X3.h>

namespace PhysBAM{

template<class T> class TETRAHEDRALIZED_VOLUME;
template<class TV> class FREE_PARTICLES;

template<class T>
class SIMULATION_LAYOUT
{
public:
    typedef VECTOR<T,3> TV;

    const STREAM_TYPE stream_type;

    const T youngs_modulus;      // Elasticity and damping coefficients
    const T poissons_ratio;
    const T rayleigh_coefficient;
    const T density;

    T mu,lambda;                 // Lame parameters
    T alpha,beta;
    ARRAY<T> mass;               // Mass (per each particle)    

    GEOMETRY_PARTICLES<TV> particles;
    DEFORMABLE_GEOMETRY_COLLECTION<TV> collection;

    const int number_of_frames;  // Total number of frames
    const T frame_time;          // Frame (snapshot) interval

    // This example does not use a CFL number (we instead specify at dt, directly)
    // const T CFL_number;          // CFL number (not to exceed 1)

private:
    TETRAHEDRALIZED_VOLUME<T>* ball_volume;
    FREE_PARTICLES<TV>* constrained_particles;
    ARRAY<TV> X_reference;
    ARRAY<MATRIX<T,3> > Dm_inverse;
    ARRAY<T> tet_rest_volume;
    
public:
    SIMULATION_LAYOUT(const STREAM_TYPE stream_type_input);
    void Initialize();
    void Add_Elastic_Forces(ARRAY_VIEW<const TV> X,ARRAY_VIEW<TV> force) const;
    void Add_Damping_Forces(ARRAY_VIEW<const TV> X,ARRAY_VIEW<const TV> V,
        ARRAY_VIEW<TV> force) const;
    void Add_External_Forces(ARRAY_VIEW<TV> force) const;
    T Maximum_Dt() const;
    void Write_Output(const int frame) const;
    void Set_Kinematic_Positions(const T time,ARRAY_VIEW<TV> X) const;
    void Set_Kinematic_Velocities(const T time,ARRAY_VIEW<TV> V) const;
    void Clear_Values_Of_Kinematic_Particles(ARRAY_VIEW<TV> array) const;
};

}
