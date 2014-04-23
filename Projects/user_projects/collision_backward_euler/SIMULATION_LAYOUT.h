// Copyright (c) 2011, Taylor Patterson, Eftychios Sifakis.
// Distributed under the FreeBSD license (see license.txt)

#pragma once

#include <PhysBAM_Tools/Matrices/DIAGONAL_MATRIX_3X3.h>
#include <PhysBAM_Tools/Matrices/MATRIX_3X3.h>
#include "ROTATED_STRESS_DERIVATIVE.h"

namespace PhysBAM{

template<class T> class TETRAHEDRALIZED_VOLUME;
template<class T> class TRIANGULATED_SURFACE;

template<class T>
class SIMULATION_LAYOUT
{
public:
    typedef VECTOR<T,3> TV;

    const STREAM_TYPE stream_type;

    const T youngs_modulus;        // Elasticity and damping coefficients
    const T collision_coefficient;
    const T poissons_ratio;
    const T rayleigh_coefficient;
    const T density;
    const T ground_height;
    const T ground_thickness;

    T mu,lambda;                  // Lame parameters
    T alpha,beta;
    ARRAY<T> mass;                // Mass (per each particle)    

    GEOMETRY_PARTICLES<TV> particles;
    DEFORMABLE_GEOMETRY_COLLECTION<TV> collection;

    const int number_of_frames;   // Total number of frames
    const T frame_time;           // Frame (snapshot) interval

    // This example does not use a CFL number (we instead specify at dt, directly)
    // const T CFL_number;          // CFL number (not to exceed 1)

private:
    TETRAHEDRALIZED_VOLUME<T>* ball_volume;
    TRIANGULATED_SURFACE<T>* ground_surface;
    ARRAY<MATRIX<T,3> > Dm_inverse;
    ARRAY<T> tet_rest_volume;
    
public:
    SIMULATION_LAYOUT(const STREAM_TYPE stream_type_input);
    void Initialize();
    void Add_Elastic_Forces(ARRAY_VIEW<const TV> X,ARRAY_VIEW<TV> force) const;
    void Add_Damping_Forces(ARRAY_VIEW<const TV> X,ARRAY_VIEW<const TV> V,
        ARRAY_VIEW<TV> force) const;
    void Add_Force_Differentials(ARRAY_VIEW<const TV> X,ARRAY_VIEW<const TV> DX,
	    ARRAY_VIEW<TV> force) const;
    void Add_External_Forces(ARRAY_VIEW<TV> force) const;
    T Maximum_Dt() const;
    void Write_Output(const int frame) const;
    void Set_Kinematic_Positions(const T time,ARRAY_VIEW<TV> X) const;
    void Set_Kinematic_Velocities(const T time,ARRAY_VIEW<TV> V) const;
    void Clear_Values_Of_Kinematic_Particles(ARRAY_VIEW<TV> array) const;

private:
    static ROTATED_STRESS_DERIVATIVE<T,3> Rotated_Stress_Derivative(const DIAGONAL_MATRIX<T,3>& Sigma,const T mu,const T lambda)
    {
        T J=Sigma.Determinant();

        T alpha11=mu;
        T alpha12=mu;
        T alpha13=mu;
        T alpha22=mu;
        T alpha23=mu;
        T alpha33=mu;

        T beta_const=(mu-lambda*log(J));
        T beta11=beta_const/sqr(Sigma(1,1));
        T beta12=beta_const/(Sigma(1,1)*Sigma(2,2));
        T beta13=beta_const/(Sigma(1,1)*Sigma(3,3));
        T beta22=beta_const/sqr(Sigma(2,2));
        T beta23=beta_const/(Sigma(2,2)*Sigma(3,3));
        T beta33=beta_const/sqr(Sigma(3,3));

        T gamma11=lambda/sqr(Sigma(1,1));
        T gamma12=lambda/(Sigma(1,1)*Sigma(2,2));
        T gamma13=lambda/(Sigma(1,1)*Sigma(3,3));
        T gamma22=lambda/sqr(Sigma(2,2));
        T gamma23=lambda/(Sigma(2,2)*Sigma(3,3));
        T gamma33=lambda/sqr(Sigma(3,3));
    
        ROTATED_STRESS_DERIVATIVE<T,3> rotated_dPdF;

        rotated_dPdF.a1111=alpha11+beta11+gamma11;
        rotated_dPdF.a1122=gamma12;
        rotated_dPdF.a1133=gamma13;
        rotated_dPdF.a2222=alpha22+beta22+gamma22;
        rotated_dPdF.a2233=gamma23;
        rotated_dPdF.a3333=alpha33+beta33+gamma33;

        rotated_dPdF.a1212=alpha12;
        rotated_dPdF.a1221=beta12;
        rotated_dPdF.a1313=alpha13;
        rotated_dPdF.a1331=beta13;
        rotated_dPdF.a2323=alpha23;
        rotated_dPdF.a2332=beta23;

        return rotated_dPdF;
    }
};

}
