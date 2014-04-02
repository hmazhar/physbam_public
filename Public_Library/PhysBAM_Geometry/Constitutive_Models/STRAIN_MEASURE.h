//#####################################################################
// Copyright 2002-2007, Zhaosheng Bao, Ronald Fedkiw, Geoffrey Irving, Neil Molino, Igor Neverov, Eftychios Sifakis, Joseph Teran.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class STRAIN_MEASURE
//#####################################################################
#ifndef __STRAIN_MEASURE__
#define __STRAIN_MEASURE__

#include <PhysBAM_Tools/Arrays/ARRAY.h>
#include <PhysBAM_Tools/Log/DEBUG_UTILITIES.h>
#include <PhysBAM_Tools/Matrices/MATRIX_2X2.h>
#include <PhysBAM_Tools/Matrices/MATRIX_3X2.h>
#include <PhysBAM_Tools/Matrices/MATRIX_3X3.h>
#include <PhysBAM_Tools/Matrices/MATRIX_BASE.h>
#include <PhysBAM_Tools/Matrices/UPPER_TRIANGULAR_MATRIX_2X2.h>
#include <PhysBAM_Tools/Matrices/UPPER_TRIANGULAR_MATRIX_3X3.h>
#include <PhysBAM_Tools/Utilities/NONCOPYABLE.h>
#include <PhysBAM_Tools/Vectors/VECTOR.h>
#include <PhysBAM_Geometry/Geometry_Particles/GEOMETRY_PARTICLES.h>
#include <PhysBAM_Geometry/Topology/TOPOLOGY_POLICY.h>
#include <PhysBAM_Geometry/Topology_Based_Geometry/TOPOLOGY_BASED_SIMPLEX_POLICY.h>
namespace PhysBAM{

template<class TV,int d>
class STRAIN_MEASURE:public NONCOPYABLE
{
    typedef typename TV::SCALAR T;
    typedef typename TOPOLOGY_BASED_SIMPLEX_POLICY<TV,d>::OBJECT T_MESH_OBJECT;
    typedef typename MESH_POLICY<d>::MESH T_MESH;
    typedef MATRIX<T,TV::m,d> T_MATRIX;
    template<int d2> struct UNUSABLE{};
public:
    T_MESH_OBJECT& mesh_object;
    T_MESH& mesh;
    GEOMETRY_PARTICLES<TV>& particles;
    ARRAY<UPPER_TRIANGULAR_MATRIX<T,d> > Dm_inverse;
private:
    ARRAY<UPPER_TRIANGULAR_MATRIX<T,d> >* Dm_inverse_save;
public:

    STRAIN_MEASURE(T_MESH_OBJECT& mesh_object);
    ~STRAIN_MEASURE();

    T_MATRIX F(const int simplex) const
    {return Ds(particles.X,simplex)*Dm_inverse(simplex);}

    T J(const int simplex) const
    {return Ds(particles.X,simplex).Parallelpiped_Measure()*Dm_inverse(simplex).Determinant();}

    T_MATRIX Velocity_Gradient(const int simplex) const
    {return Ds(particles.V,simplex)*Dm_inverse(simplex);}

    T Rest_Altitude(const int simplex) const
    {return Dm_inverse(simplex).Inverse().Simplex_Minimum_Altitude();}

    T Minimum_Rest_Altitude() const
    {T altitude=FLT_MAX;for(int t=1;t<=Dm_inverse.m;t++)altitude=min(altitude,Rest_Altitude(t));return altitude;}

    T_MATRIX Ds(ARRAY_VIEW<const TV> X,const typename IF<d==1,int,UNUSABLE<1> >::TYPE simplex) const
    {int i,j;mesh.elements(simplex).Get(i,j);
    return T_MATRIX(X(j)-X(i));}

    T_MATRIX Ds(ARRAY_VIEW<const TV> X,const typename IF<d==2,int,UNUSABLE<2> >::TYPE simplex) const
    {int i,j,k;mesh.elements(simplex).Get(i,j,k);
    return T_MATRIX(X(j)-X(i),X(k)-X(i));}

    T_MATRIX Ds(ARRAY_VIEW<const TV> X,const typename IF<d==3,int,UNUSABLE<3> >::TYPE simplex) const
    {int i,j,k,l;mesh.elements(simplex).Get(i,j,k,l);
    return T_MATRIX(X(j)-X(i),X(k)-X(i),X(l)-X(i));}

    template<class T_ARRAY>
    static T_MATRIX Ds(const T_ARRAY& X,const VECTOR<int,2>& nodes)
    {STATIC_ASSERT((IS_SAME<typename T_ARRAY::ELEMENT,TV>::value));
    int i,j;nodes.Get(i,j);return T_MATRIX(X(j)-X(i));}

    template<class T_ARRAY>
    static T_MATRIX Ds(const T_ARRAY& X,const VECTOR<int,3>& nodes)
    {STATIC_ASSERT((IS_SAME<typename T_ARRAY::ELEMENT,TV>::value));
    int i,j,k;nodes.Get(i,j,k);return T_MATRIX(X(j)-X(i),X(k)-X(i));}

    template<class T_ARRAY>
    static T_MATRIX Ds(const T_ARRAY& X,const VECTOR<int,4>& nodes)
    {STATIC_ASSERT((IS_SAME<typename T_ARRAY::ELEMENT,TV>::value));
    int i,j,k,l;nodes.Get(i,j,k,l);return T_MATRIX(X(j)-X(i),X(k)-X(i),X(l)-X(i));}

    void Distribute_Force(ARRAY_VIEW<TV> F,const int element,const T_MATRIX& forces) const
    {Distribute_Force(F,mesh.elements(element),forces);}

    static void Distribute_Force(ARRAY_VIEW<TV> F,const VECTOR<int,3>& nodes,const MATRIX<T,2>& forces)
    {int i,j,k;nodes.Get(i,j,k);
    F(i)-=TV(forces.x[0]+forces.x[2],forces.x[1]+forces.x[3]);
    F(j)+=TV(forces.x[0],forces.x[1]);
    F(k)+=TV(forces.x[2],forces.x[3]);}

    static void Distribute_Force(ARRAY_VIEW<TV> F,const VECTOR<int,3>& nodes,const MATRIX<T,3,2>& forces)
    {int i,j,k;nodes.Get(i,j,k);
    F(i)-=TV(forces.x[0]+forces.x[3],forces.x[1]+forces.x[4],forces.x[2]+forces.x[5]);
    F(j)+=TV(forces.x[0],forces.x[1],forces.x[2]);
    F(k)+=TV(forces.x[3],forces.x[4],forces.x[5]);}

    static void Distribute_Force(ARRAY_VIEW<TV> F,const VECTOR<int,4>& nodes,const MATRIX<T,3>& forces)
    {int i,j,k,l;nodes.Get(i,j,k,l);
    F(i)-=TV(forces.x[0]+forces.x[3]+forces.x[6],forces.x[1]+forces.x[4]+forces.x[7],forces.x[2]+forces.x[5]+forces.x[8]);
    F(j)+=TV(forces.x[0],forces.x[1],forces.x[2]);
    F(k)+=TV(forces.x[3],forces.x[4],forces.x[5]);
    F(l)+=TV(forces.x[6],forces.x[7],forces.x[8]);}

    static void Distribute_Impulse(ARRAY_VIEW<TV> V,const VECTOR<int,3>& nodes,const MATRIX<T,2>& impulse,ARRAY_VIEW<const T> one_over_mass)
    {int i,j,k;nodes.Get(i,j,k);
    V(i)-=one_over_mass(i)*TV(impulse.x[0]+impulse.x[2],impulse.x[1]+impulse.x[3]);
    V(j)+=one_over_mass(j)*TV(impulse.x[0],impulse.x[1]);
    V(k)+=one_over_mass(k)*TV(impulse.x[2],impulse.x[3]);}

    static void Distribute_Impulse(ARRAY_VIEW<TV> V,const VECTOR<int,3>& nodes,const MATRIX<T,3,2>& impulse,ARRAY_VIEW<const T> one_over_mass)
    {int i,j,k;nodes.Get(i,j,k);
    V(i)-=one_over_mass(i)*TV(impulse.x[0]+impulse.x[3],impulse.x[1]+impulse.x[4],impulse.x[2]+impulse.x[5]);
    V(j)+=one_over_mass(j)*TV(impulse.x[0],impulse.x[1],impulse.x[2]);
    V(k)+=one_over_mass(k)*TV(impulse.x[3],impulse.x[4],impulse.x[5]);}

    static void Distribute_Impulse(ARRAY_VIEW<TV> V,const VECTOR<int,4>& nodes,const MATRIX<T,3>& impulse,ARRAY_VIEW<const T> one_over_mass)
    {int i,j,k,l;nodes.Get(i,j,k,l);
    V(i)-=one_over_mass(i)*TV(impulse.x[0]+impulse.x[3]+impulse.x[6],impulse.x[1]+impulse.x[4]+impulse.x[7],impulse.x[2]+impulse.x[5]+impulse.x[8]);
    V(j)+=one_over_mass(j)*TV(impulse.x[0],impulse.x[1],impulse.x[2]);
    V(k)+=one_over_mass(k)*TV(impulse.x[3],impulse.x[4],impulse.x[5]);
    V(l)+=one_over_mass(l)*TV(impulse.x[6],impulse.x[7],impulse.x[8]);}

//#####################################################################
    void Initialize_Dm_Inverse(ARRAY_VIEW<const TV> X);
    void Initialize_Dm_Inverse_Save();
    void Copy_Dm_Inverse_Save_Into_Dm_Inverse(const ARRAY<int>& map);
    void Initialize_Rest_State_To_Equilateral(const T side_length);
    void Print_Altitude_Statistics();
//#####################################################################
};
}
#endif
