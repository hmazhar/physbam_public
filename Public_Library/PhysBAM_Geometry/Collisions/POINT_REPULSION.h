//#####################################################################
// Copyright 2004, Jiayi Chong, Geoffrey Irving, Igor Neverov, Andy Selle, Michael Turitzin.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class POINT_REPULSION
//#####################################################################
#ifndef __POINT_REPULSION__
#define __POINT_REPULSION__

#include <PhysBAM_Tools/Grids_Uniform_Arrays/ARRAYS_ND.h>
#include <PhysBAM_Tools/Log/LOG.h>
#include <PhysBAM_Tools/Random_Numbers/RANDOM_NUMBERS.h>
#include <PhysBAM_Tools/Read_Write/Utilities/READ_WRITE_FUNCTIONS.h>
#include <PhysBAM_Geometry/Basic_Geometry_Intersections/RAY_TRIANGLE_3D_INTERSECTION.h>
#include <PhysBAM_Geometry/Topology_Based_Geometry/TRIANGULATED_SURFACE.h>
namespace PhysBAM{

template<class T>
struct POINT_REPULSION_DATA
{
public:
    VECTOR<T,3> position;
    int triangle;
    ARRAY<int> neighbors;

    POINT_REPULSION_DATA()
    {}

    bool operator==(const POINT_REPULSION_DATA& data) const
    {return position==data.position && triangle==data.triangle && neighbors==data.neighbors;}

    template<class RW>
    void Read(std::istream& input)
    {Read_Binary<RW>(input,position,triangle,neighbors);}

    template<class RW>
    void Write(std::ostream& output) const
    {Write_Binary<RW>(output,position,triangle,neighbors);}
};


template<class T>
class POINT_REPULSION
{
    typedef VECTOR<T,3> TV;
public:
    TRIANGULATED_SURFACE<T>& surface;
    ARRAY<T> radii_of_curvature;
    ARRAY<ARRAY<int> > points_in_triangle;
    ARRAY<ARRAY<int> > vicinity_triangles;
    ARRAY<POINT_REPULSION_DATA<T> > points;
    ARRAY<T> sample_pdf;
    T surface_area;
    T surface_linear_size;
    bool use_curvature_based_repulsion_radius_restriction;
    T standard_repulsion_radius;
    T overpopulation_fraction;
    T number_points_float;
    T attenuation_factor;
    int number_points,active_point;
    T average_repulsion_radius;
    T average_number_vicinity_triangles;
    T average_number_neighbors;
    T average_distance_to_neighbor_over_repulsion_radius;
    bool draw_points,draw_active_point; T point_size;
    mutable RANDOM_NUMBERS<T> random;

    POINT_REPULSION(TRIANGULATED_SURFACE<T>& surface_in,const int point_number)
        :surface(surface_in),surface_area(0),surface_linear_size(0),use_curvature_based_repulsion_radius_restriction(false),standard_repulsion_radius(0),
        overpopulation_fraction(3),number_points_float(0),attenuation_factor(1),number_points(10),active_point(0),average_repulsion_radius(0),average_number_vicinity_triangles(0),
        average_number_neighbors(0),average_distance_to_neighbor_over_repulsion_radius(0),point_size(3)
    {
        random.Set_Seed(2345);Initialize(point_number);
    }

    T Standard_Repulsion_Radius(const int number_of_points) const
    {return sqrt((T)3.6276*surface_area/((T)pi*number_of_points));}

    T Repulsion_Radius(const T radius_of_curvature) const
    {T radius=standard_repulsion_radius;if(use_curvature_based_repulsion_radius_restriction) radius=min(radius,(T).5*radius_of_curvature);return radius;}

    T Repulsion_Radius(const int radii_input) const
    {if(use_curvature_based_repulsion_radius_restriction && radii_of_curvature.m!=surface.mesh.elements.m)
        PHYSBAM_FATAL_ERROR("Need radii_of_curvature to evaluate repulsion radius");
    return Repulsion_Radius(radii_of_curvature(radii_input));}

    static T Radius_Of_Curvature_Estimate(const TV& position_1,const TV& normal_1,const TV& position_2,const TV& normal_2)
    {return (position_1-position_2).Magnitude()/TV::Angle_Between(normal_1,normal_2);}

    void Initialize_Radii_Of_Curvature()
    {if(surface.vertex_normals==0) PHYSBAM_FATAL_ERROR("Need vertex_normals to initialize radii_of_curvature");
    radii_of_curvature.Resize(surface.mesh.elements.m);
    ARRAY_VIEW<const TV> X(surface.particles.X),normals(*surface.vertex_normals);
    for(int t=1;t<=surface.mesh.elements.m;++t){
        int i,j,k;surface.mesh.elements(t).Get(i,j,k);
        TV position=(T)one_third*(X(i)+X(j)+X(k)), normal=(T)one_third*(normals(i)+normals(j)+normals(k));
        radii_of_curvature(t)=min(Radius_Of_Curvature_Estimate(position,normal,X(i),normals(i)),
            Radius_Of_Curvature_Estimate(position,normal,X(j),normals(j)),
            Radius_Of_Curvature_Estimate(position,normal,X(k),normals(k)));}}

    T Number_Points_From_Repulsion_Radii() const
    {T reference_number=0;for(int t=1;t<=surface.mesh.elements.m;++t) reference_number+=(T)3.6276*surface.Area(t)/((T)pi*sqr(Repulsion_Radius(t)));
    return overpopulation_fraction*reference_number;}

    void Set_Standard_Repulsion_Radius(const T radius){standard_repulsion_radius=radius;number_points=(int)(number_points_float=Number_Points_From_Repulsion_Radii());}
    void Set_Standard_Number_Points(const int number){Set_Standard_Repulsion_Radius(Standard_Repulsion_Radius(number));}

    void Scale_Repulsion_Radii(const T scale)
    {if(number_points_float/sqr(scale)<=1) return;
    standard_repulsion_radius*=scale;
    number_points_float/=sqr(scale);
    number_points=(int)number_points_float;}

    void Scale_Number_Points(const T scale)
    {if(number_points_float*scale<=1) return;
    number_points_float*=scale;standard_repulsion_radius/=sqrt(scale);number_points=(int)number_points_float;}

    void Set_Number_Points(int number)
    {number=max(1,number);Scale_Number_Points(number/number_points_float);}

    T Distance_Between_Triangles(const int triangle_index_1,const int triangle_index_2) const
    {ARRAY_VIEW<const TV> X(surface.particles.X);
    int i1,j1,k1,i2,j2,k2;surface.mesh.elements(triangle_index_1).Get(i1,j1,k1);surface.mesh.elements(triangle_index_2).Get(i2,j2,k2);
    T distance_1=min((X(i1)-X(i2)).Magnitude_Squared(),(X(i1)-X(j2)).Magnitude_Squared(),(X(i1)-X(k2)).Magnitude_Squared()),
        distance_2=min((X(j1)-X(i2)).Magnitude_Squared(),(X(j1)-X(j2)).Magnitude_Squared(),(X(j1)-X(k2)).Magnitude_Squared()),
        distance_3=min((X(k1)-X(i2)).Magnitude_Squared(),(X(k1)-X(j2)).Magnitude_Squared(),(X(k1)-X(k2)).Magnitude_Squared());
    return sqrt(min(distance_1,distance_2,distance_3));}

    bool Acceptable_Vicinity_Triangle(const int triangle_index_1,const T radius,const int triangle_index_2) const
    {return Distance_Between_Triangles(triangle_index_1,triangle_index_2)<radius &&
    TV::Dot_Product((*surface.triangle_list)(triangle_index_1).Normal(),(*surface.triangle_list)(triangle_index_2).Normal())>0;}

    bool Acceptable_Neighbor_Point(const TV& point,const T radius_squared,const TV& normal,const int index) const
    {TV reference_position=points(index).position;
    return (point-reference_position).Magnitude_Squared()<radius_squared && TV::Dot_Product(normal,surface.Normal(reference_position,points(index).triangle))>0;}

    void Initialize_Sampling()
    {T sum=0;
    if(use_curvature_based_repulsion_radius_restriction)
        for(int t=1;t<=surface.mesh.elements.m;++t){sum+=surface.Area(t)/sqr(Repulsion_Radius(t));sample_pdf(t)=sum;}
    else
        for(int t=1;t<=surface.mesh.elements.m;++t){sum+=surface.Area(t);sample_pdf(t)=sum;}
    sample_pdf/=sample_pdf.Last();}

    int Triangle(const T fraction) const
    {int lo=1,hi=sample_pdf.m;
    while(lo<hi){
        int mid=(lo+hi)/2;
        if(fraction<=sample_pdf(mid)) hi=mid;
        else lo=mid+1;}
    return hi;}

    void Update_Points_In_Triangle()
    {points_in_triangle.Resize(surface.mesh.elements.m);
    for(int t=1;t<=points_in_triangle.m;++t)points_in_triangle(t).Resize(0);
    for(int i=1;i<=points.m;++i)points_in_triangle(points(i).triangle).Append(i);}

    void Seed_Points()
    {points.Resize(number_points);
    ARRAY_VIEW<const TV> X(surface.particles.X);
    for(int q=1;q<=number_points;++q){
        int triangle=Triangle((q-(T).5)/number_points);
        T a=random.Get_Uniform_Number(0,1);
        T b=random.Get_Uniform_Number(0,1);
        if(a+b>1){a=1-a;b=1-b;}
        const VECTOR<int,3>& nodes=surface.mesh.elements(triangle);
        points(q).position=(1-a-b)*X(nodes[1])+a*X(nodes[2])+b*X(nodes[3]);points(q).triangle=triangle;}
    attenuation_factor=1;Update_Points_In_Triangle();Update_Neighbor_Points();Update_Stats();}

    int Get_Subdivide_Particles_Number(int level,int triangle_number){return (int)(0.5*(pow(3.0,level)-1))*triangle_number;}

    void Set_Subdivided_Points(const TV& point1,const TV& point2,const TV& point3,const int triangle_index,int& particle_index,int level)
    {if(level<1)return;points(particle_index).position=(point1+point2+point3)/(T)3.0;points(particle_index).triangle=triangle_index;++particle_index;
    Set_Subdivided_Points(point1,point2,points(particle_index).position,triangle_index,particle_index,level-1);
    Set_Subdivided_Points(point1,point3,points(particle_index).position,triangle_index,particle_index,level-1);
    Set_Subdivided_Points(point2,point3,points(particle_index).position,triangle_index,particle_index,level-1);}

    void Seed_Subdivided_Points(int level)
    {int number_of_points=Get_Subdivide_Particles_Number(level,surface.mesh.elements.m),particle_index=1;
    Set_Standard_Number_Points(number_of_points);Set_Number_Points(number_of_points);
    points.Resize(0);points.Resize(number_of_points);
    for(int triangle=1;triangle<=surface.mesh.elements.m;++triangle){
        int i,j,k;surface.mesh.elements(triangle).Get(i,j,k);
        Set_Subdivided_Points(surface.particles.X(i),surface.particles.X(j),surface.particles.X(k),triangle,particle_index,level);}
    Update_Points_In_Triangle();Update_Neighbor_Points();Update_Stats();}

    void Flip_Use_Curvature_Based_Repulsion_Radius_Restriction()
    {use_curvature_based_repulsion_radius_restriction^=true;
    Set_Standard_Repulsion_Radius(standard_repulsion_radius);
    Update_Neighbor_Points();Update_Stats();}

    bool Intersection_With_Triangles(RAY<TV>& ray,const ARRAY<int>& triangles,int& new_triangle)
    {for(int i=1;i<=triangles.m;++i) if(INTERSECTION::Intersects(ray,(*surface.triangle_list)(triangles(i)))){new_triangle=triangles(i);return true;}
    return false;}

    void Get_Points(ARRAY<TV>& position) const
    {position.Resize(points.m);for(int i=1;i<=position.m;++i)position(i)=points(i).position;}

    void Get_Triangles(ARRAY<int>& triangles) const
    {triangles.Resize(points.m);for(int i=1;i<=triangles.m;++i) triangles(i)=points(i).triangle;}

    void Set_Points(ARRAY_VIEW<const TV> X)
    {points.Resize(X.Size());
    for(int i=1;i<=X.Size();++i){points(i).position=X(i);points(i).triangle=0;points(i).neighbors.Resize(0);}}

    template<class RW>
    void Write(std::ostream& output) const
    {if(!surface.triangle_list) surface.Update_Triangle_List();
    ARRAY<TV> point_weights(points.m);ARRAY<int> triangle_correspondences(points.m);
    ARRAY<T> triangle_areas(surface.triangle_list->m);
    for(int i=1;i<=surface.triangle_list->m;++i){
        triangle_areas(i)=(*surface.triangle_list)(i).Area();}
    for(int i=1;i<=points.m;++i){
        triangle_correspondences(i)=points(i).triangle;
        point_weights(i)=(*surface.triangle_list)(triangle_correspondences(i)).Barycentric_Coordinates(points(i).position);}
    Write_Binary<RW>(output,point_weights,triangle_correspondences,triangle_areas);}

    template<class RW>
    static void Read(std::istream& input,ARRAY<TV>& point_weights,ARRAY<int>& triangle_correspondences,ARRAY<T>& triangle_areas)
    {Read_Binary<RW>(input,point_weights,triangle_correspondences,triangle_areas);}

//#####################################################################
    void Initialize(const int point_number);
    void Initialize_Vicinity_Triangles();
    void Update_Neighbor_Points();
    void Move_Points();
    void Update_Stats();
//#####################################################################
};
// global functions
template<class T> inline std::ostream&
operator<<(std::ostream& output,const POINT_REPULSION_DATA<T>& v)
{output<<"triangle="<<v.triangle<<" "<<v.position<<std::endl;return output;}

template<class T> inline std::istream& operator>>(std::istream& input,POINT_REPULSION_DATA<T>& v)
{input>>v.triangle>>" ">>v.position;return input;}


}
#endif











