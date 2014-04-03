//#####################################################################
// Copyright 2002-2006, Robert Bridson, Ronald Fedkiw, Geoffrey Irving, Sergey Koltakov, Andrew Selle, Michael Turitzin, Jiayi Chong.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class RENDERING_TRIANGLULATED_SURFACE
//#####################################################################
#ifndef __RENDERING_TRIANGULATED_SURFACE__
#define __RENDERING_TRIANGULATED_SURFACE__

#include <PhysBAM_Tools/Arrays/PROJECTED_ARRAY.h>
#include <PhysBAM_Tools/Grids_Uniform_Interpolation/LINEAR_INTERPOLATION_UNIFORM.h>
#include <PhysBAM_Tools/Images/BMP_FILE.h>
#include <PhysBAM_Tools/Interpolation/LINEAR_INTERPOLATION.h>
#include <PhysBAM_Tools/Log/LOG.h>
#include <PhysBAM_Tools/Log/PROGRESS_INDICATOR.h>
#include <PhysBAM_Geometry/Basic_Geometry_Intersections/RAY_TRIANGLE_3D_INTERSECTION.h>
#include <PhysBAM_Geometry/Topology_Based_Geometry/TRIANGULATED_SURFACE.h>
#include <PhysBAM_Geometry/Topology_Based_Geometry_Intersections/RAY_TRIANGULATED_SURFACE_INTERSECTION.h>
#include <PhysBAM_Rendering/PhysBAM_Ray_Tracing/Rendering_Objects/RENDERING_OBJECT.h>
#include <PhysBAM_Rendering/PhysBAM_Ray_Tracing/Rendering_Shaders/RENDERING_BSSRDF_SHADER.h>
#include <PhysBAM_Rendering/PhysBAM_Ray_Tracing/Rendering_Shaders/SUBSURFACE_SCATTERING_SAMPLED_IRRADIANCE.h>
namespace PhysBAM{

template<class T> class RENDERING_BSSRDF_SHADER;

template<class T>
class RENDERING_TRIANGULATED_SURFACE:public RENDERING_OBJECT<T>
{
    typedef VECTOR<T,3> TV;
public:
    using RENDERING_OBJECT<T>::small_number;using RENDERING_OBJECT<T>::transform;using RENDERING_OBJECT<T>::material_shader;
#ifndef COMPILE_WITHOUT_DYADIC_SUPPORT
    using RENDERING_OBJECT<T>::bssrdf_tree;using RENDERING_OBJECT<T>::bssrdf_shader;
#endif
    using RENDERING_OBJECT<T>::flip_normal;using RENDERING_OBJECT<T>::World_Space_Bounding_Box;using RENDERING_OBJECT<T>::name;using RENDERING_OBJECT<T>::two_sided;
    using RENDERING_OBJECT<T>::Inside;using RENDERING_OBJECT<T>::Intersection; // silence -Woverloaded-virtual
    
    TRIANGULATED_SURFACE<T>& triangulated_surface;
    mutable TRIANGULATED_SURFACE<T>* base_triangulated_surface;
    mutable ARRAY<TRIANGLE_3D<T> > world_space_triangles;
    bool closed_volume;
    ARRAY<TV,VECTOR<int,2> > bump_map_pixels;
    LINEAR_INTERPOLATION_UNIFORM<GRID<VECTOR<T,2> >,TV> interpolation;
    GRID<VECTOR<T,2> > grid;
    ARRAY<VECTOR<T,2> >* texture_coordinates;
    ARRAY<VECTOR<int,3> >* triangle_texture_coordinates;
    ARRAY<TV>* tangent_vectors;
    std::string sample_locations_file;
    bool add_triangles_to_acceleration_structure;

    RENDERING_TRIANGULATED_SURFACE(TRIANGULATED_SURFACE<T>& triangulated_surface_input,const int triangles_per_hierarchy_group=0)
        :triangulated_surface(triangulated_surface_input),base_triangulated_surface(&triangulated_surface),closed_volume(true),texture_coordinates(0),triangle_texture_coordinates(0),tangent_vectors(0),add_triangles_to_acceleration_structure(true)
    {
        triangulated_surface.Update_Bounding_Box();if(!triangulated_surface.hierarchy)triangulated_surface.Initialize_Hierarchy(true,triangles_per_hierarchy_group);
        if(triangulated_surface.use_vertex_normals && !triangulated_surface.face_vertex_normals)triangulated_surface.Update_Vertex_Normals();
    }

    virtual ~RENDERING_TRIANGULATED_SURFACE()
    {
#ifndef COMPILE_WITHOUT_DYADIC_SUPPORT
    delete bssrdf_tree;
#endif
    delete texture_coordinates;delete triangle_texture_coordinates;delete tangent_vectors;}

    bool Intersection(RAY<TV>& ray) const PHYSBAM_OVERRIDE
    {RAY<TV> object_space_ray=Object_Space_Ray(ray);
    if(INTERSECTION::Intersects(object_space_ray,*base_triangulated_surface,small_number)){
        ray.semi_infinite=false;ray.t_max=object_space_ray.t_max;ray.aggregate_id=object_space_ray.aggregate_id;return true;}
    else return false;}

    TV Normal(const TV& location,const int aggregate=0) const PHYSBAM_OVERRIDE
    {return World_Space_Vector(triangulated_surface.Normal(Object_Space_Point(location),aggregate));}

    bool Inside(const TV& location) const PHYSBAM_OVERRIDE
    {return closed_volume && triangulated_surface.Inside(Object_Space_Point(location),small_number);}

    bool Outside(const TV& location) const PHYSBAM_OVERRIDE
    {return !closed_volume || triangulated_surface.Outside(Object_Space_Point(location),small_number);}

    bool Boundary(const TV& location) const PHYSBAM_OVERRIDE
    {return triangulated_surface.Boundary(Object_Space_Point(location),small_number);}

    TV Surface(const TV& location) const PHYSBAM_OVERRIDE
    {return World_Space_Point(triangulated_surface.Surface(Object_Space_Point(location)));}

    bool Has_Bounding_Box() const  PHYSBAM_OVERRIDE
    {return true;}
    
    RANGE<TV> Object_Space_Bounding_Box() const PHYSBAM_OVERRIDE
    {if(!triangulated_surface.bounding_box) triangulated_surface.Update_Bounding_Box();
    return *triangulated_surface.bounding_box;}

    bool Closed_Volume() const PHYSBAM_OVERRIDE
    {return closed_volume;}

    bool Close_To_Open_Surface(const TV& location,const T threshold_distance) const PHYSBAM_OVERRIDE
    {assert(!Closed_Volume());int closest_triangle;T distance;
    triangulated_surface.Surface(Object_Space_Point(location),threshold_distance,small_number,&closest_triangle,&distance);
    return distance<=threshold_distance;}

    bool Intersection(RAY<TV>& ray,const int aggregate) const PHYSBAM_OVERRIDE
    {if(!add_triangles_to_acceleration_structure) return Intersection(ray);
    if(INTERSECTION::Intersects(ray,(world_space_triangles)(aggregate),small_number)){ray.aggregate_id=aggregate;return true;}else return false;}
    
    void Get_Aggregate_World_Space_Bounding_Boxes(ARRAY<RENDERING_OBJECT_ACCELERATION_PRIMITIVE<T> >& primitives) const PHYSBAM_OVERRIDE
    {if(add_triangles_to_acceleration_structure){
        world_space_triangles.Remove_All();
        for(int i=1;i<=triangulated_surface.mesh.elements.m;i++){
            int node1,node2,node3;triangulated_surface.mesh.elements(i).Get(node1,node2,node3);
            TRIANGLE_3D<T> world_space_triangle(transform.Homogeneous_Times(triangulated_surface.particles.X(node1)),transform.Homogeneous_Times(triangulated_surface.particles.X(node2)),transform.Homogeneous_Times(triangulated_surface.particles.X(node3)));
            world_space_triangles.Append(world_space_triangle);
            primitives.Append(RENDERING_OBJECT_ACCELERATION_PRIMITIVE<T>(world_space_triangle.Bounding_Box(),this,i));}}
    else{
        triangulated_surface.Update_Bounding_Box();
        primitives.Append(RENDERING_OBJECT_ACCELERATION_PRIMITIVE<T>(World_Space_Bounding_Box(),this,1));}}

    TRIANGULATED_SURFACE<T>* Generate_Triangles() const PHYSBAM_OVERRIDE
    {return triangulated_surface.Create_Compact_Copy();}
    
    void Get_Texture_Coordinates(const TV& object_space_point,const int aggregate,T& s,T& t) const PHYSBAM_OVERRIDE
    {assert(texture_coordinates);
    int node1,node2,node3;triangulated_surface.mesh.elements(aggregate).Get(node1,node2,node3);
    int uv1,uv2,uv3;(*triangle_texture_coordinates)(aggregate).Get(uv1,uv2,uv3);
    assert(triangulated_surface.mesh.elements.m==triangle_texture_coordinates->m);
    TV weights=TRIANGLE_3D<T>::Barycentric_Coordinates(object_space_point,triangulated_surface.particles.X(node1),triangulated_surface.particles.X(node2),triangulated_surface.particles.X(node3));
    VECTOR<T,2> coordinates=(*texture_coordinates)(uv1)*weights.x+(*texture_coordinates)(uv2)*weights.y+(*texture_coordinates)(uv3)*weights.z;
    s=coordinates.x;t=coordinates.y;}

    void Get_Object_Space_Tangent_And_Bitangent(const TV& object_space_point,const TV& object_space_normal,const int aggregate,TV& object_tangent,TV& object_bitangent) const PHYSBAM_OVERRIDE
    {assert(tangent_vectors);
    int node1,node2,node3;triangulated_surface.mesh.elements(aggregate).Get(node1,node2,node3);
    TV weights=TRIANGLE_3D<T>::Barycentric_Coordinates(object_space_point,triangulated_surface.particles.X(node1),triangulated_surface.particles.X(node2),triangulated_surface.particles.X(node3));
    TV interpolated_tangent=(*tangent_vectors)(node1)*weights.x+(*tangent_vectors)(node2)*weights.y+(*tangent_vectors)(node3)*weights.z;
    object_tangent=(interpolated_tangent-object_space_normal*TV::Dot_Product(object_space_normal,interpolated_tangent)).Normalized();
    object_bitangent=-TV::Cross_Product(object_space_normal,object_tangent);}

    void Get_World_Space_Tangent_And_Bitangent(const TV& world_space_point,const TV& world_space_normal,const int aggregate,TV& world_tangent,TV& world_bitangent) const PHYSBAM_OVERRIDE
    {TV object_tangent,object_bitangent;
    Get_Object_Space_Tangent_And_Bitangent(Object_Space_Point(world_space_point),Object_Space_Vector(world_space_normal),aggregate,object_tangent,object_bitangent);
    world_tangent=World_Space_Vector(object_tangent);world_bitangent=World_Space_Vector(object_bitangent);}

    void Compute_Per_Vertex_Tangent_Vectors()
    {assert(texture_coordinates&&triangulated_surface.vertex_normals);assert(!tangent_vectors);
    tangent_vectors=new ARRAY<TV>;
    tangent_vectors->Resize(triangulated_surface.mesh.number_nodes);
    for(int i=1;i<=triangulated_surface.mesh.elements.m;++i){
        int index1,index2,index3;triangulated_surface.mesh.elements(i).Get(index1,index2,index3);
        int uv1,uv2,uv3;(*triangle_texture_coordinates)(i).Get(uv1,uv2,uv3);
        TV p1=triangulated_surface.particles.X(index1),p2=triangulated_surface.particles.X(index2),p3=triangulated_surface.particles.X(index3);
        VECTOR<T,2> st1=(*texture_coordinates)(uv1),st2=(*texture_coordinates)(uv2),st3=(*texture_coordinates)(uv3);
        T denominator=((st2.y-st1.y)*(st3.x-st1.x)-(st3.y-st1.y)*(st2.x-st1.x));
        TV dpds=((st2.y-st1.y)*(p3-p1)-(st3.y-st1.y)*(p2-p1))/denominator;
        (*tangent_vectors)(index1)+=dpds;(*tangent_vectors)(index2)+=dpds;(*tangent_vectors)(index3)+=dpds;}
    for(int i=1;i<=tangent_vectors->m;++i){
        TV cur_normal=(*triangulated_surface.vertex_normals)(i);
        // Gram-Schmidt orthogonalize and normalize each tangent vector
        (*tangent_vectors)(i)=((*tangent_vectors)(i)-cur_normal*TV::Dot_Product(cur_normal,(*tangent_vectors)(i))).Normalized();}}

    template<class RW> 
    void Read_Texture_Coordinates(const std::string& filename)
    {assert(!texture_coordinates);texture_coordinates=new ARRAY<VECTOR<T,2> >;triangle_texture_coordinates=new ARRAY<VECTOR<int,3> >;
    int backward_compatible;FILE_UTILITIES::Read_From_File<RW>(filename,*texture_coordinates,backward_compatible,*triangle_texture_coordinates);
    Compute_Per_Vertex_Tangent_Vectors();}

    void Rescale_Texture_Coordinates(T scale)
    {assert(texture_coordinates);
    for(int i=1;i<=texture_coordinates->m;i++)
    {T x_rescaled=(*texture_coordinates)(i).x*scale,y_rescaled=(*texture_coordinates)(i).y*scale;
    (*texture_coordinates)(i).x=x_rescaled-floor(x_rescaled);
    (*texture_coordinates)(i).y=y_rescaled-floor(y_rescaled);}}

#ifndef COMPILE_WITHOUT_DYADIC_SUPPORT
    void Generate_BSSRDF_Tree(RENDER_WORLD<T>& world)  PHYSBAM_OVERRIDE
    {LOG::SCOPE scope("BSSRDF Irradiance Precomputation","BSSRDF Irradiance Precomputation %s",name.c_str());
    bssrdf_tree=new SUBSURFACE_SCATTERING_SAMPLED_IRRADIANCE<T>(bssrdf_shader->samples_per_octree_cell,bssrdf_shader->error_criterion);
    if(!bssrdf_shader->use_irradiance_cache_file || !FILE_UTILITIES::File_Exists(bssrdf_shader->irradiance_cache_filename)){
        LOG::Time("Initialize irradiance samples");
        if(sample_locations_file.length())bssrdf_tree->template Initialize_Sample_Locations_From_File<float>(sample_locations_file,triangulated_surface);
        else bssrdf_tree->Initialize_Sample_Locations_From_Vertices(triangulated_surface);
        LOG::Time("Compute irradiance");
        PROGRESS_INDICATOR progress(bssrdf_tree->samples.m);
        int flip_normal_factor=(flip_normal==true)?-1:1;
        for(int i=1;i<=bssrdf_tree->samples.m;++i){
            progress.Progress();
            bssrdf_tree->samples(i).position=World_Space_Point(bssrdf_tree->samples(i).position);
            bssrdf_tree->samples(i).normal=World_Space_Vector(bssrdf_tree->samples(i).normal)*(T)(flip_normal_factor);
            bssrdf_tree->samples(i).transmitted_irradiance=bssrdf_shader->Calculate_Surface_Irradiance(bssrdf_tree->samples(i).position,bssrdf_tree->samples(i).normal,world,*this);
            bssrdf_tree->samples(i).transmitted_irradiance_product=bssrdf_tree->samples(i).transmitted_irradiance*bssrdf_tree->samples(i).area;}
        bssrdf_tree->bounding_box=RANGE<TV>::Bounding_Box(bssrdf_tree->samples.template Project<TV,&SUBSURFACE_SCATTERING_IRRADIANCE_SAMPLE<T>::position>());
        if(bssrdf_shader->use_irradiance_cache_file){
            LOG::Time("Storing Irradiance Cache File");
            FILE_UTILITIES::Write_To_File<T>(bssrdf_shader->irradiance_cache_filename,*bssrdf_tree);
            LOG::Stop_Time();}}
    else{LOG::Time(STRING_UTILITIES::string_sprintf("Reading Irradiance Cache File %s",bssrdf_shader->irradiance_cache_filename.c_str()));FILE_UTILITIES::Read_From_File<T>(bssrdf_shader->irradiance_cache_filename,*bssrdf_tree);}
    bssrdf_tree->Build_Octree();}
#endif

    void Initialize_Bump_Map(const std::string& filename)
    {IMAGE<T>::Read(filename,bump_map_pixels);grid.Initialize(bump_map_pixels.counts,RANGE<VECTOR<T,2> >::Unit_Box());} 
     
    void Do_Displacement_Map_Per_Vertex(const T perturb_factor,const T perturb_power)
    {LOG::cerr<<"Doing per vertex displacement..." << std::endl;
    // TODO: fix to average all vertex uv's
    assert(texture_coordinates->m==triangulated_surface.particles.array_collection->Size());
    for(int i=1;i<=triangulated_surface.mesh.number_nodes;++i){
        TV current_perturbation=interpolation.Clamped_To_Array_Cell(grid,bump_map_pixels,VECTOR<T,2>((*texture_coordinates)(i).x,(*texture_coordinates)(i).y));
        current_perturbation=(current_perturbation-TV(0.5,0.5,0.5))*(T)2;
        T current_perturbation_sum=current_perturbation.x+current_perturbation.y+current_perturbation.z;
        triangulated_surface.particles.X(i)+=perturb_factor*std::pow(current_perturbation_sum,perturb_power)*(*triangulated_surface.vertex_normals)(i);}
    triangulated_surface.Initialize_Hierarchy();triangulated_surface.Update_Vertex_Normals();}

//#####################################################################
};   
}
#endif
