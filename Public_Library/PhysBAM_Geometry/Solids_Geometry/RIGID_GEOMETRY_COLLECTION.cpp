//#####################################################################
// Copyright 2006-2009, Craig Schroeder, Tamar Shinar, Eftychios Sifakis.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class RIGID_GEOMETRY_COLLECTION
//#####################################################################
#include <PhysBAM_Geometry/Collisions/COLLISION_GEOMETRY_COLLECTION.h>
#include <PhysBAM_Geometry/Implicit_Objects/IMPLICIT_OBJECT.h>
#include <PhysBAM_Geometry/Implicit_Objects/IMPLICIT_OBJECT_TRANSFORMED.h>
#include <PhysBAM_Geometry/Solids_Geometry/RIGID_GEOMETRY.h>
#include <PhysBAM_Geometry/Solids_Geometry/RIGID_GEOMETRY_COLLECTION.h>
#include <PhysBAM_Geometry/Topology_Based_Geometry/STRUCTURE.h>
#include <PhysBAM_Geometry/Topology_Based_Geometry/STRUCTURE_LIST.h>
#ifndef COMPILE_WITHOUT_READ_WRITE_SUPPORT
#include <PhysBAM_Tools/Read_Write/Utilities/FILE_UTILITIES.h>
#endif
#ifndef COMPILE_WITHOUT_DYADIC_SUPPORT
#include <PhysBAM_Geometry/Implicit_Objects_Dyadic/DYADIC_IMPLICIT_OBJECT.h>
#endif
#ifndef COMPILE_WITHOUT_READ_WRITE_SUPPORT
#include <PhysBAM_Geometry/Read_Write/Geometry/READ_WRITE_STRUCTURE.h>
#endif
namespace PhysBAM{
template<class TV>
struct ALLOCATE_GEOMETRY_HELPER:public ALLOCATE_HELPER<TV>
{
    RIGID_GEOMETRY_COLLECTION<TV>& collection;
    ALLOCATE_GEOMETRY_HELPER(RIGID_GEOMETRY_COLLECTION<TV>& collection_input): collection(collection_input) {}
    RIGID_GEOMETRY<TV>* Create(int index=0) PHYSBAM_OVERRIDE {return new RIGID_GEOMETRY<TV>(collection,index);}
    virtual ~ALLOCATE_GEOMETRY_HELPER(){}
};
//#####################################################################
// Constructor
//#####################################################################
template<class TV> RIGID_GEOMETRY_COLLECTION<TV>::
RIGID_GEOMETRY_COLLECTION(RIGID_GEOMETRY_PARTICLES<TV>& particles_input,RIGID_GEOMETRY_EXAMPLE_VELOCITIES<TV>* rigid_geometry_example_velocities_input,
    COLLISION_GEOMETRY_COLLECTION<TV>* collision_body_list_input,ALLOCATE_HELPER<TV>* allocate_helper_input)
    :particles(particles_input),collision_body_list(collision_body_list_input),structure_list(*new STRUCTURE_LIST<TV,int>),always_create_structure(false),structure_hash(*new HASHTABLE<std::string,int>),frame_list_key(0),
    frame_list_active(0),check_stale(false),last_read_key(-1),last_read_active(-1),allocate_helper(allocate_helper_input),rigid_geometry_example_velocities(rigid_geometry_example_velocities_input),owns_particles(false)
{
    if(!allocate_helper) allocate_helper=new ALLOCATE_GEOMETRY_HELPER<TV>(*this);
    if(!collision_body_list) collision_body_list=new COLLISION_GEOMETRY_COLLECTION<TV>();
}
//#####################################################################
// Constructor
//#####################################################################
template<class TV> RIGID_GEOMETRY_COLLECTION<TV>::
RIGID_GEOMETRY_COLLECTION(RIGID_GEOMETRY_EXAMPLE_VELOCITIES<TV>* rigid_geometry_example_velocities_input,COLLISION_GEOMETRY_COLLECTION<TV>* collision_body_list_input,ALLOCATE_HELPER<TV>* allocate_helper_input)
    :particles(*new RIGID_GEOMETRY_PARTICLES<TV>()),collision_body_list(collision_body_list_input),structure_list(*new STRUCTURE_LIST<TV,int>),always_create_structure(false),structure_hash(*new HASHTABLE<std::string,int>),
    frame_list_key(0),frame_list_active(0),check_stale(false),last_read_key(-1),last_read_active(-1),allocate_helper(allocate_helper_input),rigid_geometry_example_velocities(rigid_geometry_example_velocities_input),owns_particles(true)
{
    if(!allocate_helper) allocate_helper=new ALLOCATE_GEOMETRY_HELPER<TV>(*this);
    if(!collision_body_list) collision_body_list=new COLLISION_GEOMETRY_COLLECTION<TV>();
}
//#####################################################################
// Destructor
//#####################################################################
template<class TV> RIGID_GEOMETRY_COLLECTION<TV>::
~RIGID_GEOMETRY_COLLECTION()
{
    if(owns_particles) delete &particles;
    delete &structure_hash;delete allocate_helper;delete &structure_list;
}
//#####################################################################
// Function Exists
//#####################################################################
template<class TV> bool RIGID_GEOMETRY_COLLECTION<TV>::
Exists(const int particle) const
{
    return particle>0 && particle<=particles.array_collection->Size() && particles.rigid_geometry(particle);
}
//#####################################################################
// Function Exists
//#####################################################################
template<class TV> bool RIGID_GEOMETRY_COLLECTION<TV>::
Is_Active(const int particle) const
{
    return Exists(particle) && particles.rigid_geometry(particle)->particle_index>0;
}
//#####################################################################
// Function Update_Kinematic_Particles
//#####################################################################
template<class TV> void RIGID_GEOMETRY_COLLECTION<TV>::
Update_Kinematic_Particles()
{
    static_rigid_geometry.Remove_All();kinematic_rigid_geometry.Remove_All();
    for(int p=1;p<=particles.array_collection->Size();p++) if(Is_Active(p)){RIGID_GEOMETRY<TV>& rigid_geometry=Rigid_Geometry(p);
        if(rigid_geometry.is_static) static_rigid_geometry.Append(p); else kinematic_rigid_geometry.Append(p);}
}
//#####################################################################
// Function Add_Rigid_Geometry
//#####################################################################
template<class TV> int RIGID_GEOMETRY_COLLECTION<TV>::
Add_Rigid_Geometry(STREAM_TYPE stream_type,const std::string& basename,const T scaling_factor,const bool read_simplicial_boundary,const bool read_implicit_object,const bool read_simplicial_interior,const bool read_rgd_file)
{
    RIGID_GEOMETRY<TV>* rigid_geometry=new RIGID_GEOMETRY<TV>(*this,true);
    return Add_Rigid_Geometry(rigid_geometry,stream_type,basename,scaling_factor,read_simplicial_boundary,read_implicit_object,read_simplicial_interior,read_rgd_file);
}
//#####################################################################
// Function Add_Rigid_Geometry
//#####################################################################
template<class TV> int RIGID_GEOMETRY_COLLECTION<TV>::
Add_Rigid_Geometry(RIGID_GEOMETRY<TV>* rigid_geometry,STREAM_TYPE stream_type,const std::string& basename,const T scaling_factor,
    const bool read_simplicial_boundary,const bool read_implicit_object,const bool read_simplicial_interior,const bool read_rgd_file)
{
    int id=rigid_geometry->particle_index;

    // structures
    ARRAY<int> structure_ids;
    TV structure_center=particles.X(id);
    if(TV::dimension==2){
        if(read_simplicial_boundary && !Find_Or_Read_Structure(stream_type,structure_ids,basename+".curve2d",scaling_factor,structure_center))
            LOG::cout<<"Note: No curve2d file for "<<basename<<std::endl;
        if(read_implicit_object && !Find_Or_Read_Structure(stream_type,structure_ids,basename+".phi2d",scaling_factor,structure_center))
            LOG::cout<<"Note: No phi2d file for "<<basename<<std::endl;
        if(read_simplicial_interior && !Find_Or_Read_Structure(stream_type,structure_ids,basename+".tri2d",scaling_factor,structure_center))
            LOG::cout<<"Note: No tri2d file for "<<basename<<std::endl;}
    else{
        if(read_simplicial_boundary && !Find_Or_Read_Structure(stream_type,structure_ids,basename+".tri",scaling_factor,structure_center))
            LOG::cout<<"Note: No tri file for "<<basename<<std::endl;
        if(read_implicit_object && !Find_Or_Read_Structure(stream_type,structure_ids,basename+".phi",scaling_factor,structure_center) && 
            !Find_Or_Read_Structure(stream_type,structure_ids,basename+".oct",scaling_factor,structure_center))
            LOG::cout<<"Note: No phi or oct file for "<<basename<<std::endl;
        if(read_simplicial_interior && !Find_Or_Read_Structure(stream_type,structure_ids,basename+".tet",scaling_factor,structure_center))
            LOG::cout<<"Note: No tet file for "<<basename<<std::endl;}
    assert(structure_ids.m<=3);
    particles.structure_ids(id)=VECTOR<int,3>();
    for(int i=1;i<=structure_ids.m;i++){
        if(structure_ids(i)){
            particles.structure_ids(id)(i)=structure_ids(i);
            Rigid_Geometry(id).Add_Structure(*structure_list.Element(structure_ids(i)));}}

    return id;
}
//#####################################################################
// Function Register_Analytic_Replacement_Structure
//#####################################################################
template<class TV> bool RIGID_GEOMETRY_COLLECTION<TV>::
Register_Analytic_Replacement_Structure(const std::string& filename,const T scaling_factor,STRUCTURE<TV>* structure)
{
    std::string hashname=STRING_UTILITIES::string_sprintf("%s@%.6f",filename.c_str(),scaling_factor); // mangle hash name
    if(structure_hash.Contains(hashname)) return false;
    int id=structure?structure_list.Add_Element(structure):0;
    structure_hash.Insert(hashname,id);
    return true;
}
template<class T> void
Wrap_Structure_Helper(STRUCTURE<VECTOR<T,1> >*& structure,const VECTOR<T,1>& center)
{}
template<class TV> void
Wrap_Structure_Helper(STRUCTURE<TV>*& structure,const TV& center)
{   // TODO(jontg): This shouldn't be necessary
#ifndef COMPILE_WITHOUT_DYADIC_SUPPORT
    if(DYADIC_IMPLICIT_OBJECT<TV>* ptr=dynamic_cast<DYADIC_IMPLICIT_OBJECT<TV>*>(structure))
        structure=new IMPLICIT_OBJECT_TRANSFORMED<TV,typename TV::SCALAR>(ptr,true,center,1);    
#endif
}
#ifndef COMPILE_WITHOUT_READ_WRITE_SUPPORT
//#####################################################################
// Function Find_Or_Read_Structure
//#####################################################################
template<class TV> bool RIGID_GEOMETRY_COLLECTION<TV>::
Find_Or_Read_Structure(const STREAM_TYPE stream_type,ARRAY<int>& structure_ids,const std::string& filename,const T scaling_factor,const TV& center)
{
    int id;
    if(!FILE_UTILITIES::File_Exists(filename)) return false;
    std::string hashname=STRING_UTILITIES::string_sprintf("%s@%.6f",filename.c_str(),scaling_factor); // mangle hash name
    if(!always_create_structure&&structure_hash.Get(hashname,id)){ // already read in
        if(!structure_list.Is_Active(id)) PHYSBAM_FATAL_ERROR();} // // only works if the referenced geometry is still in memory
    else{ // read in for the first time
        STRUCTURE<TV>* structure=0;
        if(!stream_type.use_doubles)
            structure=Read_Write<STRUCTURE<TV>,float>::Create_From_File(filename);
#ifndef COMPILE_WITHOUT_DOUBLE_SUPPORT
        else
            structure=Read_Write<STRUCTURE<TV>,double>::Create_From_File(filename);
#endif
        if(scaling_factor!=1){
            Wrap_Structure_Helper(structure,center);
            structure->Rescale(scaling_factor);}
        id=structure_list.Add_Element(structure);
        if(!always_create_structure) structure_hash.Insert(hashname,id);}
    structure_ids.Append(id);
    return true;
}
#endif
//#####################################################################
// Function Destroy_Unreferenced_Geometry
//#####################################################################
template<class TV> void RIGID_GEOMETRY_COLLECTION<TV>::
Destroy_Unreferenced_Geometry() 
{
    ARRAY<bool> referenced(structure_list.Number_Of_Active_Elements());
    for(int i=1;i<=particles.array_collection->Size();i++) for(int j=1;j<=particles.structure_ids(i).m;j++) if(particles.structure_ids(i)(j))
        referenced(structure_list.Element_Index(particles.structure_ids(i)(j)))=true;
    for(int i=structure_list.Number_Of_Active_Elements();i>=1;i--) if(!referenced(i)) structure_list.Remove_Element(structure_list.Active_Element_Id(i));
}
//#####################################################################
// Function Destroy_Unreferenced_Geometry
//#####################################################################
template<class TV> RIGID_GEOMETRY<TV>* RIGID_GEOMETRY_COLLECTION<TV>::
New_Body(int index)
{
    return allocate_helper->Create(index);
}
//#####################################################################
template class RIGID_GEOMETRY_COLLECTION<VECTOR<float,1> >;
template class RIGID_GEOMETRY_COLLECTION<VECTOR<float,2> >;
template class RIGID_GEOMETRY_COLLECTION<VECTOR<float,3> >;
#ifndef COMPILE_WITHOUT_DOUBLE_SUPPORT
template class RIGID_GEOMETRY_COLLECTION<VECTOR<double,1> >;
template class RIGID_GEOMETRY_COLLECTION<VECTOR<double,2> >;
template class RIGID_GEOMETRY_COLLECTION<VECTOR<double,3> >;
#endif
}
