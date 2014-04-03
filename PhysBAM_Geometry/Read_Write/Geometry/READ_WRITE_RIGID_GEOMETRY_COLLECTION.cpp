//#####################################################################
// Copyright 2004-2009, Ron Fedkiw, Eran Guendelman, Tamar Shinar, Joseph Teran, Rachel Weinstein.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#include <PhysBAM_Tools/Arrays/INDIRECT_ARRAY.h>
#include <PhysBAM_Tools/Arrays_Computations/ARRAY_COPY.h>
#include <PhysBAM_Geometry/Read_Write/Geometry/READ_WRITE_RIGID_GEOMETRY_COLLECTION.h>
#include <PhysBAM_Geometry/Read_Write/Geometry/READ_WRITE_STRUCTURE_LIST.h>
#include <PhysBAM_Geometry/Solids_Geometry/RIGID_GEOMETRY.h>
#include <PhysBAM_Geometry/Solids_Geometry/RIGID_GEOMETRY_COLLECTION.h>
#include <climits>
using namespace PhysBAM;
//#####################################################################
// Function Read
//#####################################################################
template<class RW,class TV> void Read_Write<RIGID_GEOMETRY_COLLECTION<TV>,RW>::
Read(const STREAM_TYPE stream_type,const std::string& directory,const int frame,RIGID_GEOMETRY_COLLECTION<TV>& object,ARRAY<int>* needs_init,ARRAY<int>* needs_destroy)
{
    Read_Write<STRUCTURE_LIST<TV,int>,RW>::Read(directory,"rigid_body_structure_",frame,object.structure_list);
    ARRAY<RIGID_GEOMETRY<TV>*> bodies(object.particles.rigid_geometry);
    ARRAYS_COMPUTATIONS::Fill(object.particles.rigid_geometry,(RIGID_GEOMETRY<TV>*)0);
    FILE_UTILITIES::Read_From_File(stream_type,STRING_UTILITIES::string_sprintf("%s/%d/rigid_geometry_particles",directory.c_str(),frame),object.particles);
    while(object.particles.rigid_geometry.m<bodies.m) delete bodies.Pop();
    object.particles.rigid_geometry.Prefix(bodies.m)=bodies;

    ARRAY<int> needs_init_default;
    if(needs_init) needs_init->Remove_All();
    if(needs_destroy) needs_destroy->Remove_All();
    char version;int last_id=0;ARRAY<int> active_ids;int local_frame=frame;
    std::string active_list_name=STRING_UTILITIES::string_sprintf("%s/common/rigid_body_active_ids_list",directory.c_str(),frame);
    std::string active_name=STRING_UTILITIES::string_sprintf("%s/%d/rigid_body_active_ids",directory.c_str(),frame);
    if(FILE_UTILITIES::File_Exists(active_list_name)){
        if(!object.frame_list_active){object.frame_list_active=new ARRAY<int>;FILE_UTILITIES::Read_From_File(stream_type,active_list_name,*object.frame_list_active);}
        local_frame=(*object.frame_list_active)(object.frame_list_active->Binary_Search(frame));}
    if(object.last_read_active!=local_frame && FILE_UTILITIES::File_Exists(active_name)){
        FILE_UTILITIES::Read_From_File(stream_type,active_name,version,last_id,active_ids);
        object.last_read_active=local_frame;PHYSBAM_ASSERT(version==1);
        if(needs_destroy) for(int i=last_id+1;i<=object.particles.array_collection->Size();i++) if(!object.particles.rigid_geometry(i)) needs_destroy->Append(i);
        object.particles.Resize(last_id);}
    else for(int id(1);id<=object.particles.array_collection->Size();id++) if(object.Is_Active(id)){active_ids.Append(id);}
    if(object.particles.rigid_geometry.Subset(active_ids).Contains(0)){ // don't need to re-read these things if we will not be initializing any newly-active bodies
        std::string key_file_list=STRING_UTILITIES::string_sprintf("%s/common/rigid_body_key_list",directory.c_str());
        std::string key_file=STRING_UTILITIES::string_sprintf("%s/%d/rigid_body_key",directory.c_str(),frame);
        char version;
        if(FILE_UTILITIES::File_Exists(key_file_list)){
            if(!object.frame_list_key){object.frame_list_key=new ARRAY<int>;FILE_UTILITIES::Read_From_File(stream_type,key_file_list,*object.frame_list_key);}
            local_frame=(*object.frame_list_key)(object.frame_list_key->Binary_Search(frame));}
        if(object.last_read_key!=local_frame && FILE_UTILITIES::File_Exists(key_file)){
            FILE_UTILITIES::Read_From_File(stream_type,key_file,version,object.particles.structure_ids);
            object.last_read_active=local_frame;PHYSBAM_ASSERT(version==2 || version==3);}

        try{
            std::istream* input=FILE_UTILITIES::Safe_Open_Input(directory+"/common/rigid_body_names",false);
            int num;*input>>num;input->ignore(INT_MAX,'\n');
            object.rigid_body_names.Resize(num);
            for(int i(1);i<=object.rigid_body_names.Size();i++) std::getline(*input,object.rigid_body_names(i));
            delete input;}
        catch(FILESYSTEM_ERROR&){
            LOG::cerr<<"Did not find rigid body names."<<std::endl;
            object.rigid_body_names.Clean_Memory();}}

    ARRAY<bool> exists(last_id);INDIRECT_ARRAY<ARRAY<bool>,ARRAY<int>&> subset=exists.Subset(active_ids);ARRAYS_COMPUTATIONS::Fill(subset,true);
    for(int i=1;i<=exists.m;i++) if(!exists(i) && object.Is_Active(i)) object.Deactivate_Geometry(i);
    if(active_ids.m>0){
        for(int i=1;i<=active_ids.m;i++){int p=active_ids(i);
            // initialize new rigid body with given id, and initialize geometry
            if(!object.particles.rigid_geometry(p)){
                if(needs_init) needs_init->Append(p);
                RIGID_GEOMETRY<TV>* rigid_geometry=object.New_Body(p);
                if(p<=object.rigid_body_names.Size()) rigid_geometry->Set_Name(object.rigid_body_names(p));
                for(int s=1;s<=object.particles.structure_ids(p).m;s++)
                    if(object.particles.structure_ids(p)(s))
                        rigid_geometry->Add_Structure(*object.structure_list.Element(object.particles.structure_ids(p)(s)));}
            if(object.Is_Active(p) && object.particles.structure_ids(p)==VECTOR<int,3>()) object.Deactivate_Geometry(p);}}
}
//#####################################################################
// Function Write
//#####################################################################
template<class RW,class TV> void Read_Write<RIGID_GEOMETRY_COLLECTION<TV>,RW>::
Write(const STREAM_TYPE stream_type,const std::string& directory,const int frame,const RIGID_GEOMETRY_COLLECTION<TV>& object)
{
    Read_Write<STRUCTURE_LIST<TV,int>,RW>::Write(directory,"rigid_body_structure_",frame,object.structure_list);
    FILE_UTILITIES::Write_To_File(stream_type,STRING_UTILITIES::string_sprintf("%s/%d/rigid_geometry_particles",directory.c_str(),frame),object.particles);

    // update names
    object.rigid_body_names.Resize(object.particles.array_collection->Size());
    ARRAY<int> active_ids;
    for(int id(1);id<=object.particles.array_collection->Size();id++) if(object.Is_Active(id)){active_ids.Append(id);object.rigid_body_names(id)=object.Rigid_Geometry(id).name;}
    if(active_ids.m>0 && !(object.check_stale && object.is_stale_active)){
        if(object.check_stale){
            if(!object.frame_list_active) object.frame_list_active=new ARRAY<int>;object.frame_list_active->Append(frame);
            FILE_UTILITIES::Write_To_File(stream_type,STRING_UTILITIES::string_sprintf("%s/common/rigid_body_active_ids_list",directory.c_str()),*object.frame_list_active);
            object.is_stale_active=false;}
        FILE_UTILITIES::Write_To_File(stream_type,STRING_UTILITIES::string_sprintf("%s/%d/rigid_body_active_ids",directory.c_str(),frame),(char)1,object.particles.array_collection->Size(),active_ids);}
    std::ostream* output=FILE_UTILITIES::Safe_Open_Output(directory+"/common/rigid_body_names",false);
    *output<<object.rigid_body_names.Size()<<std::endl;
    for(int i(1);i<=object.rigid_body_names.Size();i++) *output<<object.rigid_body_names(i)<<std::endl;
    delete output;
    char version=3;
    if(object.particles.structure_ids.m>0 && !(object.check_stale && object.is_stale_key)){
        if(object.check_stale){
            if(!object.frame_list_key) object.frame_list_key=new ARRAY<int>;object.frame_list_key->Append(frame);
            FILE_UTILITIES::Write_To_File(stream_type,STRING_UTILITIES::string_sprintf("%s/common/rigid_body_key_list",directory.c_str()),*object.frame_list_key);
            object.is_stale_key=false;}
        FILE_UTILITIES::Write_To_File(stream_type,STRING_UTILITIES::string_sprintf("%s/%d/rigid_body_key",directory.c_str(),frame),version,object.particles.structure_ids);}
}
template class Read_Write<RIGID_GEOMETRY_COLLECTION<VECTOR<float,1> >,float>;
template class Read_Write<RIGID_GEOMETRY_COLLECTION<VECTOR<float,2> >,float>;
template class Read_Write<RIGID_GEOMETRY_COLLECTION<VECTOR<float,3> >,float>;
#ifndef COMPILE_WITHOUT_DOUBLE_SUPPORT
template class Read_Write<RIGID_GEOMETRY_COLLECTION<VECTOR<double,1> >,double>;
template class Read_Write<RIGID_GEOMETRY_COLLECTION<VECTOR<double,2> >,double>;
template class Read_Write<RIGID_GEOMETRY_COLLECTION<VECTOR<double,3> >,double>;
template class Read_Write<RIGID_GEOMETRY_COLLECTION<VECTOR<double,1> >,float>;
template class Read_Write<RIGID_GEOMETRY_COLLECTION<VECTOR<double,2> >,float>;
template class Read_Write<RIGID_GEOMETRY_COLLECTION<VECTOR<double,3> >,float>;
template class Read_Write<RIGID_GEOMETRY_COLLECTION<VECTOR<float,1> >,double>;
template class Read_Write<RIGID_GEOMETRY_COLLECTION<VECTOR<float,2> >,double>;
template class Read_Write<RIGID_GEOMETRY_COLLECTION<VECTOR<float,3> >,double>;
#endif
