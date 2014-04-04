//#####################################################################
// Copyright 2003-2009, Zhaosheng Bao, Kevin Der, Ron Fedkiw, Eran Guendelman, Geoffrey Irving, Nipun Kwatra, Cynthia Lau, Michael Lentine, Sergey Levine, Frank Losasso, Nick Rasmussen, Craig Schroeder, Andrew Selle, Tamar Shinar, Eftychios Sifakis, Jonathan Su, Jerry Talton, Joseph Teran, Rachel Weinstein.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#include <PhysBAM_Tools/Grids_RLE/RLE_GRID_ITERATOR_CELL_3D.h>
#include <PhysBAM_Tools/Grids_Uniform/GRID.h>
#include <PhysBAM_Tools/Parsing/PARAMETER_LIST.h>
#include <PhysBAM_Tools/Parsing/PARSE_ARGS.h>
#include <PhysBAM_Tools/Read_Write/Arrays/READ_WRITE_ARRAY.h>
#include <PhysBAM_Tools/Read_Write/Grids_Uniform/READ_WRITE_GRID.h>
#include <PhysBAM_Tools/Read_Write/Grids_Uniform_Arrays/READ_WRITE_ARRAYS.h>
#include <PhysBAM_Tools/Read_Write/Grids_Uniform_Arrays/READ_WRITE_FACE_ARRAYS.h>
#include <PhysBAM_Tools/Read_Write/Utilities/FILE_UTILITIES.h>
#include <PhysBAM_Tools/Utilities/PROCESS_UTILITIES.h>
#include <PhysBAM_Geometry/Basic_Geometry/BASIC_GEOMETRY_FORWARD.h>
#include <PhysBAM_Geometry/Geometry_Particles/GEOMETRY_PARTICLES_FORWARD.h>
#include <PhysBAM_Geometry/Geometry_Particles/REGISTER_GEOMETRY_READ_WRITE.h>
#include <PhysBAM_Geometry/Grids_Uniform_Level_Sets/LEVELSET_3D.h>
#include <PhysBAM_Geometry/Read_Write/Grids_Uniform_Level_Sets/READ_WRITE_LEVELSET_3D.h>
#include <PhysBAM_Rendering/PhysBAM_OpenGL/OpenGL/ANIMATED_VISUALIZATION.h>
#include <PhysBAM_Rendering/PhysBAM_OpenGL/OpenGL/OPENGL_BOOL_COLOR_MAP.h>
#include <PhysBAM_Rendering/PhysBAM_OpenGL/OpenGL/OPENGL_BOX_3D.h>
#include <PhysBAM_Rendering/PhysBAM_OpenGL/OpenGL/OPENGL_COLOR_RAMP.h>
#include <PhysBAM_Rendering/PhysBAM_OpenGL/OpenGL/OPENGL_CONSTANT_COLOR_MAP.h>
#include <PhysBAM_Rendering/PhysBAM_OpenGL/OpenGL/OPENGL_GRID_3D.h>
#include <PhysBAM_Rendering/PhysBAM_OpenGL/OpenGL/OPENGL_INDEXED_COLOR_MAP.h>
#include <PhysBAM_Rendering/PhysBAM_OpenGL/OpenGL/OPENGL_LEVELSET_COLOR_MAP.h>
#include <PhysBAM_Rendering/PhysBAM_OpenGL/OpenGL/OPENGL_OCTREE_SLICE.h>
#include <PhysBAM_Rendering/PhysBAM_OpenGL/OpenGL/OPENGL_RLE_GRID_2D.h>
#include <PhysBAM_Rendering/PhysBAM_OpenGL/OpenGL/OPENGL_RLE_GRID_3D.h>
#include <PhysBAM_Rendering/PhysBAM_OpenGL/OpenGL/OPENGL_RLE_SLICE.h>
#include <PhysBAM_Rendering/PhysBAM_OpenGL/OpenGL/OPENGL_SELECTION.h>
#include <PhysBAM_Rendering/PhysBAM_OpenGL/OpenGL/OPENGL_SLICE_MANAGER.h>
#include <PhysBAM_Rendering/PhysBAM_OpenGL/OpenGL/OPENGL_TRIANGULATED_SURFACE.h>
#include <PhysBAM_Rendering/PhysBAM_OpenGL/OpenGL/OPENGL_UNIFORM_SLICE.h>
#include <PhysBAM_Rendering/PhysBAM_OpenGL/OpenGL/OPENGL_WORLD.h>
#include <PhysBAM_Rendering/PhysBAM_OpenGL/OpenGL_Components/OPENGL_COMPONENT_BASIC.h>
#include <PhysBAM_Rendering/PhysBAM_OpenGL/OpenGL_Components/OPENGL_COMPONENT_DEFORMABLE_GEOMETRY_COLLECTION_3D.h>
#include <PhysBAM_Rendering/PhysBAM_OpenGL/OpenGL_Components/OPENGL_COMPONENT_DIAGNOSTICS.h>
#include <PhysBAM_Rendering/PhysBAM_OpenGL/OpenGL_Components/OPENGL_COMPONENT_FACE_SCALAR_FIELD_3D.h>
#include <PhysBAM_Rendering/PhysBAM_OpenGL/OpenGL_Components/OPENGL_COMPONENT_GRID_BASED_VECTOR_FIELD_3D.h>
#include <PhysBAM_Rendering/PhysBAM_OpenGL/OpenGL_Components/OPENGL_COMPONENT_LEVELSET_3D.h>
#include <PhysBAM_Rendering/PhysBAM_OpenGL/OpenGL_Components/OPENGL_COMPONENT_MAC_VELOCITY_FIELD_3D.h>
#include <PhysBAM_Rendering/PhysBAM_OpenGL/OpenGL_Components/OPENGL_COMPONENT_OCTREE_CELL_SCALAR_FIELD.h>
#include <PhysBAM_Rendering/PhysBAM_OpenGL/OpenGL_Components/OPENGL_COMPONENT_OCTREE_FACE_SCALAR_FIELD.h>
#include <PhysBAM_Rendering/PhysBAM_OpenGL/OpenGL_Components/OPENGL_COMPONENT_OCTREE_GRID.h>
#include <PhysBAM_Rendering/PhysBAM_OpenGL/OpenGL_Components/OPENGL_COMPONENT_OCTREE_NODE_SCALAR_FIELD.h>
#include <PhysBAM_Rendering/PhysBAM_OpenGL/OpenGL_Components/OPENGL_COMPONENT_OCTREE_NODE_VECTOR_FIELD.h>
#include <PhysBAM_Rendering/PhysBAM_OpenGL/OpenGL_Components/OPENGL_COMPONENT_PARTICLES_3D.h>
#include <PhysBAM_Rendering/PhysBAM_OpenGL/OpenGL_Components/OPENGL_COMPONENT_PSEUDO_DIRICHLET_3D.h>
#include <PhysBAM_Rendering/PhysBAM_OpenGL/OpenGL_Components/OPENGL_COMPONENT_RIGID_GEOMETRY_COLLECTION_3D.h>
#include <PhysBAM_Rendering/PhysBAM_OpenGL/OpenGL_Components/OPENGL_COMPONENT_RLE_CELL_SCALAR_FIELD_3D.h>
#include <PhysBAM_Rendering/PhysBAM_OpenGL/OpenGL_Components/OPENGL_COMPONENT_RLE_FACE_SCALAR_FIELD_3D.h>
#include <PhysBAM_Rendering/PhysBAM_OpenGL/OpenGL_Components/OPENGL_COMPONENT_RLE_GRID_3D.h>
#include <PhysBAM_Rendering/PhysBAM_OpenGL/OpenGL_Components/OPENGL_COMPONENT_SCALAR_FIELD_3D.h>
#include <PhysBAM_Rendering/PhysBAM_OpenGL/OpenGL_Components/OPENGL_COMPONENT_SYMMETRIC_MATRIX_FIELD_3D.h>
#include <PhysBAM_Rendering/PhysBAM_OpenGL/OpenGL_Components/OPENGL_COMPONENT_THIN_SHELLS_DEBUGGING_3D.h>
#include <PhysBAM_Rendering/PhysBAM_OpenGL/OpenGL_Components/OPENGL_COMPONENT_TRIANGULATED_SURFACE.h>
#include <fstream>
#include <sstream>
#include <string.h>

using namespace PhysBAM;
using namespace std;

template<class T,class RW=T>
class VISUALIZATION:public ANIMATED_VISUALIZATION
{
    typedef VECTOR<T,3> TV;
public:
    VISUALIZATION();
    ~VISUALIZATION();

protected:
    virtual void Add_Arguments(PARSE_ARGS &parse_args);
    virtual void Parse_Arguments(PARSE_ARGS &parse_args);
    virtual void Initialize_Components_And_Key_Bindings();
    virtual void Update_OpenGL_Strings();

private:
    void Read_Grid();
    virtual void Pre_Frame_Extra();
    virtual void Set_Frame_Extra();

    // callbacks
    void Command_Prompt();
    void Command_Prompt_Response();
    void Slice_Has_Changed();
    DEFINE_CALLBACK_CREATOR(VISUALIZATION,Command_Prompt);
    DEFINE_CALLBACK_CREATOR(VISUALIZATION,Command_Prompt_Response);
    DEFINE_CALLBACK_CREATOR(VISUALIZATION,Slice_Has_Changed);

    // TODO: need better grid control 
#ifndef COMPILE_WITHOUT_RLE_SUPPORT
    OPENGL_COMPONENT_RLE_GRID_3D<T>* rle_grid_component;
#endif
    OPENGL_COMPONENT_BASIC<OPENGL_GRID_3D<T> >* grid_component;
    OPENGL_COMPONENT_BASIC<OPENGL_GRID_3D<T> >* coarse_grid_component;
#ifndef COMPILE_WITHOUT_DYADIC_SUPPORT
    OPENGL_COMPONENT_OCTREE_GRID<T>* octree_grid_component;
#endif

    // Options
    std::string basedir;
    bool has_diagnostics;

    GRID<TV> grid,mac_grid,regular_grid,coarse_grid,coarse_mac_grid,coarse_regular_grid;
    bool has_valid_grid,has_valid_rle_grid,has_valid_octree_grid,has_valid_coarse_grid;
    bool node_based,coarse_node_based;
    OPENGL_SLICE_MANAGER slice_manager;

    ARRAY<int> rigid_bodies_no_draw_list;
    ARRAY<int> deformable_no_draw_list;

    bool allow_caching;
    bool always_add_mac_velocities;
};

// ------------------------------------------------------------------

template<class T,class RW> VISUALIZATION<T,RW>::
VISUALIZATION()
    :ANIMATED_VISUALIZATION()
#ifndef COMPILE_WITHOUT_RLE_SUPPORT
    ,rle_grid_component(0)
#endif
    ,grid_component(0),coarse_grid_component(0)
#ifndef COMPILE_WITHOUT_DYADIC_SUPPORT
    ,octree_grid_component(0)
#endif
{
}

template<class T,class RW> VISUALIZATION<T,RW>::
~VISUALIZATION()
{
    if(grid_component){delete &grid_component->object;delete grid_component;}
    if(coarse_grid_component){delete &coarse_grid_component->object;delete coarse_grid_component;}
#ifndef COMPILE_WITHOUT_RLE_SUPPORT
    delete rle_grid_component;
#endif
}

template<class T,class RW> void VISUALIZATION<T,RW>::
Add_Arguments(PARSE_ARGS &parse_args)
{
    basedir=".";

    ANIMATED_VISUALIZATION::Add_Arguments(parse_args);

    parse_args.Add_Option_Argument("-float");
    parse_args.Add_Option_Argument("-double");
    parse_args.Add_Option_Argument("-no_caching");
    parse_args.Add_Option_Argument("-macvel","force adding mac velocities component");
    parse_args.Add_String_Argument("-rigid_bodies_no_draw","");
    parse_args.Add_String_Argument("-deformable_no_draw","");
    parse_args.Set_Extra_Arguments(-1,"[<basedir>]");
}

template<class T,class RW> void VISUALIZATION<T,RW>::
Parse_Arguments(PARSE_ARGS &parse_args)
{
    ANIMATED_VISUALIZATION::Parse_Arguments(parse_args);

    if(parse_args.Num_Extra_Args() > 1)
        parse_args.Print_Usage(true);
    else if(parse_args.Num_Extra_Args()==1)
        basedir=parse_args.Extra_Arg(1);

    if(parse_args.Is_Value_Set("-rigid_bodies_no_draw")){ARRAY<int> int_list;
        STRING_UTILITIES::Parse_Integer_List(parse_args.Get_String_Value("-rigid_bodies_no_draw"),int_list);
        rigid_bodies_no_draw_list.Resize(int_list.Size());
        for(int i=1;i<=rigid_bodies_no_draw_list.m;i++) rigid_bodies_no_draw_list(i)=int(int_list(i));}

    STRING_UTILITIES::Parse_Integer_List(parse_args.Get_String_Value("-deformable_no_draw"),deformable_no_draw_list);
    for(int i=deformable_no_draw_list.m;i;i--)
        if(deformable_no_draw_list(i)<1)
            deformable_no_draw_list.Remove_Index_Lazy(i);

    last_frame_filename=basedir+"/common/last_frame";

    if(FILE_UTILITIES::File_Exists(basedir+"/common/first_frame") && !parse_args.Is_Value_Set("-start_frame")) FILE_UTILITIES::Read_From_Text_File(basedir+"/common/first_frame",start_frame);

    allow_caching=!parse_args.Get_Option_Value("-no_caching");
    always_add_mac_velocities=parse_args.Is_Value_Set("-macvel");

    // don't override camera script filename if it was already set in base class based on command line argument
    if(camera_script_filename.empty()) camera_script_filename=basedir+"/camera_script";

#ifdef __linux__
    if(!parse_args.Is_Value_Set("-window_title")) opengl_window_title="opengl_3d: " + FILE_UTILITIES::Real_Path(basedir);
#endif
}

template<class T,class RW> void VISUALIZATION<T,RW>::
Read_Grid()
{
    has_valid_grid=false;
    has_valid_coarse_grid=false;
#ifndef COMPILE_WITHOUT_RLE_SUPPORT
    has_valid_rle_grid=false;
#endif
    has_valid_octree_grid=false;
#ifndef COMPILE_WITHOUT_RLE_SUPPORT
    std::string rle_filename;
#endif
    std::string octree_filename;
    std::string filename,coarse_filename;

    std::string octree_grid_prefix="";
#ifndef COMPILE_WITHOUT_RLE_SUPPORT
    rle_filename=STRING_UTILITIES::string_sprintf("%s/%d/rle_grid",basedir.c_str(),start_frame);
#endif
    octree_filename=STRING_UTILITIES::string_sprintf("%s/%d/octree_grid",basedir.c_str(),start_frame);
    filename=STRING_UTILITIES::string_sprintf("%s/%d/levelset",basedir.c_str(),start_frame);
    coarse_filename=STRING_UTILITIES::string_sprintf("%s/%d/coarse_levelset",basedir.c_str(),start_frame);
    // For backwards compatibility
    if(!FILE_UTILITIES::File_Exists(filename)) filename=STRING_UTILITIES::string_sprintf("%s/%d/levelset.phi",basedir.c_str(),start_frame);

#ifndef COMPILE_WITHOUT_RLE_SUPPORT
    if(FILE_UTILITIES::File_Exists(rle_filename)){
        std::cout<<"Reading rle grid '"<<rle_filename<<"'"<<std::endl<<std::flush;
        rle_filename=STRING_UTILITIES::string_sprintf("%s/%d/rle_grid",basedir.c_str());
        rle_grid_component=new OPENGL_COMPONENT_RLE_GRID_3D<T>(rle_filename);
        // Add_Component(rle_grid_component); // add later instead of here, but Pre_Frame_Extra will ensure rle grid is read in
        has_valid_rle_grid=true;}
#endif
#ifndef COMPILE_WITHOUT_DYADIC_SUPPORT
    if(FILE_UTILITIES::File_Exists(octree_filename)){
        std::cout<<"Reading octree '"<<octree_filename<<"'"<<std::endl<<std::flush;
        octree_filename=STRING_UTILITIES::string_sprintf("%s/%d/octree_grid",basedir.c_str());
        octree_grid_component=new OPENGL_COMPONENT_OCTREE_GRID<T>(octree_filename);
        //Add_Component(octree_grid_component); // add later instead of here, but Pre_Frame_Extra will ensure octree grid is read in
        has_valid_octree_grid=true;}
#endif

    if(FILE_UTILITIES::File_Exists(coarse_filename)){
        std::cout<<"Reading coarse_grid from '"<<coarse_filename<<"'..."<<std::endl;
        ARRAY<T,VECTOR<int,3> > phi;
        LEVELSET_3D<GRID<TV> > levelset(coarse_grid,phi);
        FILE_UTILITIES::Read_From_File<RW>(coarse_filename,levelset);
        has_valid_coarse_grid=true;}
    else if(FILE_UTILITIES::File_Exists(basedir+"/common/coarse_grid")){
        coarse_filename=basedir+"/common/coarse_grid";
        LOG::cout<<"Reading coarse grid from '"<<coarse_filename<<"'..."<<std::flush;
        FILE_UTILITIES::Read_From_File<RW>(coarse_filename,coarse_grid);
        has_valid_coarse_grid=true;}

    if(FILE_UTILITIES::File_Exists(filename)){
        std::cout<<"Reading grid from '"<<filename<<"'..."<<std::endl;
        ARRAY<T,VECTOR<int,3> > phi;
        LEVELSET_3D<GRID<TV> > levelset(grid,phi);
        FILE_UTILITIES::Read_From_File<RW>(filename,levelset);
        has_valid_grid=true;}
    else if(FILE_UTILITIES::File_Exists(basedir+"/common/grid")){
        filename=basedir+"/common/grid";
        std::cout<<"Reading grid from '"<<filename<<"'..."<<std::endl;
        FILE_UTILITIES::Read_From_File<RW>(filename,grid);
        has_valid_grid=true;}

    if(has_valid_grid){
        node_based=!grid.Is_MAC_Grid();
        mac_grid=grid.Get_MAC_Grid();regular_grid=grid.Get_Regular_Grid();}

    if(has_valid_coarse_grid){
        coarse_node_based=!coarse_grid.Is_MAC_Grid();
        coarse_mac_grid=coarse_grid.Get_MAC_Grid();coarse_regular_grid=coarse_grid.Get_Regular_Grid();}
}

//#####################################################################
// Initialize_Components_And_Key_Bindings
//#####################################################################
template<class T,class RW> void VISUALIZATION<T,RW>::
Initialize_Components_And_Key_Bindings()
{
    ANIMATED_VISUALIZATION::Initialize_Components_And_Key_Bindings();
    opengl_world.Set_Key_Binding_Category_Priority(1);
    opengl_world.Unbind_Keys("abBCdDEeFjJjkKlLMotTvV 1!2@3#4$5%67&89 ^=-`{}\b\\[]~\t");

    std::string filename,filename2,coarse_filename;
    filename=basedir+"/common/sim.param";
    PARAMETER_LIST parameter_list;
    parameter_list.Read(filename);

    frame_rate=parameter_list.template Get_Parameter<int>("frame_rate",frame_rate);
    std::cout<<"Using frame rate "<<frame_rate<<std::endl;


    filename=basedir+"/%d/diagnostics";
    if(FILE_UTILITIES::Frame_File_Exists(filename,start_frame))
        Add_Component(new OPENGL_COMPONENT_DIAGNOSTICS(filename),"Diagnostics",'\0',BASIC_VISUALIZATION::OWNED);

    Read_Grid();
    if(has_valid_grid){
        OPENGL_UNIFORM_SLICE* slice=new OPENGL_UNIFORM_SLICE(opengl_world);slice->Initialize(grid);
        slice_manager.slice=slice;
        std::cout<<"Using uniform grid slice"<<std::endl;}

#ifndef COMPILE_WITHOUT_RLE_SUPPORT
    if(has_valid_rle_grid){
        OPENGL_RLE_SLICE* slice=new OPENGL_RLE_SLICE(opengl_world);slice->Initialize(rle_grid_component->opengl_grid.grid);
        slice_manager.slice=slice;slice_manager.Add_Object(&rle_grid_component->opengl_grid);
        std::cout<<"Using rle grid slice"<<std::endl;}
#endif

#ifndef COMPILE_WITHOUT_DYADIC_SUPPORT
    if(has_valid_octree_grid){
        OPENGL_OCTREE_SLICE* slice=new OPENGL_OCTREE_SLICE(opengl_world);slice->Initialize(octree_grid_component->opengl_grid.grid);
        slice_manager.slice=slice;slice_manager.Add_Object(&octree_grid_component->opengl_grid);
        std::cout<<"Using octree grid slice"<<std::endl;}
#endif

    if(slice_manager.slice) slice_manager.Set_Slice_Has_Changed_Callback(Slice_Has_Changed_CB());

    if(has_valid_grid){
        OPENGL_BOX_3D<T>* opengl_box=new OPENGL_BOX_3D<T>(*(new BOX<TV>(grid.Domain())),OPENGL_COLOR::Gray(0.5));
        OPENGL_COMPONENT_BASIC<OPENGL_BOX_3D<T> >* domain_box_component=new OPENGL_COMPONENT_BASIC<OPENGL_BOX_3D<T> >(*opengl_box);
        Add_Component(domain_box_component,"Domain box",'6',BASIC_VISUALIZATION::OWNED|BASIC_VISUALIZATION::START_HIDDEN);}

    {OPENGL_COMPONENT_DEFORMABLE_GEOMETRY_COLLECTION_3D<T>* deformable_objects_component=0;
    std::string deformable_object_filename=basedir+STRING_UTILITIES::string_sprintf("/%d/deformable_object_structures",start_frame);
    if(FILE_UTILITIES::File_Exists(basedir+"/common/deformable_object_structures") || FILE_UTILITIES::File_Exists(deformable_object_filename)){
        OPENGL_MATERIAL front_material=OPENGL_MATERIAL::Matte(OPENGL_COLOR::Yellow());
        OPENGL_MATERIAL back_material=OPENGL_MATERIAL::Matte(OPENGL_COLOR::Yellow(0.5));
        deformable_objects_component=new OPENGL_COMPONENT_DEFORMABLE_GEOMETRY_COLLECTION_3D<T>(basedir+"/",start_frame);
        deformable_objects_component->Set_All_Materials(front_material,front_material,back_material);
        if(deformable_no_draw_list.m){
            deformable_objects_component->active_list.Resize(PhysBAM::max(deformable_objects_component->active_list.m,deformable_no_draw_list.Last()),true,true,true);
            INDIRECT_ARRAY<ARRAY<bool>,ARRAY<int>&> subset=deformable_objects_component->active_list.Subset(deformable_no_draw_list);ARRAYS_COMPUTATIONS::Fill(subset,false);}
        opengl_world.Set_Key_Binding_Category("Deformable Objects");
        Add_Component(deformable_objects_component,"Deformable Objects",'8',BASIC_VISUALIZATION::OWNED|BASIC_VISUALIZATION::SELECTABLE);
        opengl_world.Append_Bind_Key('9',deformable_objects_component->Toggle_Active_Value_CB());
        opengl_world.Append_Bind_Key('(',deformable_objects_component->Show_Only_First_CB());
        opengl_world.Append_Bind_Key('w',deformable_objects_component->Toggle_Selection_Mode_CB());
        opengl_world.Append_Bind_Key('h',deformable_objects_component->Toggle_Hide_Unselected_CB());
        opengl_world.Append_Bind_Key('b',deformable_objects_component->Toggle_Draw_Interior_CB());
        opengl_world.Append_Bind_Key('t',deformable_objects_component->Toggle_Draw_Subsets_CB());
        opengl_world.Append_Bind_Key('e',deformable_objects_component->Cycle_Display_Mode_CB());
        opengl_world.Append_Bind_Key('$',deformable_objects_component->Toggle_Velocity_Vectors_CB());
        opengl_world.Append_Bind_Key('=',deformable_objects_component->Increase_Vector_Size_CB());
        opengl_world.Append_Bind_Key('-',deformable_objects_component->Decrease_Vector_Size_CB());
        opengl_world.Append_Bind_Key(OPENGL_KEY(OPENGL_KEY::F7),deformable_objects_component->Cycle_Relative_Velocity_Mode_CB());
        if(deformable_objects_component->has_tetrahedralized_volumes){
            opengl_world.Append_Bind_Key('C',deformable_objects_component->Cycle_Cutaway_Mode_CB());
            opengl_world.Append_Bind_Key('{',deformable_objects_component->Decrease_Cutaway_Fraction_CB());
            opengl_world.Append_Bind_Key('}',deformable_objects_component->Increase_Cutaway_Fraction_CB());}
        opengl_world.Append_Bind_Key('Z',deformable_objects_component->Highlight_Particle_CB());
        if(slice_manager.slice) slice_manager.Add_Object(deformable_objects_component);}

    if(FILE_UTILITIES::Frame_File_Exists(basedir+"/%d/rigid_geometry_particles",start_frame)){
        OPENGL_COMPONENT_RIGID_GEOMETRY_COLLECTION_3D<T>* rigid_bodies_component=0;
        //if(deformable_objects_component) rigid_bodies_component=new OPENGL_COMPONENT_RIGID_BODY_COLLECTION_3D<T>(deformable_objects_component->solid_body_collection.rigid_body_collection,basedir,true);
        //else rigid_bodies_component=new OPENGL_COMPONENT_RIGID_BODY_COLLECTION_3D<T>(basedir,true);
        rigid_bodies_component=new OPENGL_COMPONENT_RIGID_GEOMETRY_COLLECTION_3D<T>(basedir,true);
        rigid_bodies_component->Set_Vector_Size(0.01);
        rigid_bodies_component->Set_Frame(start_frame); // needed before reinitialize so that no draw list will work
        rigid_bodies_component->Reinitialize();
        for(int i=1;i<=rigid_bodies_no_draw_list.m;i++){
            std::cout<<"Rigid bodies: not drawing object "<<rigid_bodies_no_draw_list(i)<<std::endl;
            rigid_bodies_component->Set_Draw_Object(rigid_bodies_no_draw_list(i),false);}
        // Read hints if available
        std::string hints_filename=basedir+"/common/opengl_hints";
        if(FILE_UTILITIES::File_Exists(hints_filename)){
            std::cout<<"Using opengl rigid body hints"<<std::endl;
            rigid_bodies_component->Read_Hints(hints_filename);}
        rigid_bodies_component->Reinitialize();
        opengl_world.Set_Key_Binding_Category("Rigid Bodies");
        Add_Component(rigid_bodies_component,"Rigid Bodies",'\0',BASIC_VISUALIZATION::OWNED|BASIC_VISUALIZATION::SELECTABLE);
        opengl_world.Append_Bind_Key('5',rigid_bodies_component->Toggle_Draw_Mode_CB());
        opengl_world.Append_Bind_Key('%',rigid_bodies_component->Toggle_Velocity_Vectors_CB());
        opengl_world.Append_Bind_Key('a',rigid_bodies_component->Toggle_Individual_Axes_CB());
        opengl_world.Append_Bind_Key('%',rigid_bodies_component->Toggle_Show_Object_Names_CB());
        opengl_world.Append_Bind_Key('=',rigid_bodies_component->Increase_Vector_Size_CB());
        opengl_world.Append_Bind_Key('-',rigid_bodies_component->Decrease_Vector_Size_CB());
        opengl_world.Append_Bind_Key('M',rigid_bodies_component->Toggle_Draw_Particles_CB());
        opengl_world.Append_Bind_Key('o',rigid_bodies_component->Toggle_One_Sided_CB());
        if(slice_manager.slice) slice_manager.Add_Object(rigid_bodies_component);}}

    std::string soft_constraints_deformable_object_filename=basedir+"/soft_constraints_deformable_object_particles";
    if(FILE_UTILITIES::File_Exists(soft_constraints_deformable_object_filename) // TODO(jontg): Not sure what to do here...
        || FILE_UTILITIES::File_Exists(soft_constraints_deformable_object_filename+STRING_UTILITIES::string_sprintf(".%d",start_frame))){
        OPENGL_COMPONENT_DEFORMABLE_GEOMETRY_COLLECTION_3D<T>* soft_constraints_deformable_objects_component=new OPENGL_COMPONENT_DEFORMABLE_GEOMETRY_COLLECTION_3D<T>(basedir+"/soft_constraints_",start_frame);
        soft_constraints_deformable_objects_component->selectable=true;
        Add_Component(soft_constraints_deformable_objects_component,"Soft Constraints Deformable Objects",'\0',BASIC_VISUALIZATION::OWNED|BASIC_VISUALIZATION::SELECTABLE);
        if(slice_manager.slice) slice_manager.Add_Object(soft_constraints_deformable_objects_component);}

    filename=basedir+"/%d/levelset";
    coarse_filename=basedir+"/%d/coarse_levelset";
    filename2=allow_caching?(basedir+"/%d/levelset_tesselated"):"";
    // for backwards compatiblity
    if(!FILE_UTILITIES::Frame_File_Exists(filename,start_frame) && !FILE_UTILITIES::Frame_File_Exists(basedir+"/%d/levelset_1",start_frame)){
        filename=basedir+"/levelset_%d.phi";
        filename2=allow_caching?(basedir+"/levelset_%d.tri"):"";}
    if(FILE_UTILITIES::Frame_File_Exists(filename,start_frame)||FILE_UTILITIES::Frame_File_Exists(basedir+"/%d/levelset_1",start_frame)){
        OPENGL_COMPONENT_LEVELSET_3D<T>* levelset_component=new OPENGL_COMPONENT_LEVELSET_3D<T>(filename,filename2,basedir+"/levelset_%d.%d",basedir+"/levelset_tesselated_%d.%d",true);
        opengl_world.Set_Key_Binding_Category("Levelset");
        Add_Component(levelset_component,"Levelset",'l',BASIC_VISUALIZATION::OWNED);
//        levelset_component->Set_Surface_Material(OPENGL_MATERIAL::Plastic(OPENGL_COLOR((T).6,(T).65,1)),OPENGL_MATERIAL::Plastic(OPENGL_COLOR((T).6,(T).65,1)));
        if(slice_manager.slice) slice_manager.Add_Object(levelset_component);
        opengl_world.Append_Bind_Key('L',levelset_component->Toggle_Slice_Color_Mode_CB());
        opengl_world.Append_Bind_Key("^l",levelset_component->Toggle_Display_Overlay_CB());
        opengl_world.Append_Bind_Key('`',levelset_component->Toggle_Smooth_Slice_CB());
        if(levelset_component->Use_Sets()){
            opengl_world.Append_Bind_Key('M',levelset_component->Toggle_Draw_Multiple_Levelsets_CB());
            opengl_world.Append_Bind_Key('>',levelset_component->Next_Set_CB());
            opengl_world.Append_Bind_Key('<',levelset_component->Previous_Set_CB());}}
    if(FILE_UTILITIES::Frame_File_Exists(coarse_filename,start_frame)){
        OPENGL_COMPONENT_LEVELSET_3D<T>* levelset_component=new OPENGL_COMPONENT_LEVELSET_3D<T>(coarse_filename,filename2,basedir+"/coarse_levelset_%d.%d",basedir+"/coarse_levelset_tesselated_%d.%d",true);
        levelset_component->ghost_cells=0;
        opengl_world.Set_Key_Binding_Category("Levelset");
        Add_Component(levelset_component,"Levelset coarse",'.',BASIC_VISUALIZATION::OWNED);
        levelset_component->Set_Surface_Material(OPENGL_MATERIAL::Plastic(OPENGL_COLOR((T).6,(T)1,(T).65)),OPENGL_MATERIAL::Plastic(OPENGL_COLOR((T).6,(T)1,(T).65)));
        if(slice_manager.slice) slice_manager.Add_Object(levelset_component);
        opengl_world.Append_Bind_Key('>',levelset_component->Toggle_Slice_Color_Mode_CB());
        opengl_world.Append_Bind_Key("^.",levelset_component->Toggle_Display_Overlay_CB());}

    filename=basedir+"/%d/object_levelset";
    filename2=basedir+"/%d/object_levelset_tesselated";
    if(FILE_UTILITIES::Frame_File_Exists(filename,start_frame)){
        OPENGL_COMPONENT_LEVELSET_3D<T>* object_levelset_component=new OPENGL_COMPONENT_LEVELSET_3D<T>(filename,filename2,"","",true);
        Add_Component(object_levelset_component,"Object Levelset",'\0',BASIC_VISUALIZATION::OWNED|BASIC_VISUALIZATION::START_HIDDEN);
        if(slice_manager.slice) slice_manager.Add_Object(object_levelset_component);
        opengl_world.Set_Key_Binding_Category("Object Levelset");
        opengl_world.Append_Bind_Key("^l",object_levelset_component->Toggle_Draw_CB());
        opengl_world.Append_Bind_Key('L',object_levelset_component->Toggle_Slice_Color_Mode_CB());
        opengl_world.Append_Bind_Key('`',object_levelset_component->Toggle_Smooth_Slice_CB());}

#ifndef COMPILE_WITHOUT_RLE_SUPPORT
    filename=basedir+"/%d/rle_levelset";
    if(has_valid_rle_grid && FILE_UTILITIES::Frame_File_Exists(filename,start_frame)){
        OPENGL_COMPONENT_RLE_CELL_SCALAR_FIELD_3D<T>* rle_levelset_component=new OPENGL_COMPONENT_RLE_CELL_SCALAR_FIELD_3D<T>(rle_grid_component->opengl_grid.grid,filename,
            new OPENGL_LEVELSET_COLOR_MAP<T>(OPENGL_COLOR::Blue(),OPENGL_COLOR::Red()));
        rle_levelset_component->opengl_rle_cell_scalar_field.draw_filled=true;
        opengl_world.Set_Key_Binding_Category("RLE Levelset");
        Add_Component(rle_levelset_component,"RLE Levelset",'l',BASIC_VISUALIZATION::OWNED);
        opengl_world.Append_Bind_Key('L',rle_levelset_component->Draw_Surface_CB());
        slice_manager.Add_Object(rle_levelset_component);}

    filename=basedir+"/%d/rle_object_levelset";
    if(has_valid_rle_grid && FILE_UTILITIES::Frame_File_Exists(filename,start_frame)){
        OPENGL_COMPONENT_RLE_CELL_SCALAR_FIELD_3D<T>* rle_object_levelset_component=new OPENGL_COMPONENT_RLE_CELL_SCALAR_FIELD_3D<T>(rle_grid_component->opengl_grid.grid,filename,
            new OPENGL_LEVELSET_COLOR_MAP<T>(OPENGL_COLOR::Blue(),OPENGL_COLOR::Red()));
        rle_object_levelset_component->opengl_rle_cell_scalar_field.draw_filled=true;
        rle_object_levelset_component->Set_Draw(false);
        Add_Component(rle_object_levelset_component,"RLE Object Levelset",'\0',BASIC_VISUALIZATION::OWNED);
        opengl_world.Set_Key_Binding_Category("RLE Object Levelset");
        opengl_world.Append_Bind_Key("^l",rle_object_levelset_component->Toggle_Draw_CB());
        opengl_world.Append_Bind_Key('L',rle_object_levelset_component->Draw_Surface_CB());
        slice_manager.Add_Object(rle_object_levelset_component);}
#endif

#ifndef COMPILE_WITHOUT_DYADIC_SUPPORT
    filename=basedir+"/%d/octree_levelset";
    if(has_valid_octree_grid && FILE_UTILITIES::Frame_File_Exists(filename,start_frame)){
        OPENGL_COMPONENT_OCTREE_CELL_SCALAR_FIELD<T>* octree_levelset_component=new OPENGL_COMPONENT_OCTREE_CELL_SCALAR_FIELD<T>(octree_grid_component->opengl_grid.grid,filename,
            new OPENGL_LEVELSET_COLOR_MAP<T>(OPENGL_COLOR::Blue(),OPENGL_COLOR::Red()));
        opengl_world.Set_Key_Binding_Category("Octree Level Set");
        Add_Component(octree_levelset_component,"Octree Levelset",'l',BASIC_VISUALIZATION::OWNED);
        opengl_world.Append_Bind_Key('+',octree_levelset_component->Increase_Point_Size_CB());
        opengl_world.Append_Bind_Key('_',octree_levelset_component->Decrease_Point_Size_CB());
        opengl_world.Append_Bind_Key('L',octree_levelset_component->Draw_Surface_CB());
        slice_manager.Add_Object(octree_levelset_component);}

    filename=basedir+"/%d/octree_levelset_nodes";
    if(has_valid_octree_grid && FILE_UTILITIES::Frame_File_Exists(filename,start_frame)){
        OPENGL_COMPONENT_OCTREE_NODE_SCALAR_FIELD<T>* octree_levelset_nodes_component=new OPENGL_COMPONENT_OCTREE_NODE_SCALAR_FIELD<T>(octree_grid_component->opengl_grid.grid,filename,
            new OPENGL_LEVELSET_COLOR_MAP<T>(OPENGL_COLOR::Blue(),OPENGL_COLOR::Red()));
        opengl_world.Set_Key_Binding_Category("Octree Nodal Level Set");
        Add_Component(octree_levelset_nodes_component,"Nodal Octree Levelset",';',BASIC_VISUALIZATION::OWNED);
        opengl_world.Append_Bind_Key('+',octree_levelset_nodes_component->Increase_Point_Size_CB());
        opengl_world.Append_Bind_Key('_',octree_levelset_nodes_component->Decrease_Point_Size_CB());
        opengl_world.Append_Bind_Key(':',octree_levelset_nodes_component->Draw_Surface_CB());
        slice_manager.Add_Object(octree_levelset_nodes_component);}
#endif

    opengl_world.Set_Key_Binding_Category("Density");
    filename=basedir+"/%d/density";
    if(has_valid_grid && FILE_UTILITIES::Frame_File_Exists(filename,start_frame)){
        OPENGL_COMPONENT_SCALAR_FIELD_3D<T>* density_component=new OPENGL_COMPONENT_SCALAR_FIELD_3D<T>(grid,filename,OPENGL_COLOR_RAMP<T>::Two_Color_Ramp(0,1,OPENGL_COLOR::Gray(0,0),OPENGL_COLOR::Gray(1,1)));
        density_component->opengl_scalar_field.Update();
        Add_Component(density_component,"Density",'d',BASIC_VISUALIZATION::OWNED);
        opengl_world.Append_Bind_Key('D',density_component->Toggle_Color_Map_CB());
        opengl_world.Append_Bind_Key('`',density_component->Toggle_Smooth_Slice_CB());
        slice_manager.Add_Object(density_component);}

    filename=basedir+"/%d/density_gradient";
    if(has_valid_grid && FILE_UTILITIES::Frame_File_Exists(filename,start_frame)){
        OPENGL_COMPONENT_SCALAR_FIELD_3D<T>* density_gradient_component=new OPENGL_COMPONENT_SCALAR_FIELD_3D<T>(grid,filename,
            OPENGL_COLOR_RAMP<T>::Two_Color_Ramp(0,1,OPENGL_COLOR::Gray(0,0),OPENGL_COLOR::Gray(1,1)));
        density_gradient_component->opengl_scalar_field.Update();
        Add_Component(density_gradient_component,"Density Gradient",'8',BASIC_VISUALIZATION::OWNED|BASIC_VISUALIZATION::START_HIDDEN);
        slice_manager.Add_Object(density_gradient_component);}
        
    opengl_world.Set_Key_Binding_Category("Soot");
    filename=basedir+"/%d/soot";
    if(has_valid_grid && FILE_UTILITIES::Frame_File_Exists(filename,start_frame)){
        OPENGL_COMPONENT_SCALAR_FIELD_3D<T>* soot_component=new OPENGL_COMPONENT_SCALAR_FIELD_3D<T>(grid,filename,OPENGL_COLOR_RAMP<T>::Two_Color_Ramp(0,.01,OPENGL_COLOR::Gray(0,0),OPENGL_COLOR::Gray(1,1)));
        soot_component->opengl_scalar_field.Update();
        Add_Component(soot_component,"Soot",'i',BASIC_VISUALIZATION::OWNED|BASIC_VISUALIZATION::START_HIDDEN);
        opengl_world.Append_Bind_Key('D',soot_component->Toggle_Color_Map_CB());
        opengl_world.Append_Bind_Key('`',soot_component->Toggle_Smooth_Slice_CB());
        slice_manager.Add_Object(soot_component);}

    filename=basedir+"/%d/soot_fuel";
    if(has_valid_grid && FILE_UTILITIES::Frame_File_Exists(filename,start_frame)){
        OPENGL_COMPONENT_SCALAR_FIELD_3D<T>* soot_fuel_component=new OPENGL_COMPONENT_SCALAR_FIELD_3D<T>(grid,filename,OPENGL_COLOR_RAMP<T>::Two_Color_Ramp(0,.01,OPENGL_COLOR::Gray(0,0),OPENGL_COLOR::Gray(1,1)));
        soot_fuel_component->opengl_scalar_field.Update();
        Add_Component(soot_fuel_component,"soot_fuel",'f',BASIC_VISUALIZATION::OWNED|BASIC_VISUALIZATION::START_HIDDEN);
        opengl_world.Append_Bind_Key('D',soot_fuel_component->Toggle_Color_Map_CB());
        opengl_world.Append_Bind_Key('`',soot_fuel_component->Toggle_Smooth_Slice_CB());
        slice_manager.Add_Object(soot_fuel_component);}

    opengl_world.Set_Key_Binding_Category("Internal Energy");
    filename=basedir+"/%d/internal_energy";
    if(has_valid_grid && FILE_UTILITIES::Frame_File_Exists(filename,start_frame)){
        OPENGL_COMPONENT_SCALAR_FIELD_3D<T>* internal_energy_component=new OPENGL_COMPONENT_SCALAR_FIELD_3D<T>(grid,filename,
            OPENGL_COLOR_RAMP<T>::Two_Color_Ramp(0,1,OPENGL_COLOR::Gray(0,0),OPENGL_COLOR::Gray(1,1)));
        internal_energy_component->opengl_scalar_field.Update();
        Add_Component(internal_energy_component,"Internal Energy",'E',BASIC_VISUALIZATION::OWNED|BASIC_VISUALIZATION::START_HIDDEN);
        slice_manager.Add_Object(internal_energy_component);}

#ifndef COMPILE_WITHOUT_DYADIC_SUPPORT
    filename=basedir+"/%d/octree_density";
    if(has_valid_octree_grid && FILE_UTILITIES::Frame_File_Exists(filename,start_frame)){
        OPENGL_COMPONENT_OCTREE_NODE_SCALAR_FIELD<T>* octree_density_component=new OPENGL_COMPONENT_OCTREE_NODE_SCALAR_FIELD<T>(octree_grid_component->opengl_grid.grid,filename,
            OPENGL_COLOR_RAMP<T>::Two_Color_Ramp(0,1,OPENGL_COLOR::Gray(0,0),OPENGL_COLOR::Gray(1,1)));
        Add_Component(octree_density_component,"Octree Density",'d',BASIC_VISUALIZATION::OWNED);
        opengl_world.Append_Bind_Key('+',octree_density_component->Increase_Point_Size_CB());
        opengl_world.Append_Bind_Key('_',octree_density_component->Decrease_Point_Size_CB());
        slice_manager.Add_Object(octree_density_component);}
#endif

    opengl_world.Set_Key_Binding_Category("Temperature");
    filename=basedir+"/%d/temperature";
    if(has_valid_grid && FILE_UTILITIES::Frame_File_Exists(filename,start_frame)){
        ARRAY<T,VECTOR<int,3> > temp0;FILE_UTILITIES::Read_From_File<RW>(FILE_UTILITIES::Get_Frame_Filename(filename,start_frame),temp0);
        T min_temp=temp0.Min(),max_temp=temp0.Max();min_temp=283.15;
        std::cout<<"Using temperature: [min="<<min_temp<<" max="<<max_temp<<"]"<<std::endl;
        OPENGL_COMPONENT_SCALAR_FIELD_3D<T>* temperature_component=new OPENGL_COMPONENT_SCALAR_FIELD_3D<T>(grid,filename,OPENGL_COLOR_RAMP<T>::Two_Color_Ramp(0,1,OPENGL_COLOR::Gray(0,0),OPENGL_COLOR::Red(1,1)));
        temperature_component->opengl_scalar_field.Set_Scale_Range(temp0.Min(),temp0.Max());
        Add_Component(temperature_component,"Temperature",'t',BASIC_VISUALIZATION::OWNED|BASIC_VISUALIZATION::START_HIDDEN);
        opengl_world.Append_Bind_Key('T',temperature_component->Toggle_Color_Map_CB());
        slice_manager.Add_Object(temperature_component);}

#ifndef COMPILE_WITHOUT_DYADIC_SUPPORT
    filename=basedir+"/%d/octree_temperature";
    if(has_valid_octree_grid && FILE_UTILITIES::Frame_File_Exists(filename,start_frame)){
        ARRAY<T> temp0;FILE_UTILITIES::Read_From_File<RW>(FILE_UTILITIES::Get_Frame_Filename(filename,start_frame),temp0);
        T min_temp=temp0.Min(),max_temp=temp0.Max();
        std::cout<<"Using octree temperature: [min="<<min_temp<<" max="<<max_temp<<"]"<<std::endl;
        if(max_temp-min_temp<1e-5){min_temp=283.15;max_temp=3000;std::cout<<"Switching to default temperature range 283.15 to 3000"<<std::endl;}
        OPENGL_COMPONENT_OCTREE_NODE_SCALAR_FIELD<T>* octree_temperature_component=new OPENGL_COMPONENT_OCTREE_NODE_SCALAR_FIELD<T>(octree_grid_component->opengl_grid.grid,filename,
            OPENGL_COLOR_RAMP<T>::Two_Color_Ramp(0,1,OPENGL_COLOR::Gray(0,0),OPENGL_COLOR::Red(1,1)));
        octree_temperature_component->Toggle_Draw();
        Add_Component(octree_temperature_component,"Octree Temperature",'t',BASIC_VISUALIZATION::OWNED|BASIC_VISUALIZATION::START_HIDDEN);
        opengl_world.Append_Bind_Key('+',octree_temperature_component->Increase_Point_Size_CB());
        opengl_world.Append_Bind_Key('_',octree_temperature_component->Decrease_Point_Size_CB());
        slice_manager.Add_Object(octree_temperature_component);}
#endif

    opengl_world.Set_Key_Binding_Category("SPH");
    filename=basedir+"/%d/sph_cell_weights";
    if(has_valid_grid && FILE_UTILITIES::Frame_File_Exists(filename,start_frame)){
        OPENGL_COMPONENT_SCALAR_FIELD_3D<T>* sph_cell_weight_component=new OPENGL_COMPONENT_SCALAR_FIELD_3D<T>(grid,filename,OPENGL_COLOR_RAMP<T>::Two_Color_Ramp(0,100,OPENGL_COLOR::Cyan(0,0),OPENGL_COLOR::Cyan(1)));
        sph_cell_weight_component->opengl_scalar_field.Update();
        Add_Component(sph_cell_weight_component,"SPH Cell Weights",'H',BASIC_VISUALIZATION::OWNED|BASIC_VISUALIZATION::START_HIDDEN);
        slice_manager.Add_Object(sph_cell_weight_component);}

    opengl_world.Set_Key_Binding_Category("Pressure");
    OPENGL_COLOR_MAP<T>* pressure_color_map=OPENGL_COLOR_RAMP<T>::Two_Color_Ramp(101635,1e7,OPENGL_COLOR::Cyan(0,0),OPENGL_COLOR::Cyan(1));
    filename=basedir+"/%d/pressure";
    if(has_valid_grid && FILE_UTILITIES::Frame_File_Exists(filename,start_frame)){
        OPENGL_COMPONENT_SCALAR_FIELD_3D<T>* pressure_component=new OPENGL_COMPONENT_SCALAR_FIELD_3D<T>(mac_grid,filename,pressure_color_map);
        Add_Component(pressure_component,"Pressure",'7',BASIC_VISUALIZATION::OWNED|BASIC_VISUALIZATION::START_HIDDEN);
        opengl_world.Append_Bind_Key('D',pressure_component->Toggle_Color_Map_CB());
        opengl_world.Append_Bind_Key('`',pressure_component->Toggle_Smooth_Slice_CB());
        slice_manager.Add_Object(pressure_component);}

    filename=basedir+"/%d/pressure_gradient";
    if(has_valid_grid && FILE_UTILITIES::Frame_File_Exists(filename,start_frame)){
        OPENGL_COMPONENT_SCALAR_FIELD_3D<T>* pressure_gradient_component=new OPENGL_COMPONENT_SCALAR_FIELD_3D<T>(grid,filename,
            OPENGL_COLOR_RAMP<T>::Two_Color_Ramp(0,1000,OPENGL_COLOR::Magenta(0,0),OPENGL_COLOR::Magenta(1)));
        pressure_gradient_component->opengl_scalar_field.Update();
        Add_Component(pressure_gradient_component,"Pressure Gradient",'1',BASIC_VISUALIZATION::OWNED|BASIC_VISUALIZATION::START_HIDDEN);
        slice_manager.Add_Object(pressure_gradient_component);}

#ifndef COMPILE_WITHOUT_RLE_SUPPORT
    filename=basedir+"/%d/rle_pressure";
    if(has_valid_rle_grid && FILE_UTILITIES::Frame_File_Exists(filename,start_frame)){
        OPENGL_COMPONENT_RLE_CELL_SCALAR_FIELD_3D<T>* rle_pressure_component=new OPENGL_COMPONENT_RLE_CELL_SCALAR_FIELD_3D<T>(rle_grid_component->opengl_grid.grid,filename,pressure_color_map);
        Add_Component(rle_pressure_component,"RLE Pressure",'7',BASIC_VISUALIZATION::OWNED|BASIC_VISUALIZATION::START_HIDDEN);
        slice_manager.Add_Object(rle_pressure_component);}
#endif

#ifndef COMPILE_WITHOUT_DYADIC_SUPPORT
    filename=basedir+"/%d/octree_pressure";
    if(has_valid_octree_grid && FILE_UTILITIES::Frame_File_Exists(filename,start_frame)){
        OPENGL_COMPONENT_OCTREE_CELL_SCALAR_FIELD<T>* octree_pressure_component=new OPENGL_COMPONENT_OCTREE_CELL_SCALAR_FIELD<T>(octree_grid_component->opengl_grid.grid,filename,pressure_color_map);
        Add_Component(octree_pressure_component,"Octree Pressure",'7',BASIC_VISUALIZATION::OWNED|BASIC_VISUALIZATION::START_HIDDEN);
        slice_manager.Add_Object(octree_pressure_component);}

    filename=basedir+"/%d/octree_occupied_cells";
    if(has_valid_octree_grid && FILE_UTILITIES::Frame_File_Exists(filename,start_frame)){
        OPENGL_COMPONENT_OCTREE_CELL_SCALAR_FIELD<T,bool>* octree_occupied_cells=new OPENGL_COMPONENT_OCTREE_CELL_SCALAR_FIELD<T,bool>(octree_grid_component->opengl_grid.grid,filename,
            new OPENGL_CONSTANT_COLOR_MAP<bool>(OPENGL_COLOR::Green()));
        Add_Component(octree_occupied_cells,"Octree Occupied Cells",'\0',BASIC_VISUALIZATION::OWNED|BASIC_VISUALIZATION::START_HIDDEN);
        opengl_world.Append_Bind_Key(OPENGL_KEY(OPENGL_KEY::F11),octree_occupied_cells->Toggle_Draw_CB());
        slice_manager.Add_Object(octree_occupied_cells);}
#endif

    opengl_world.Set_Key_Binding_Category("Velocity");
    { // regular and mac velocity TODO: kill node based stuff?
        OPENGL_COMPONENT_GRID_BASED_VECTOR_FIELD_3D<T>* vector_velocity_component=0;
        OPENGL_COMPONENT_MAC_VELOCITY_FIELD_3D<T>* mac_velocity_component=0;
        filename=basedir+"/%d/velocities";
        if(has_valid_grid && FILE_UTILITIES::Frame_File_Exists(filename,start_frame)){
            if(node_based){
                vector_velocity_component=new OPENGL_COMPONENT_GRID_BASED_VECTOR_FIELD_3D<T>(regular_grid,filename);
                vector_velocity_component->opengl_grid_based_vector_field.size=0.01;
                vector_velocity_component->opengl_grid_based_vector_field.vector_color=OPENGL_COLOR::Green();
                Add_Component(vector_velocity_component,"Node velocities",'\0',BASIC_VISUALIZATION::OWNED|BASIC_VISUALIZATION::START_HIDDEN);
                slice_manager.Add_Object(vector_velocity_component);}
            else{
                mac_velocity_component=new OPENGL_COMPONENT_MAC_VELOCITY_FIELD_3D<T>(mac_grid,filename);
                mac_velocity_component->opengl_mac_velocity_field.size=0.1;
                mac_velocity_component->opengl_mac_velocity_field.vector_color=OPENGL_COLOR::Green();
                Add_Component(mac_velocity_component,"MAC velocities",'\0',BASIC_VISUALIZATION::OWNED|BASIC_VISUALIZATION::START_HIDDEN);
                slice_manager.Add_Object(mac_velocity_component);}}
        
        filename=basedir+"/%d/mac_velocities";
        if(has_valid_grid && (always_add_mac_velocities || FILE_UTILITIES::Frame_File_Exists(filename,start_frame))){
            mac_velocity_component=new OPENGL_COMPONENT_MAC_VELOCITY_FIELD_3D<T>(mac_grid,filename);
            mac_velocity_component->opengl_mac_velocity_field.size=0.01;
            mac_velocity_component->opengl_mac_velocity_field.vector_color=OPENGL_COLOR::Magenta();
            mac_velocity_component->opengl_mac_velocity_field.Set_Velocity_Mode(OPENGL_MAC_VELOCITY_FIELD_3D<T>::FACE_CENTERED);
            Add_Component(mac_velocity_component,"MAC velocities",'\0',BASIC_VISUALIZATION::OWNED|BASIC_VISUALIZATION::START_HIDDEN);
            slice_manager.Add_Object(mac_velocity_component);}
        
        if(mac_velocity_component){
            opengl_world.Append_Bind_Key('C',mac_velocity_component->Toggle_Draw_CB());
            opengl_world.Append_Bind_Key('V',mac_velocity_component->Toggle_Velocity_Mode_And_Draw_CB());
            opengl_world.Append_Bind_Key('v',mac_velocity_component->Toggle_Draw_Vorticity_CB());
            opengl_world.Append_Bind_Key('N',mac_velocity_component->Normalize_Vorticity_Color_Map_CB());
            opengl_world.Append_Bind_Key('=',mac_velocity_component->Increase_Vector_Size_CB());
            opengl_world.Append_Bind_Key('-',mac_velocity_component->Decrease_Vector_Size_CB());
            opengl_world.Append_Bind_Key('h',mac_velocity_component->Toggle_Arrowhead_CB());}
        if(vector_velocity_component){
            opengl_world.Append_Bind_Key('v',vector_velocity_component->Toggle_Draw_CB());
            opengl_world.Append_Bind_Key('=',vector_velocity_component->Increase_Vector_Size_CB());
            opengl_world.Append_Bind_Key('-',vector_velocity_component->Decrease_Vector_Size_CB());
            opengl_world.Append_Bind_Key('h',vector_velocity_component->Toggle_Arrowhead_CB());}

        OPENGL_COMPONENT_MAC_VELOCITY_FIELD_3D<T>* coarse_mac_velocity_component=0;
        opengl_world.Set_Key_Binding_Category("Coarse Velocity");
        filename=basedir+"/%d/coarse_mac_velocities";
        if(has_valid_coarse_grid && FILE_UTILITIES::Frame_File_Exists(filename,start_frame)){
            coarse_mac_velocity_component=new OPENGL_COMPONENT_MAC_VELOCITY_FIELD_3D<T>(coarse_mac_grid,filename);
            coarse_mac_velocity_component->opengl_mac_velocity_field.scale=grid.Domain_Indices().max_corner(1)/coarse_grid.Domain_Indices().max_corner(1);
            coarse_mac_velocity_component->opengl_mac_velocity_field.size=.01;
            coarse_mac_velocity_component->opengl_mac_velocity_field.vector_color=OPENGL_COLOR::Blue();
            coarse_mac_velocity_component->opengl_mac_velocity_field.Set_Velocity_Mode(OPENGL_MAC_VELOCITY_FIELD_3D<T>::FACE_CENTERED);
            Add_Component(coarse_mac_velocity_component,"coarse MAC velocities",'\0',BASIC_VISUALIZATION::OWNED|BASIC_VISUALIZATION::START_HIDDEN);
            slice_manager.Add_Object(coarse_mac_velocity_component);}
        if(coarse_mac_velocity_component){
            opengl_world.Append_Bind_Key('x',coarse_mac_velocity_component->Toggle_Draw_CB());
            opengl_world.Append_Bind_Key('V',coarse_mac_velocity_component->Toggle_Velocity_Mode_And_Draw_CB());
            opengl_world.Append_Bind_Key('=',coarse_mac_velocity_component->Increase_Vector_Size_CB());
            opengl_world.Append_Bind_Key('-',coarse_mac_velocity_component->Decrease_Vector_Size_CB());
            opengl_world.Append_Bind_Key('h',coarse_mac_velocity_component->Toggle_Arrowhead_CB());}
        filename=basedir+"/%d/centered_velocities";
        if(has_valid_grid && FILE_UTILITIES::Frame_File_Exists(filename,start_frame)){
            OPENGL_COMPONENT_GRID_BASED_VECTOR_FIELD_3D<T>* center_velocity_component=new OPENGL_COMPONENT_GRID_BASED_VECTOR_FIELD_3D<T>(mac_grid,filename);
            center_velocity_component->opengl_grid_based_vector_field.size=0.1;
            center_velocity_component->opengl_grid_based_vector_field.vector_color=OPENGL_COLOR::Green();
            Add_Component(center_velocity_component,"Centered velocities",'B',BASIC_VISUALIZATION::OWNED|BASIC_VISUALIZATION::START_HIDDEN);
            opengl_world.Append_Bind_Key('=',center_velocity_component->Increase_Vector_Size_CB());
            opengl_world.Append_Bind_Key('-',center_velocity_component->Decrease_Vector_Size_CB());
            slice_manager.Add_Object(center_velocity_component);}
    }

    filename=basedir+"/%d/object_velocities";
    if(has_valid_grid && FILE_UTILITIES::Frame_File_Exists(filename,start_frame)){
        OPENGL_COMPONENT_GRID_BASED_VECTOR_FIELD_3D<T>* object_velocity_component=new OPENGL_COMPONENT_GRID_BASED_VECTOR_FIELD_3D<T>(regular_grid,filename);
        object_velocity_component->opengl_grid_based_vector_field.size=0.01;
        object_velocity_component->opengl_grid_based_vector_field.vector_color=OPENGL_COLOR::Yellow();
        Add_Component(object_velocity_component,"Object velocities",'%',BASIC_VISUALIZATION::OWNED|BASIC_VISUALIZATION::START_HIDDEN);
        opengl_world.Append_Bind_Key('=',object_velocity_component->Increase_Vector_Size_CB());
        opengl_world.Append_Bind_Key('-',object_velocity_component->Decrease_Vector_Size_CB());
        opengl_world.Append_Bind_Key('h',object_velocity_component->Toggle_Arrowhead_CB());
        slice_manager.Add_Object(object_velocity_component);}

    filename=basedir+"/%d/forces";
    if(has_valid_grid && FILE_UTILITIES::Frame_File_Exists(filename,start_frame)){
        OPENGL_COMPONENT_MAC_VELOCITY_FIELD_3D<T>* force_component=new OPENGL_COMPONENT_MAC_VELOCITY_FIELD_3D<T>(mac_grid,filename);
        force_component->opengl_mac_velocity_field.size=0.01;
        force_component->opengl_mac_velocity_field.vector_color=OPENGL_COLOR::Yellow();
        force_component->opengl_mac_velocity_field.Set_Velocity_Mode(OPENGL_MAC_VELOCITY_FIELD_3D<T>::FACE_CENTERED);
        Add_Component(force_component,"fluid control force",'G',BASIC_VISUALIZATION::OWNED|BASIC_VISUALIZATION::START_HIDDEN);
        opengl_world.Append_Bind_Key('=',force_component->Increase_Vector_Size_CB());
        opengl_world.Append_Bind_Key('-',force_component->Decrease_Vector_Size_CB());
        opengl_world.Append_Bind_Key('h',force_component->Toggle_Arrowhead_CB());
        slice_manager.Add_Object(force_component);}

#ifndef COMPILE_WITHOUT_RLE_SUPPORT
    opengl_world.Set_Key_Binding_Category("RLE Velocity");
    filename=basedir+"/%d/rle_face_velocities";
    if(has_valid_rle_grid && FILE_UTILITIES::Frame_File_Exists(filename,start_frame)){
        OPENGL_COMPONENT_RLE_FACE_SCALAR_FIELD_3D<T,T>* rle_face_velocity_component=new OPENGL_COMPONENT_RLE_FACE_SCALAR_FIELD_3D<T,T>(rle_grid_component->opengl_grid.grid,filename,
            new OPENGL_CONSTANT_COLOR_MAP<T>(OPENGL_COLOR::Green()),false);
        rle_face_velocity_component->opengl_rle_face_scalar_field.line_size=0.01;
        Add_Component(rle_face_velocity_component,"RLE Face Velocities",'V',BASIC_VISUALIZATION::OWNED|BASIC_VISUALIZATION::START_HIDDEN);
        opengl_world.Append_Bind_Key('C',rle_face_velocity_component->Toggle_Draw_CB());
        opengl_world.Append_Bind_Key('=',rle_face_velocity_component->Increase_Vector_Size_CB());
        opengl_world.Append_Bind_Key('-',rle_face_velocity_component->Decrease_Vector_Size_CB());
        slice_manager.Add_Object(rle_face_velocity_component);}
#endif

#ifndef COMPILE_WITHOUT_DYADIC_SUPPORT
    filename=basedir+"/%d/octree_face_velocities";
    if(has_valid_octree_grid && FILE_UTILITIES::Frame_File_Exists(filename,start_frame)){
        OPENGL_COMPONENT_OCTREE_FACE_SCALAR_FIELD<T,T>* octree_face_velocity_component=new OPENGL_COMPONENT_OCTREE_FACE_SCALAR_FIELD<T,T>(octree_grid_component->opengl_grid.grid,filename,
            new OPENGL_CONSTANT_COLOR_MAP<T>(OPENGL_COLOR::Magenta()),false);
        octree_face_velocity_component->opengl_octree_face_scalar_field.line_size=0.01;
        Add_Component(octree_face_velocity_component,"Octree Face Velocities",'V',BASIC_VISUALIZATION::OWNED|BASIC_VISUALIZATION::START_HIDDEN);
        opengl_world.Append_Bind_Key('=',octree_face_velocity_component->Increase_Vector_Size_CB());
        opengl_world.Append_Bind_Key('-',octree_face_velocity_component->Decrease_Vector_Size_CB());
        slice_manager.Add_Object(octree_face_velocity_component);}
#endif

    filename=basedir+"/%d/velocities_ghost_fuel";
    if(has_valid_grid && FILE_UTILITIES::Frame_File_Exists(filename,start_frame)){
        if(node_based){
            OPENGL_COMPONENT_GRID_BASED_VECTOR_FIELD_3D<T>* vector_velocity_ghost_plus_component=new OPENGL_COMPONENT_GRID_BASED_VECTOR_FIELD_3D<T>(regular_grid,filename);
            vector_velocity_ghost_plus_component->opengl_grid_based_vector_field.size=0.01;
            vector_velocity_ghost_plus_component->opengl_grid_based_vector_field.vector_color=OPENGL_COLOR::Green();
            Add_Component(vector_velocity_ghost_plus_component,"Node ghost plus velocities",'b',BASIC_VISUALIZATION::OWNED|BASIC_VISUALIZATION::START_HIDDEN);
            opengl_world.Append_Bind_Key('h',vector_velocity_ghost_plus_component->Toggle_Arrowhead_CB());
            opengl_world.Append_Bind_Key('=',vector_velocity_ghost_plus_component->Increase_Vector_Size_CB());
            opengl_world.Append_Bind_Key('-',vector_velocity_ghost_plus_component->Decrease_Vector_Size_CB());
            slice_manager.Add_Object(vector_velocity_ghost_plus_component);}}

    filename=basedir+"/%d/velocities_ghost";
    if(has_valid_grid && FILE_UTILITIES::Frame_File_Exists(filename,start_frame)){
        if(node_based){
            OPENGL_COMPONENT_GRID_BASED_VECTOR_FIELD_3D<T>* vector_velocity_ghost_minus_component=new OPENGL_COMPONENT_GRID_BASED_VECTOR_FIELD_3D<T>(regular_grid,filename);
            vector_velocity_ghost_minus_component->opengl_grid_based_vector_field.size=0.01;
            vector_velocity_ghost_minus_component->opengl_grid_based_vector_field.vector_color=OPENGL_COLOR::Green();
            vector_velocity_ghost_minus_component->Set_Draw(false);
            Add_Component(vector_velocity_ghost_minus_component,"Node ghost minus velocities",'c',BASIC_VISUALIZATION::OWNED|BASIC_VISUALIZATION::START_HIDDEN);
            opengl_world.Append_Bind_Key('h',vector_velocity_ghost_minus_component->Toggle_Arrowhead_CB());
            opengl_world.Append_Bind_Key('=',vector_velocity_ghost_minus_component->Increase_Vector_Size_CB());
            opengl_world.Append_Bind_Key('-',vector_velocity_ghost_minus_component->Decrease_Vector_Size_CB());
            slice_manager.Add_Object(vector_velocity_ghost_minus_component);}}

    filename=basedir+"/%d/center_velocities";
    if(has_valid_grid && FILE_UTILITIES::Frame_File_Exists(filename,start_frame) && parameter_list.template Get_Parameter<bool>("use_cell_centered_velocities",false)){
        OPENGL_COMPONENT_GRID_BASED_VECTOR_FIELD_3D<T>* center_velocity_component=new OPENGL_COMPONENT_GRID_BASED_VECTOR_FIELD_3D<T>(mac_grid,filename);
        center_velocity_component->opengl_grid_based_vector_field.size=0.1;
        center_velocity_component->opengl_grid_based_vector_field.vector_color=OPENGL_COLOR::Green();
        Add_Component(center_velocity_component,"Centered velocities",'V',BASIC_VISUALIZATION::OWNED|BASIC_VISUALIZATION::START_HIDDEN);
        opengl_world.Append_Bind_Key('=',center_velocity_component->Increase_Vector_Size_CB());
        opengl_world.Append_Bind_Key('-',center_velocity_component->Decrease_Vector_Size_CB());
        slice_manager.Add_Object(center_velocity_component);}

    opengl_world.Set_Key_Binding_Category("Pressure Jump");
    filename=basedir+"/%d/pressure_jumps";
    if(has_valid_grid && FILE_UTILITIES::Frame_File_Exists(filename,start_frame)){
        OPENGL_COMPONENT_GRID_BASED_VECTOR_FIELD_3D<T>* pressure_jump_component=new OPENGL_COMPONENT_GRID_BASED_VECTOR_FIELD_3D<T>(grid,filename);
        pressure_jump_component->opengl_grid_based_vector_field.size=0.001;
        pressure_jump_component->opengl_grid_based_vector_field.vector_color=OPENGL_COLOR::Magenta();
        Add_Component(pressure_jump_component,"Pressure jumps",'&',BASIC_VISUALIZATION::OWNED|BASIC_VISUALIZATION::START_HIDDEN);
        opengl_world.Append_Bind_Key('=',pressure_jump_component->Increase_Vector_Size_CB());
        opengl_world.Append_Bind_Key('-',pressure_jump_component->Decrease_Vector_Size_CB());
        opengl_world.Append_Bind_Key('h',pressure_jump_component->Toggle_Arrowhead_CB());
        slice_manager.Add_Object(pressure_jump_component);}

    filename=basedir+"/%d/beta_face";
    if(has_valid_grid && FILE_UTILITIES::Frame_File_Exists(basedir+"/%d/beta_face",start_frame)){
        OPENGL_COMPONENT_FACE_SCALAR_FIELD_3D<T,T>* beta_face_component=new OPENGL_COMPONENT_FACE_SCALAR_FIELD_3D<T,T>(mac_grid,filename,OPENGL_COLOR_RAMP<T>::Two_Color_Ramp(0,(T).002,OPENGL_COLOR::Gray(1),OPENGL_COLOR::Gray(0)));
        Add_Component(beta_face_component,"Beta Face",'\0',BASIC_VISUALIZATION::OWNED|BASIC_VISUALIZATION::START_HIDDEN);
        opengl_world.Append_Bind_Key(OPENGL_KEY(OPENGL_KEY::F6),beta_face_component->Toggle_Draw_CB());
        slice_manager.Add_Object(beta_face_component);}

    opengl_world.Set_Key_Binding_Category("Pressure Jump");
#ifndef COMPILE_WITHOUT_DYADIC_SUPPORT
    filename=basedir+"/%d/octree_pressure_jumps";
    if(has_valid_octree_grid && FILE_UTILITIES::Frame_File_Exists(filename,start_frame)){
        OPENGL_COMPONENT_OCTREE_NODE_VECTOR_FIELD<T>* octree_pressure_jump_component=new OPENGL_COMPONENT_OCTREE_NODE_VECTOR_FIELD<T>(octree_grid_component->opengl_grid.grid,filename);
        octree_pressure_jump_component->opengl_octree_node_vector_field.size=0.01;
        octree_pressure_jump_component->opengl_octree_node_vector_field.vector_color=OPENGL_COLOR::Magenta();
        Add_Component(octree_pressure_jump_component,"Octree pressure jumps",'&',BASIC_VISUALIZATION::OWNED|BASIC_VISUALIZATION::START_HIDDEN);
        opengl_world.Append_Bind_Key('=',octree_pressure_jump_component->Increase_Vector_Size_CB());
        opengl_world.Append_Bind_Key('-',octree_pressure_jump_component->Decrease_Vector_Size_CB());
        opengl_world.Append_Bind_Key('h',octree_pressure_jump_component->Toggle_Arrowhead_CB());
        slice_manager.Add_Object(octree_pressure_jump_component);}
#endif

    filename=basedir+"/%d/pseudo_dirichlet";
    if(has_valid_grid && FILE_UTILITIES::Frame_File_Exists(filename,start_frame)){
        OPENGL_COMPONENT_PSEUDO_DIRICHLET_3D<T>* pseudo_dirichlet_component=new OPENGL_COMPONENT_PSEUDO_DIRICHLET_3D<T>(grid,filename);
        pseudo_dirichlet_component->Set_Vector_Size(0.1);
        Add_Component(pseudo_dirichlet_component,"pseudo dirichlet",'\0',BASIC_VISUALIZATION::OWNED|BASIC_VISUALIZATION::START_HIDDEN);
        opengl_world.Append_Bind_Key(OPENGL_KEY(OPENGL_KEY::F3),pseudo_dirichlet_component->Toggle_Draw_CB());
        opengl_world.Append_Bind_Key('=',pseudo_dirichlet_component->Increase_Vector_Size_CB());
        opengl_world.Append_Bind_Key('-',pseudo_dirichlet_component->Decrease_Vector_Size_CB());
        slice_manager.Add_Object(pseudo_dirichlet_component);}

    filename=basedir+"/%d/thin_shells_grid_visibility";
    if(has_valid_grid && FILE_UTILITIES::Frame_File_Exists(filename,start_frame)){
        OPENGL_COMPONENT_THIN_SHELLS_DEBUGGING_3D<T>* thin_shells_debugging_component=new OPENGL_COMPONENT_THIN_SHELLS_DEBUGGING_3D<T>(grid,basedir);
        Add_Component(thin_shells_debugging_component,"thin shells debugging",'\0',BASIC_VISUALIZATION::OWNED);
        opengl_world.Append_Bind_Key(OPENGL_KEY(OPENGL_KEY::F4),thin_shells_debugging_component->Toggle_Draw_Grid_Visibility_Mode_CB());
        opengl_world.Append_Bind_Key(OPENGL_KEY(OPENGL_KEY::F5),thin_shells_debugging_component->Toggle_Draw_Density_Valid_Mask_CB());
        slice_manager.Add_Object(thin_shells_debugging_component);}

    if(has_valid_grid){
        OPENGL_GRID_3D<T>* opengl_grid=new OPENGL_GRID_3D<T>(*(new GRID<TV>(grid)),OPENGL_COLOR::Gray(0.5));
        grid_component=new OPENGL_COMPONENT_BASIC<OPENGL_GRID_3D<T> >(*opengl_grid);
        opengl_world.Set_Key_Binding_Category("Grid");
        Add_Component(grid_component,"Grid",'6',BASIC_VISUALIZATION::START_HIDDEN|BASIC_VISUALIZATION::SELECTABLE);
        opengl_world.Append_Bind_Key('^',grid_component->object.Toggle_Draw_Ghost_Values_CB());
        slice_manager.Add_Object(grid_component);}
    if(has_valid_coarse_grid){
        OPENGL_GRID_3D<T>* opengl_grid=new OPENGL_GRID_3D<T>(*(new GRID<TV>(coarse_grid)),OPENGL_COLOR::Ground_Tan(.5));opengl_grid->scale=grid.Domain_Indices().max_corner(1)/coarse_grid.Domain_Indices().max_corner(1);
        coarse_grid_component=new OPENGL_COMPONENT_BASIC<OPENGL_GRID_3D<T> >(*opengl_grid);
        Add_Component(coarse_grid_component,"Coarse Grid",'y',BASIC_VISUALIZATION::START_HIDDEN|BASIC_VISUALIZATION::SELECTABLE);
        opengl_world.Append_Bind_Key('^',coarse_grid_component->object.Toggle_Draw_Ghost_Values_CB());
        slice_manager.Add_Object(coarse_grid_component);}


#ifndef COMPILE_WITHOUT_RLE_SUPPORT
    if(rle_grid_component){
        opengl_world.Set_Key_Binding_Category("Grid");
        Add_Component(rle_grid_component,"RLE Grid",'6',BASIC_VISUALIZATION::START_HIDDEN|BASIC_VISUALIZATION::SELECTABLE);}
#endif
#ifndef COMPILE_WITHOUT_DYADIC_SUPPORT
    if(octree_grid_component){
        opengl_world.Set_Key_Binding_Category("Grid");
        Add_Component(octree_grid_component,"Octree Grid",'6',BASIC_VISUALIZATION::START_HIDDEN|BASIC_VISUALIZATION::SELECTABLE);}
#endif

    opengl_world.Set_Key_Binding_Category("Particles");

    bool particles_stored_per_cell_uniform=false;
#ifndef COMPILE_WITHOUT_DYADIC_SUPPORT
    if((
#ifndef COMPILE_WITHOUT_RLE_SUPPORT
           has_valid_rle_grid ||
#endif
           has_valid_octree_grid) && !has_valid_grid) particles_stored_per_cell_adaptive=true;
    else
#endif
        if(has_valid_grid) particles_stored_per_cell_uniform=true;


    OPENGL_COLOR_MAP<bool>* psi_N_color_map=new OPENGL_CONSTANT_COLOR_MAP<bool>(OPENGL_COLOR::Cyan());
    OPENGL_COLOR_MAP<bool>* psi_D_color_map=new OPENGL_CONSTANT_COLOR_MAP<bool>(OPENGL_COLOR::Magenta());
    if(has_valid_coarse_grid){
        // TODO: this is legacy output form, probably can be removed now
        OPENGL_COMPONENT_FACE_SCALAR_FIELD_3D<T,bool>* psi_N_component=0;
        if(FILE_UTILITIES::Frame_File_Exists(basedir+"/%d/psi_N_u",start_frame) &&
           FILE_UTILITIES::Frame_File_Exists(basedir+"/%d/psi_N_v",start_frame) &&
           FILE_UTILITIES::Frame_File_Exists(basedir+"/%d/psi_N_w",start_frame))
            psi_N_component=new OPENGL_COMPONENT_FACE_SCALAR_FIELD_3D<T,bool>(coarse_grid,basedir+"/%d/psi_N_u",basedir+"/%d/psi_N_v",basedir+"/%d/psi_N_w",psi_N_color_map);
        else if(FILE_UTILITIES::Frame_File_Exists(basedir+"/%d/coarse_psi_N",start_frame))
            psi_N_component=new OPENGL_COMPONENT_FACE_SCALAR_FIELD_3D<T,bool>(coarse_grid,basedir+"/%d/coarse_psi_N",psi_N_color_map);
        if(psi_N_component){
            psi_N_component->opengl_scalar_field.scale=grid.Domain_Indices().max_corner(1)/coarse_grid.Domain_Indices().max_corner(1);
            Add_Component(psi_N_component,"Coarse Psi_N points",'\0',BASIC_VISUALIZATION::OWNED|BASIC_VISUALIZATION::START_HIDDEN);
            opengl_world.Append_Bind_Key(OPENGL_KEY(OPENGL_KEY::F2),psi_N_component->Toggle_Draw_CB());
            slice_manager.Add_Object(psi_N_component);}}
    if(has_valid_grid){
        // TODO: this is legacy output form, probably can be removed now
        OPENGL_COMPONENT_FACE_SCALAR_FIELD_3D<T,bool>* psi_N_component=0;
        if(FILE_UTILITIES::Frame_File_Exists(basedir+"/%d/psi_N_u",start_frame) &&
           FILE_UTILITIES::Frame_File_Exists(basedir+"/%d/psi_N_v",start_frame) &&
           FILE_UTILITIES::Frame_File_Exists(basedir+"/%d/psi_N_w",start_frame))
            psi_N_component=new OPENGL_COMPONENT_FACE_SCALAR_FIELD_3D<T,bool>(grid,basedir+"/%d/psi_N_u",basedir+"/%d/psi_N_v",basedir+"/%d/psi_N_w",psi_N_color_map);
        else if(FILE_UTILITIES::Frame_File_Exists(basedir+"/%d/psi_N",start_frame))
            psi_N_component=new OPENGL_COMPONENT_FACE_SCALAR_FIELD_3D<T,bool>(grid,basedir+"/%d/psi_N",psi_N_color_map);
        if(psi_N_component){
            Add_Component(psi_N_component,"Psi_N points",'\0',BASIC_VISUALIZATION::OWNED|BASIC_VISUALIZATION::START_HIDDEN);
            opengl_world.Append_Bind_Key(OPENGL_KEY(OPENGL_KEY::F1),psi_N_component->Toggle_Draw_CB());
            slice_manager.Add_Object(psi_N_component);}}

#ifndef COMPILE_WITHOUT_RLE_SUPPORT
    filename=basedir+"/%d/rle_psi_N";
    if(has_valid_rle_grid && FILE_UTILITIES::Frame_File_Exists(filename,start_frame)){
        OPENGL_COMPONENT_RLE_FACE_SCALAR_FIELD_3D<T,bool>* rle_psi_N_component=new OPENGL_COMPONENT_RLE_FACE_SCALAR_FIELD_3D<T,bool>(rle_grid_component->opengl_grid.grid,filename,psi_N_color_map,true);
        Add_Component(rle_psi_N_component,"RLE Psi_N points",'\0',BASIC_VISUALIZATION::OWNED|BASIC_VISUALIZATION::START_HIDDEN);
        opengl_world.Append_Bind_Key(OPENGL_KEY(OPENGL_KEY::F1),rle_psi_N_component->Toggle_Draw_CB());
        slice_manager.Add_Object(rle_psi_N_component);}
#endif

#ifndef COMPILE_WITHOUT_DYADIC_SUPPORT
    filename=basedir+"/%d/octree_psi_N";
    if(has_valid_octree_grid && FILE_UTILITIES::Frame_File_Exists(filename,start_frame)){
        OPENGL_COMPONENT_OCTREE_FACE_SCALAR_FIELD<T,bool>* octree_psi_N_component=new OPENGL_COMPONENT_OCTREE_FACE_SCALAR_FIELD<T,bool>(octree_grid_component->opengl_grid.grid,filename,psi_N_color_map,true);
        Add_Component(octree_psi_N_component,"Octree Psi_N points",'\0',BASIC_VISUALIZATION::OWNED|BASIC_VISUALIZATION::START_HIDDEN);
        opengl_world.Append_Bind_Key(OPENGL_KEY(OPENGL_KEY::F1),octree_psi_N_component->Toggle_Draw_CB());
        slice_manager.Add_Object(octree_psi_N_component);}
#endif

    filename=basedir+"/%d/coarse_psi_D";
    if(has_valid_coarse_grid && FILE_UTILITIES::Frame_File_Exists(filename,start_frame)){
        OPENGL_COMPONENT_SCALAR_FIELD_3D<T,bool>* psi_D_component=new OPENGL_COMPONENT_SCALAR_FIELD_3D<T,bool>(coarse_mac_grid,filename,psi_D_color_map,OPENGL_SCALAR_FIELD_3D<T,bool>::DRAW_POINTS);
        Add_Component(psi_D_component,"Psi_D points",'\0',BASIC_VISUALIZATION::OWNED|BASIC_VISUALIZATION::START_HIDDEN);
        opengl_world.Append_Bind_Key(OPENGL_KEY(OPENGL_KEY::F2),psi_D_component->Toggle_Draw_CB());
        slice_manager.Add_Object(psi_D_component);}
    filename=basedir+"/%d/psi_D";
    if(has_valid_grid && FILE_UTILITIES::Frame_File_Exists(filename,start_frame)){
        OPENGL_COMPONENT_SCALAR_FIELD_3D<T,bool>* psi_D_component=new OPENGL_COMPONENT_SCALAR_FIELD_3D<T,bool>(mac_grid,filename,psi_D_color_map,OPENGL_SCALAR_FIELD_3D<T,bool>::DRAW_POINTS);
        Add_Component(psi_D_component,"Psi_D points",'\0',BASIC_VISUALIZATION::OWNED|BASIC_VISUALIZATION::START_HIDDEN);
        opengl_world.Append_Bind_Key(OPENGL_KEY(OPENGL_KEY::F1),psi_D_component->Toggle_Draw_CB());
        slice_manager.Add_Object(psi_D_component);}

    filename=basedir+"/%d/maccormack_cell_mask";
    if(has_valid_grid && FILE_UTILITIES::Frame_File_Exists(filename,start_frame)){
        OPENGL_COMPONENT_SCALAR_FIELD_3D<T,bool>* maccormack_cell_mask_component=new OPENGL_COMPONENT_SCALAR_FIELD_3D<T,bool>(mac_grid,filename,psi_D_color_map,OPENGL_SCALAR_FIELD_3D<T,bool>::DRAW_POINTS);
        Add_Component(maccormack_cell_mask_component,"Maccormack cell mask points",'\0',BASIC_VISUALIZATION::OWNED|BASIC_VISUALIZATION::START_HIDDEN);
        opengl_world.Append_Bind_Key(OPENGL_KEY(OPENGL_KEY::F9),maccormack_cell_mask_component->Toggle_Draw_CB());
        slice_manager.Add_Object(maccormack_cell_mask_component);}

    filename=basedir+"/%d/maccormack_face_mask";
    if(has_valid_grid && FILE_UTILITIES::Frame_File_Exists(filename,start_frame)){
        OPENGL_COMPONENT_FACE_SCALAR_FIELD_3D<T,bool>* maccormack_face_mask_component=new OPENGL_COMPONENT_FACE_SCALAR_FIELD_3D<T,bool>(grid,filename,psi_N_color_map);
        Add_Component(maccormack_face_mask_component,"Maccormack face mask points",'\0',BASIC_VISUALIZATION::OWNED|BASIC_VISUALIZATION::START_HIDDEN);
        opengl_world.Append_Bind_Key(OPENGL_KEY(OPENGL_KEY::F10),maccormack_face_mask_component->Toggle_Draw_CB());
        slice_manager.Add_Object(maccormack_face_mask_component);}

#ifndef COMPILE_WITHOUT_RLE_SUPPORT
    filename=basedir+"/%d/rle_psi_D";
    if(has_valid_rle_grid && FILE_UTILITIES::Frame_File_Exists(filename,start_frame)){
        OPENGL_COMPONENT_RLE_CELL_SCALAR_FIELD_3D<T,bool>* rle_psi_D_component=new OPENGL_COMPONENT_RLE_CELL_SCALAR_FIELD_3D<T,bool>(rle_grid_component->opengl_grid.grid,filename,psi_D_color_map);
        Add_Component(rle_psi_D_component,"RLE Psi_D points",'\0',BASIC_VISUALIZATION::OWNED|BASIC_VISUALIZATION::START_HIDDEN);
        opengl_world.Append_Bind_Key(OPENGL_KEY(OPENGL_KEY::F1),rle_psi_D_component->Toggle_Draw_CB());
        slice_manager.Add_Object(rle_psi_D_component);}
#endif

#ifndef COMPILE_WITHOUT_DYADIC_SUPPORT
    filename=basedir+"/%d/octree_psi_D";
    if(has_valid_octree_grid && FILE_UTILITIES::Frame_File_Exists(filename,start_frame)){
        OPENGL_COMPONENT_OCTREE_CELL_SCALAR_FIELD<T,bool>* octree_psi_D_component=new OPENGL_COMPONENT_OCTREE_CELL_SCALAR_FIELD<T,bool>(octree_grid_component->opengl_grid.grid,filename,psi_D_color_map);
        Add_Component(octree_psi_D_component,"Octree Psi_D points",'\0',BASIC_VISUALIZATION::OWNED|BASIC_VISUALIZATION::START_HIDDEN);
        opengl_world.Append_Bind_Key(OPENGL_KEY(OPENGL_KEY::F1),octree_psi_D_component->Toggle_Draw_CB());
        slice_manager.Add_Object(octree_psi_D_component);}
#endif

    OPENGL_INDEXED_COLOR_MAP* colors_color_map=OPENGL_INDEXED_COLOR_MAP::Basic_16_Color_Map();colors_color_map->Set_Index_Mode(OPENGL_INDEXED_COLOR_MAP::PERIODIC);
    filename=basedir+"/%d/colors";
    if(has_valid_grid && FILE_UTILITIES::Frame_File_Exists(filename,start_frame)){
        OPENGL_COMPONENT_SCALAR_FIELD_3D<T,int>* psi_colors_component=new OPENGL_COMPONENT_SCALAR_FIELD_3D<T,int>(mac_grid,filename,colors_color_map,OPENGL_SCALAR_FIELD_3D<T,int>::DRAW_POINTS);
        Add_Component(psi_colors_component,"Psi colors",'\0',BASIC_VISUALIZATION::OWNED|BASIC_VISUALIZATION::START_HIDDEN);
        opengl_world.Append_Bind_Key(OPENGL_KEY(OPENGL_KEY::F2),psi_colors_component->Toggle_Draw_CB());
        slice_manager.Add_Object(psi_colors_component);}

#ifndef COMPILE_WITHOUT_RLE_SUPPORT
    filename=basedir+"/%d/rle_colors";
    if(has_valid_rle_grid && FILE_UTILITIES::Frame_File_Exists(filename,start_frame)){
        OPENGL_COMPONENT_RLE_CELL_SCALAR_FIELD_3D<T,int>* rle_psi_colors_component=new OPENGL_COMPONENT_RLE_CELL_SCALAR_FIELD_3D<T,int>(rle_grid_component->opengl_grid.grid,filename,colors_color_map);
        Add_Component(rle_psi_colors_component,"RLE Psi colors",'\0',BASIC_VISUALIZATION::OWNED|BASIC_VISUALIZATION::START_HIDDEN);
        opengl_world.Append_Bind_Key(OPENGL_KEY(OPENGL_KEY::F2),rle_psi_colors_component->Toggle_Draw_CB());
        slice_manager.Add_Object(rle_psi_colors_component);}
#endif

#ifndef COMPILE_WITHOUT_DYADIC_SUPPORT
    filename=basedir+"/%d/octree_colors";
    if(has_valid_octree_grid && FILE_UTILITIES::Frame_File_Exists(filename,start_frame)){
        OPENGL_COMPONENT_OCTREE_CELL_SCALAR_FIELD<T,int>* octree_psi_colors_component=new OPENGL_COMPONENT_OCTREE_CELL_SCALAR_FIELD<T,int>(octree_grid_component->opengl_grid.grid,filename,colors_color_map);
        Add_Component(octree_psi_colors_component,"Octree Psi colors",'\0',BASIC_VISUALIZATION::OWNED|BASIC_VISUALIZATION::START_HIDDEN);
        opengl_world.Append_Bind_Key(OPENGL_KEY(OPENGL_KEY::F2),octree_psi_colors_component->Toggle_Draw_CB());
        slice_manager.Add_Object(octree_psi_colors_component);}
#endif

{filename=basedir+"/%d/strain";
    OPENGL_COMPONENT_SYMMETRIC_MATRIX_FIELD_3D<T>* strain_component=0; // TODO: make this not a hack for multiphase
    if(has_valid_grid && FILE_UTILITIES::Frame_File_Exists(filename,start_frame)){
        strain_component=new OPENGL_COMPONENT_SYMMETRIC_MATRIX_FIELD_3D<T>(grid,filename);
        Add_Component(strain_component,"Strain",'e',BASIC_VISUALIZATION::OWNED|BASIC_VISUALIZATION::START_HIDDEN);
        opengl_world.Append_Bind_Key('+',strain_component->Increase_Size_CB());
        opengl_world.Append_Bind_Key('_',strain_component->Decrease_Size_CB());
        slice_manager.Add_Object(strain_component);}

    filename=basedir+"/%d/strain_1";
    if(!strain_component && has_valid_grid && FILE_UTILITIES::Frame_File_Exists(filename,start_frame)){
        strain_component=new OPENGL_COMPONENT_SYMMETRIC_MATRIX_FIELD_3D<T>(grid,filename);
        Add_Component(strain_component,"Strain",'e',BASIC_VISUALIZATION::OWNED|BASIC_VISUALIZATION::START_HIDDEN);
        opengl_world.Append_Bind_Key('+',strain_component->Increase_Size_CB());
        opengl_world.Append_Bind_Key('_',strain_component->Decrease_Size_CB());
        slice_manager.Add_Object(strain_component);}}

#ifndef COMPILE_WITHOUT_DYADIC_SUPPORT
    filename=basedir+"/%d/octree_velocities";
    if(has_valid_octree_grid && FILE_UTILITIES::Frame_File_Exists(filename,start_frame)){
        OPENGL_COMPONENT_OCTREE_NODE_VECTOR_FIELD<T>* octree_velocity_component=new OPENGL_COMPONENT_OCTREE_NODE_VECTOR_FIELD<T>(octree_grid_component->opengl_grid.grid,filename);
        octree_velocity_component->opengl_octree_node_vector_field.vector_color=OPENGL_COLOR::Green();
        octree_velocity_component->opengl_octree_node_vector_field.size=0.01;
        Add_Component(octree_velocity_component,"Octree Velocities",'v',BASIC_VISUALIZATION::OWNED|BASIC_VISUALIZATION::START_HIDDEN);
        opengl_world.Append_Bind_Key('=',octree_velocity_component->Increase_Vector_Size_CB());
        opengl_world.Append_Bind_Key('-',octree_velocity_component->Decrease_Vector_Size_CB());
        opengl_world.Append_Bind_Key('h',octree_velocity_component->Toggle_Arrowhead_CB());
        slice_manager.Add_Object(octree_velocity_component);}

    filename=basedir+"/%d/octree_velocities_ghost";
    if(has_valid_octree_grid && FILE_UTILITIES::Frame_File_Exists(filename,start_frame)){
        OPENGL_COMPONENT_OCTREE_NODE_VECTOR_FIELD<T>* octree_velocity_minus_component=new OPENGL_COMPONENT_OCTREE_NODE_VECTOR_FIELD<T>(octree_grid_component->opengl_grid.grid,filename);
        octree_velocity_minus_component->opengl_octree_node_vector_field.vector_color=OPENGL_COLOR::Green();
        Add_Component(octree_velocity_minus_component,"Octree Minus Velocities",'c',BASIC_VISUALIZATION::OWNED|BASIC_VISUALIZATION::START_HIDDEN);
        opengl_world.Append_Bind_Key('=',octree_velocity_minus_component->Increase_Vector_Size_CB());
        opengl_world.Append_Bind_Key('-',octree_velocity_minus_component->Decrease_Vector_Size_CB());
        opengl_world.Append_Bind_Key('h',octree_velocity_minus_component->Toggle_Arrowhead_CB());
        slice_manager.Add_Object(octree_velocity_minus_component);}

    filename=basedir+"/%d/octree_velocities_ghost_fuel";
    if(has_valid_octree_grid && FILE_UTILITIES::Frame_File_Exists(filename,start_frame)){
        OPENGL_COMPONENT_OCTREE_NODE_VECTOR_FIELD<T>* octree_velocity_plus_component=new OPENGL_COMPONENT_OCTREE_NODE_VECTOR_FIELD<T>(octree_grid_component->opengl_grid.grid,filename);
        octree_velocity_plus_component->opengl_octree_node_vector_field.vector_color=OPENGL_COLOR::Green();
        Add_Component(octree_velocity_plus_component,"Octree Plus Velocities",'b',BASIC_VISUALIZATION::OWNED|BASIC_VISUALIZATION::START_HIDDEN);
        opengl_world.Append_Bind_Key('=',octree_velocity_plus_component->Increase_Vector_Size_CB());
        opengl_world.Append_Bind_Key('-',octree_velocity_plus_component->Decrease_Vector_Size_CB());
        opengl_world.Append_Bind_Key('h',octree_velocity_plus_component->Toggle_Arrowhead_CB());
        slice_manager.Add_Object(octree_velocity_plus_component);}
#endif

    opengl_world.Set_Key_Binding_Category("Motion");

    filename=basedir+"/%d/surface.tri";
    if(FILE_UTILITIES::Frame_File_Exists(filename,start_frame)){
        OPENGL_COMPONENT_TRIANGULATED_SURFACE<T,T>* triangulated_surface_component=new OPENGL_COMPONENT_TRIANGULATED_SURFACE<T,T>(filename);
        Add_Component(triangulated_surface_component,"Triangulated Surface",',',BASIC_VISUALIZATION::OWNED);}

    if(slice_manager.slice){
        opengl_world.Set_Key_Binding_Category("Slice Control");
        opengl_world.Append_Bind_Key("^h",slice_manager.Toggle_Slice_Mode_CB("Toggle 3D/Slice mode"));
        opengl_world.Append_Bind_Key('\\',slice_manager.Toggle_Slice_Axis_CB("Toggle slice axis"));
        opengl_world.Append_Bind_Key(']',slice_manager.Increment_Slice_CB("Increment slice"));
        opengl_world.Append_Bind_Key('[',slice_manager.Decrement_Slice_CB("Decrement slice"));}

    opengl_world.Set_Key_Binding_Category("Misc.");
    opengl_world.Append_Bind_Key('~',Command_Prompt_CB("Command prompt"));
    opengl_world.Append_Bind_Key(OPENGL_KEY(OPENGL_KEY::F5),Draw_All_Objects_CB());

    // initialize selection priority (highest on top)
    Selection_Priority(OPENGL_SELECTION::POINTS_3D)=100;
    Selection_Priority(OPENGL_SELECTION::COMPONENT_PARTICLES_3D)=100;
    Selection_Priority(OPENGL_SELECTION::ARTICULATED_RIGID_BODIES_JOINT_3D)=95;
    Selection_Priority(OPENGL_SELECTION::MUSCLE_3D)=95;
    Selection_Priority(OPENGL_SELECTION::MUSCLE_SURFACE_3D)=95;
    Selection_Priority(OPENGL_SELECTION::TRIANGULATED_SURFACE_VERTEX)=90;
    Selection_Priority(OPENGL_SELECTION::TRIANGULATED_SURFACE_SEGMENT)=89;
    Selection_Priority(OPENGL_SELECTION::TRIANGULATED_SURFACE_TRIANGLE)=88;
    Selection_Priority(OPENGL_SELECTION::TETRAHEDRALIZED_VOLUME_VERTEX)=85;
    Selection_Priority(OPENGL_SELECTION::TETRAHEDRALIZED_VOLUME_TETRAHEDRON)=84;
    Selection_Priority(OPENGL_SELECTION::COMPONENT_RIGID_BODIES_3D)=80;
    Selection_Priority(OPENGL_SELECTION::COMPONENT_DEFORMABLE_COLLECTION_3D)=80;
    Selection_Priority(OPENGL_SELECTION::SEGMENTED_CURVE_VERTEX_3D)=79;
    Selection_Priority(OPENGL_SELECTION::SEGMENTED_CURVE_SEGMENT_3D)=78;
    Selection_Priority(OPENGL_SELECTION::GRID_NODE_3D)=70;
#ifndef COMPILE_WITHOUT_DYADIC_SUPPORT
    Selection_Priority(OPENGL_SELECTION::OCTREE_NODE)=70;
    Selection_Priority(OPENGL_SELECTION::OCTREE_FACE)=65;
#endif
    Selection_Priority(OPENGL_SELECTION::GRID_CELL_3D)=60;
#ifndef COMPILE_WITHOUT_RLE_SUPPORT
    Selection_Priority(OPENGL_SELECTION::RLE_CELL_3D)=60;
#endif
#ifndef COMPILE_WITHOUT_DYADIC_SUPPORT
    Selection_Priority(OPENGL_SELECTION::OCTREE_CELL)=60;
#endif
}

template<class T,class RW> void VISUALIZATION<T,RW>::
Update_OpenGL_Strings()
{
    ANIMATED_VISUALIZATION::Update_OpenGL_Strings();

    // TODO: slice manager should be a component
    std::ostringstream output_stream;
    if(slice_manager.slice) slice_manager.slice->Print_Slice_Info(output_stream);
    opengl_world.Add_String(output_stream.str());
}

template<class T,class RW> void VISUALIZATION<T,RW>::
Pre_Frame_Extra()
{
#ifndef COMPILE_WITHOUT_RLE_SUPPORT
    if(rle_grid_component) rle_grid_component->Set_Frame(frame);
#endif
#ifndef COMPILE_WITHOUT_DYADIC_SUPPORT
    if(octree_grid_component) octree_grid_component->Set_Frame(frame);
#endif
    if(grid_component) grid_component->Set_Frame(frame);
}

template<class T,class RW> void VISUALIZATION<T,RW>::
Set_Frame_Extra()
{
    std::string filename=STRING_UTILITIES::string_sprintf("%s/%d/frame_title",basedir.c_str(),frame);
    if(FILE_UTILITIES::File_Exists(filename)){std::ifstream input(filename.c_str());getline(input,frame_title);}
    else frame_title="";
    filename=STRING_UTILITIES::string_sprintf("%s/%d/time",basedir.c_str(),frame);
    if(FILE_UTILITIES::File_Exists(filename)){T time;FILE_UTILITIES::template Read_From_File<RW>(filename,time);frame_title=STRING_UTILITIES::string_sprintf("(%.05f) ",time)+frame_title;}
}

template<class T,class RW> void VISUALIZATION<T,RW>::
Command_Prompt_Response()
{
    if(!opengl_world.prompt_response.empty()){
        std::string command;
        std::istringstream sstream(opengl_world.prompt_response);
        sstream>>command;
        if(command=="s"){
            int id;
            if(sstream>>id){
                OPENGL_SELECTION* selection=0;
                if(selection){Set_Current_Selection(selection); Update_OpenGL_Strings();}}}
        else if(command=="n"){
            int id;
            if(sstream>>id){
                OPENGL_SELECTION* selection=0;
#ifndef COMPILE_WITHOUT_DYADIC_SUPPORT
                if(octree_grid_component){selection=octree_grid_component->opengl_grid.Get_Node_Selection(id);}
#endif
                if(selection){Set_Current_Selection(selection); Update_OpenGL_Strings();}}}
        else if(command=="c"){
            int id;
            if(sstream>>id){
                OPENGL_SELECTION* selection=0;
#ifndef COMPILE_WITHOUT_DYADIC_SUPPORT
                if(octree_grid_component){selection=octree_grid_component->opengl_grid.Get_Cell_Selection(id);}
#endif
                if(selection){Set_Current_Selection(selection); Update_OpenGL_Strings();}}}
    }
}

template<class T,class RW> void VISUALIZATION<T,RW>::
Command_Prompt()
{
    opengl_world.Prompt_User("Command: ",Command_Prompt_Response_CB());
}

template<class T,class RW> void VISUALIZATION<T,RW>::
Slice_Has_Changed()
{
    Update_OpenGL_Strings();
}

// ==========================================================================

int main(int argc,char* argv[])
{
    Initialize_Geometry_Particle();
    Initialize_Read_Write_Structures();
    PROCESS_UTILITIES::Set_Floating_Point_Exception_Handling(true);
    PROCESS_UTILITIES::Set_Backtrace(true);
    bool type_double=false; // float by default
    if(PARSE_ARGS::Find_And_Remove("-float",argc,argv)) type_double=false;
    if(PARSE_ARGS::Find_And_Remove("-double",argc,argv)) type_double=true;

    ANIMATED_VISUALIZATION* visualization=0;
    if(!type_double) visualization=new VISUALIZATION<float>();
#ifndef COMPILE_WITHOUT_DOUBLE_SUPPORT
    else visualization=new VISUALIZATION<double>();
#else
    else{std::cerr<<"Double support not enabled."<<std::endl;exit(1);}
#endif
    visualization->Initialize_And_Run(argc,argv);
    return 0;
}
