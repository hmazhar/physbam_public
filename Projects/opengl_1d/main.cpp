//#####################################################################
// Copyright 2006-2009, Jon Gretarsson, Nipun Kwatra, Andrew Selle, Jonathan Su.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#include <PhysBAM_Tools/Grids_Uniform/GRID.h>
#include <PhysBAM_Tools/Log/LOG.h>
#include <PhysBAM_Tools/Parsing/PARSE_ARGS.h>
#include <PhysBAM_Tools/Read_Write/Grids_Uniform/READ_WRITE_GRID.h>
#include <PhysBAM_Tools/Read_Write/Utilities/FILE_UTILITIES.h>
#include <PhysBAM_Tools/Utilities/PROCESS_UTILITIES.h>
#include <PhysBAM_Geometry/Geometry_Particles/REGISTER_GEOMETRY_READ_WRITE.h>
#include <PhysBAM_Geometry/Grids_Uniform_Level_Sets/LEVELSET_1D.h>
#include <PhysBAM_Geometry/Read_Write/Grids_Uniform_Level_Sets/READ_WRITE_LEVELSET_2D.h>
#include <PhysBAM_Rendering/PhysBAM_OpenGL/OpenGL/ANIMATED_VISUALIZATION.h>
#include <PhysBAM_Rendering/PhysBAM_OpenGL/OpenGL/OPENGL_CONSTANT_COLOR_MAP.h>
#include <PhysBAM_Rendering/PhysBAM_OpenGL/OpenGL/OPENGL_GRID_1D.h>
#include <PhysBAM_Rendering/PhysBAM_OpenGL/OpenGL_Components/OPENGL_COMPONENT_BASIC.h>
#include <PhysBAM_Rendering/PhysBAM_OpenGL/OpenGL_Components/OPENGL_COMPONENT_BINTREE_CELL_SCALAR_FIELD.h>
#include <PhysBAM_Rendering/PhysBAM_OpenGL/OpenGL_Components/OPENGL_COMPONENT_BINTREE_FACE_SCALAR_FIELD.h>
#include <PhysBAM_Rendering/PhysBAM_OpenGL/OpenGL_Components/OPENGL_COMPONENT_BINTREE_GRID.h>
#include <PhysBAM_Rendering/PhysBAM_OpenGL/OpenGL_Components/OPENGL_COMPONENT_DEFORMABLE_GEOMETRY_COLLECTION_1D.h>
#include <PhysBAM_Rendering/PhysBAM_OpenGL/OpenGL_Components/OPENGL_COMPONENT_FACE_SCALAR_FIELD_1D.h>
#include <PhysBAM_Rendering/PhysBAM_OpenGL/OpenGL_Components/OPENGL_COMPONENT_LEVELSET_1D.h>
#include <PhysBAM_Rendering/PhysBAM_OpenGL/OpenGL_Components/OPENGL_COMPONENT_RIGID_GEOMETRY_COLLECTION_1D.h>
#include <PhysBAM_Rendering/PhysBAM_OpenGL/OpenGL_Components/OPENGL_COMPONENT_SCALAR_FIELD_1D.h>
#include <fstream>
#include <sstream>
#include <string.h>

using namespace PhysBAM;

const char* DEFAULT_BASEDIR=".";

template<class T,class RW=T>
class OPENGL_1D_VISUALIZATION:public ANIMATED_VISUALIZATION
{
    typedef VECTOR<T,1> TV;
public:
    OPENGL_1D_VISUALIZATION();
    ~OPENGL_1D_VISUALIZATION();

private:
    void Read_Grid();
protected:
    virtual void Add_Arguments(PARSE_ARGS& parse_args);
    virtual void Parse_Arguments(PARSE_ARGS& parse_args);
    virtual void Initialize_Components_And_Key_Bindings();
    virtual void Add_OpenGL_Initialization();
    virtual void Pre_Frame_Extra();
    virtual void Set_Frame_Extra();

    // Options
    std::string basedir;

    bool has_valid_grid;

    bool node_based;

    GRID<TV> grid,mac_grid,regular_grid;
    OPENGL_COMPONENT_BASIC<OPENGL_GRID_1D<T> >* grid_component;
#if COMPILE_WITH_BINTREE_SUPPORT
    bool has_valid_bintree_grid;
    OPENGL_COMPONENT_BINTREE_GRID<T>* bintree_grid_component;
#endif
};

//#####################################################################
// Function OPENGL_1D_VISUALIZATION
//#####################################################################
template<class T,class RW> OPENGL_1D_VISUALIZATION<T,RW>::
OPENGL_1D_VISUALIZATION()
    :grid_component(0)
#if COMPILE_WITH_BINTREE_SUPPORT
     ,bintree_grid_component(0)
#endif
{
    add_axes=false;
}
//#####################################################################
// Function OPENGL_1D_VISUALIZATION
//#####################################################################
template<class T,class RW> OPENGL_1D_VISUALIZATION<T,RW>::
~OPENGL_1D_VISUALIZATION()
{
    delete &grid_component->object;
    delete grid_component;
#if COMPILE_WITH_BINTREE_SUPPORT
     delete bintree_grid_component;
#endif
}
//#####################################################################
// Function Add_Arguments
//#####################################################################
template<class T,class RW> void OPENGL_1D_VISUALIZATION<T,RW>::
Add_Arguments(PARSE_ARGS& parse_args)
{
    basedir=DEFAULT_BASEDIR;   // default basedir

    ANIMATED_VISUALIZATION::Add_Arguments(parse_args);

    parse_args.Add_Option_Argument("-float");
    parse_args.Add_Option_Argument("-double");
    parse_args.Add_Option_Argument("-no_ghost");
    parse_args.Set_Extra_Arguments(-1,"[<basedir>]");
}
//#####################################################################
// Function Parse_Arguments
//#####################################################################
template<class T,class RW> void OPENGL_1D_VISUALIZATION<T,RW>::
Parse_Arguments(PARSE_ARGS& parse_args)
{
    ANIMATED_VISUALIZATION::Parse_Arguments(parse_args);

    if(parse_args.Num_Extra_Args()>1)
        parse_args.Print_Usage(true);
    else if(parse_args.Num_Extra_Args()==1)
        basedir=parse_args.Extra_Arg(1);

    last_frame_filename=basedir+"/common/last_frame";

    if(FILE_UTILITIES::File_Exists(basedir+"/common/first_frame") && !parse_args.Is_Value_Set("-start_frame")) FILE_UTILITIES::Read_From_Text_File(basedir+"/common/first_frame",start_frame);

    // don't override camera script filename if it was already set in base class based on command line argument
    if(camera_script_filename.empty()) camera_script_filename=basedir+"/camera_script";

#ifdef __linux__
    if(!parse_args.Is_Value_Set("-window_title")) opengl_window_title="opengl_1d: " + FILE_UTILITIES::Real_Path(basedir);
#endif
}
//#####################################################################
// Function Read_Grid
//#####################################################################
template<class T,class RW> void OPENGL_1D_VISUALIZATION<T,RW>::
Read_Grid()
{
    has_valid_grid=false;
    std::string filename;
    if(FILE_UTILITIES::File_Exists(basedir+"/common/grid")){
        filename=basedir+"/common/grid";
        LOG::cout<<"Reading grid from '"<<filename<<"'..."<<std::flush;
        FILE_UTILITIES::Read_From_File<RW>(filename,grid);
        has_valid_grid=true;}
    if(has_valid_grid){
        if(!grid.MAC_offset){
            regular_grid=grid;
            mac_grid.Initialize(grid.counts.x-1,grid.Domain(),true);
            node_based=true;}
        else{
            mac_grid=grid;
            regular_grid.Initialize(grid.counts.x+1,grid.Domain(),false);
            node_based=false;}
        LOG::cout<<"regular grid "<<regular_grid<<" mac grid "<<mac_grid<<" node_based "<<node_based<<std::endl;}

#if COMPILE_WITH_BINTREE_SUPPORT
        std::string bintree_filename=STRING_UTILITIES::string_sprintf("%s/%d/dyadic_grid",basedir.c_str(),start_frame);
        if(FILE_UTILITIES::File_Exists(bintree_filename)){
            LOG::cout<<"Reading bintree '"<<bintree_filename<<"'"<<std::endl<<std::flush;
            bintree_filename=STRING_UTILITIES::string_sprintf("%s/%%d/dyadic_grid",basedir.c_str());
            bintree_grid_component=new OPENGL_COMPONENT_BINTREE_GRID<T>(bintree_filename);
            LOG::cout<<"Reading bintree '"<<bintree_filename<<"'"<<std::endl<<std::flush;
            // Add_Component(bintree_grid_component); // add later instead of here, but Set_Frame_Extra will ensure bintree grid is read in
            has_valid_bintree_grid=true;}
#endif
}
//#####################################################################
// Function Initialize_Components
//#####################################################################
template<class T,class RW> void OPENGL_1D_VISUALIZATION<T,RW>::
Initialize_Components_And_Key_Bindings()
{
    ANIMATED_VISUALIZATION::Initialize_Components_And_Key_Bindings();
    std::string filename;

    Read_Grid();

    opengl_world.Unbind_Keys("delmo7vDEMO&V+-");


    opengl_world.Set_Key_Binding_Category_Priority(1);
    opengl_world.Set_Key_Binding_Category("Compressible");

    if(has_valid_grid){
        OPENGL_GRID_1D<T>* opengl_grid=new OPENGL_GRID_1D<T>(grid,OPENGL_COLOR::Gray(.5),basedir,frame);
        grid_component=new OPENGL_COMPONENT_BASIC<OPENGL_GRID_1D<T> >(*opengl_grid);
        Add_Component(grid_component,"Grid",'6',BASIC_VISUALIZATION::SELECTABLE);}

    filename=basedir+"/%d/u";
    if(FILE_UTILITIES::Frame_File_Exists(filename,start_frame))
        Add_Component(new OPENGL_COMPONENT_SCALAR_FIELD_1D<T,T,RW>(grid,filename,OPENGL_COLOR::Red(),OPENGL_COLOR::Magenta()),
            "u",'\0',BASIC_VISUALIZATION::OWNED | BASIC_VISUALIZATION::SELECTABLE);
    filename=basedir+"/%d/rigid_geometry_particles";
    if(FILE_UTILITIES::Frame_File_Exists(filename,start_frame)){
        OPENGL_COMPONENT_RIGID_GEOMETRY_COLLECTION_1D<T>* rigid_body_component=new OPENGL_COMPONENT_RIGID_GEOMETRY_COLLECTION_1D<T>(basedir,true);
        Add_Component(rigid_body_component,"Rigid Bodies",'5',BASIC_VISUALIZATION::OWNED);
        opengl_world.Append_Bind_Key('n',rigid_body_component->Toggle_Show_Object_Names_CB());}
    filename=basedir+"/%d/deformable_object_particles";
    if(FILE_UTILITIES::Frame_File_Exists(filename,start_frame)){
        OPENGL_COMPONENT_DEFORMABLE_GEOMETRY_COLLECTION_1D<T>* deformable_geometry_component=new OPENGL_COMPONENT_DEFORMABLE_GEOMETRY_COLLECTION_1D<T>(basedir,start_frame);
        Add_Component(deformable_geometry_component,"Deformable Bodies",'8',BASIC_VISUALIZATION::OWNED);}
    filename=basedir+"/%d/u_exact";
    if(FILE_UTILITIES::Frame_File_Exists(filename,start_frame))
        Add_Component(new OPENGL_COMPONENT_SCALAR_FIELD_1D<T,T,RW>(grid,filename,OPENGL_COLOR::Blue(),OPENGL_COLOR::Cyan()),
            "u_exact",'\0',BASIC_VISUALIZATION::OWNED | BASIC_VISUALIZATION::SELECTABLE);
    filename=basedir+"/%d/levelset";
    if(FILE_UTILITIES::Frame_File_Exists(filename,start_frame))
        Add_Component(new OPENGL_COMPONENT_LEVELSET_1D<T>(grid,filename,OPENGL_COLOR::Yellow(),OPENGL_COLOR::Yellow()),
            "levelset",'l',BASIC_VISUALIZATION::OWNED | BASIC_VISUALIZATION::SELECTABLE);
    // Compressible
    filename=basedir+"/%d/density";
    if(FILE_UTILITIES::Frame_File_Exists(filename,start_frame))
        Add_Component(new OPENGL_COMPONENT_SCALAR_FIELD_1D<T,T,RW>(grid,filename,OPENGL_COLOR::Yellow(),OPENGL_COLOR::Yellow()),
            "density",'1',BASIC_VISUALIZATION::OWNED | BASIC_VISUALIZATION::SELECTABLE);
    filename=basedir+"/%d/centered_velocities";
    if(FILE_UTILITIES::Frame_File_Exists(filename,start_frame))
        Add_Component(new OPENGL_COMPONENT_SCALAR_FIELD_1D<T,T,RW>(grid,filename,OPENGL_COLOR::Magenta(),OPENGL_COLOR::Magenta()),
            "velocity",'v',BASIC_VISUALIZATION::OWNED | BASIC_VISUALIZATION::SELECTABLE);
    filename=basedir+"/%d/mac_velocities";
    if(FILE_UTILITIES::Frame_File_Exists(filename,start_frame))
        Add_Component(new OPENGL_COMPONENT_FACE_SCALAR_FIELD_1D<T,T,RW>(grid,filename,OPENGL_COLOR::Red(),OPENGL_COLOR::Red()),
            "mac_velocities",'a',BASIC_VISUALIZATION::OWNED | BASIC_VISUALIZATION::SELECTABLE | BASIC_VISUALIZATION::START_HIDDEN);
    filename=basedir+"/%d/momentum";
    if(FILE_UTILITIES::Frame_File_Exists(filename,start_frame))
        Add_Component(new OPENGL_COMPONENT_SCALAR_FIELD_1D<T,T,RW>(grid,filename,OPENGL_COLOR::Magenta(),OPENGL_COLOR::Magenta()),
            "momentum",'2',BASIC_VISUALIZATION::OWNED | BASIC_VISUALIZATION::SELECTABLE | BASIC_VISUALIZATION::START_HIDDEN);
    filename=basedir+"/%d/energy";
    if(FILE_UTILITIES::Frame_File_Exists(filename,start_frame))
        Add_Component(new OPENGL_COMPONENT_SCALAR_FIELD_1D<T,T,RW>(grid,filename,OPENGL_COLOR::Cyan(),OPENGL_COLOR::Cyan()),
            "energy",'3',BASIC_VISUALIZATION::OWNED | BASIC_VISUALIZATION::SELECTABLE | BASIC_VISUALIZATION::START_HIDDEN);
    filename=basedir+"/%d/internal_energy";
    if(FILE_UTILITIES::Frame_File_Exists(filename,start_frame))
        Add_Component(new OPENGL_COMPONENT_SCALAR_FIELD_1D<T,T,RW>(grid,filename,OPENGL_COLOR::Cyan(),OPENGL_COLOR::Cyan()),
            "internal energy",'4',BASIC_VISUALIZATION::OWNED | BASIC_VISUALIZATION::SELECTABLE | BASIC_VISUALIZATION::START_HIDDEN);
    filename=basedir+"/%d/pressure";
    if(FILE_UTILITIES::Frame_File_Exists(filename,start_frame))
        Add_Component(new OPENGL_COMPONENT_SCALAR_FIELD_1D<T,T,RW>(grid,filename,OPENGL_COLOR::Cyan(),OPENGL_COLOR::Cyan()),
            "pressure",'7',BASIC_VISUALIZATION::OWNED | BASIC_VISUALIZATION::SELECTABLE);
    filename=basedir+"/%d/compressible_implicit_pressure";
    if(FILE_UTILITIES::Frame_File_Exists(filename,start_frame)){
        OPENGL_COMPONENT_SCALAR_FIELD_1D<T,T,RW>* compressible_implicit_pressure_component=
            new OPENGL_COMPONENT_SCALAR_FIELD_1D<T,T,RW>(grid,filename,OPENGL_COLOR::Red(),OPENGL_COLOR::Red());
        Add_Component(compressible_implicit_pressure_component,
            "compressible_implicit_pressure",'9',BASIC_VISUALIZATION::OWNED | BASIC_VISUALIZATION::SELECTABLE | BASIC_VISUALIZATION::START_HIDDEN);}
    filename=basedir+"/%d/entropy";
    if(FILE_UTILITIES::Frame_File_Exists(filename,start_frame))
        Add_Component(new OPENGL_COMPONENT_SCALAR_FIELD_1D<T,T,RW>(grid,filename,OPENGL_COLOR::Green(),OPENGL_COLOR::Green()),
            "entropy",'e',BASIC_VISUALIZATION::OWNED | BASIC_VISUALIZATION::SELECTABLE | BASIC_VISUALIZATION::START_HIDDEN);
    filename=basedir+"/%d/speedofsound";
    if(FILE_UTILITIES::Frame_File_Exists(filename,start_frame))
        Add_Component(new OPENGL_COMPONENT_SCALAR_FIELD_1D<T,T,RW>(grid,filename,OPENGL_COLOR::Red(),OPENGL_COLOR::Red()),
            "speedofsound",'o',BASIC_VISUALIZATION::OWNED | BASIC_VISUALIZATION::SELECTABLE | BASIC_VISUALIZATION::START_HIDDEN);
    filename=basedir+"/%d/machnumber";
    if(FILE_UTILITIES::Frame_File_Exists(filename,start_frame))
        Add_Component(new OPENGL_COMPONENT_SCALAR_FIELD_1D<T,T,RW>(grid,filename,OPENGL_COLOR::Blue(),OPENGL_COLOR::Magenta()),
            "machnumber",'m',BASIC_VISUALIZATION::OWNED | BASIC_VISUALIZATION::SELECTABLE | BASIC_VISUALIZATION::START_HIDDEN);
    filename=basedir+"/%d/velocity_plus_c";
    if(FILE_UTILITIES::Frame_File_Exists(filename,start_frame))
        Add_Component(new OPENGL_COMPONENT_SCALAR_FIELD_1D<T,T,RW>(grid,filename,OPENGL_COLOR::Green(),OPENGL_COLOR::Green()),
            "velocityplusc",'+',BASIC_VISUALIZATION::OWNED | BASIC_VISUALIZATION::SELECTABLE | BASIC_VISUALIZATION::START_HIDDEN);
    filename=basedir+"/%d/velocity_minus_c";
    if(FILE_UTILITIES::Frame_File_Exists(filename,start_frame))
        Add_Component(new OPENGL_COMPONENT_SCALAR_FIELD_1D<T,T,RW>(grid,filename,OPENGL_COLOR::Blue(),OPENGL_COLOR::Blue()),
            "velocityminusc",'-',BASIC_VISUALIZATION::OWNED | BASIC_VISUALIZATION::SELECTABLE | BASIC_VISUALIZATION::START_HIDDEN);
    filename=basedir+"_exact/%d/density";
    if(FILE_UTILITIES::Frame_File_Exists(filename,start_frame))
        Add_Component(new OPENGL_COMPONENT_SCALAR_FIELD_1D<T,T,RW>(grid,filename,OPENGL_COLOR::Blue(),OPENGL_COLOR::Cyan()),
            "density_exact",'!',BASIC_VISUALIZATION::OWNED | BASIC_VISUALIZATION::SELECTABLE | BASIC_VISUALIZATION::START_HIDDEN);
    filename=basedir+"_exact/%d/centered_velocities";
    if(FILE_UTILITIES::Frame_File_Exists(filename,start_frame))
        Add_Component(new OPENGL_COMPONENT_SCALAR_FIELD_1D<T,T,RW>(grid,filename,OPENGL_COLOR::Green(),OPENGL_COLOR::Yellow()),
            "velocity_exact",'V',BASIC_VISUALIZATION::OWNED | BASIC_VISUALIZATION::SELECTABLE | BASIC_VISUALIZATION::START_HIDDEN);
    filename=basedir+"_exact/%d/energy";
    if(FILE_UTILITIES::Frame_File_Exists(filename,start_frame))
        Add_Component(new OPENGL_COMPONENT_SCALAR_FIELD_1D<T,T,RW>(grid,filename,OPENGL_COLOR::Red(),OPENGL_COLOR::Red()),
            "energy_exact",'#',BASIC_VISUALIZATION::OWNED | BASIC_VISUALIZATION::SELECTABLE | BASIC_VISUALIZATION::START_HIDDEN);
    filename=basedir+"_exact/%d/internal_energy";
    if(FILE_UTILITIES::Frame_File_Exists(filename,start_frame))
        Add_Component(new OPENGL_COMPONENT_SCALAR_FIELD_1D<T,T,RW>(grid,filename,OPENGL_COLOR::Red(),OPENGL_COLOR::Red()),
            "internal energy_exact",'$',BASIC_VISUALIZATION::OWNED | BASIC_VISUALIZATION::SELECTABLE | BASIC_VISUALIZATION::START_HIDDEN);
    filename=basedir+"_exact/%d/pressure";
    if(FILE_UTILITIES::Frame_File_Exists(filename,start_frame))
        Add_Component(new OPENGL_COMPONENT_SCALAR_FIELD_1D<T,T,RW>(grid,filename,OPENGL_COLOR::Red(),OPENGL_COLOR::Magenta()),
            "pressure_exact",'&',BASIC_VISUALIZATION::OWNED | BASIC_VISUALIZATION::SELECTABLE | BASIC_VISUALIZATION::START_HIDDEN);
    filename=basedir+"_exact/%d/entropy";
    if(FILE_UTILITIES::Frame_File_Exists(filename,start_frame))
        Add_Component(new OPENGL_COMPONENT_SCALAR_FIELD_1D<T,T,RW>(grid,filename,OPENGL_COLOR::Blue(),OPENGL_COLOR::Cyan()),
            "entropy_exact",'E',BASIC_VISUALIZATION::OWNED | BASIC_VISUALIZATION::SELECTABLE | BASIC_VISUALIZATION::START_HIDDEN);
    filename=basedir+"_exact/%d/speedofsound";
    if(FILE_UTILITIES::Frame_File_Exists(filename,start_frame))
        Add_Component(new OPENGL_COMPONENT_SCALAR_FIELD_1D<T,T,RW>(grid,filename,OPENGL_COLOR::Green(),OPENGL_COLOR::Yellow()),
            "speedofsound_exact",'O',BASIC_VISUALIZATION::OWNED | BASIC_VISUALIZATION::SELECTABLE | BASIC_VISUALIZATION::START_HIDDEN);
    filename=basedir+"_exact/%d/machnumber";
    if(FILE_UTILITIES::Frame_File_Exists(filename,start_frame)){
        Add_Component(new OPENGL_COMPONENT_SCALAR_FIELD_1D<T,T,RW>(grid,filename,OPENGL_COLOR::Red(),OPENGL_COLOR::Magenta()),
            "machnumber_exact",'M',BASIC_VISUALIZATION::OWNED | BASIC_VISUALIZATION::SELECTABLE | BASIC_VISUALIZATION::START_HIDDEN);}
    filename=basedir+"/%d/p_cavitation";
    if(FILE_UTILITIES::Frame_File_Exists(filename,start_frame)){
        Add_Component(new OPENGL_COMPONENT_SCALAR_FIELD_1D<T,T,RW>(grid,filename,OPENGL_COLOR::Red(),OPENGL_COLOR::Magenta()),
            "p_cavitation",'[',BASIC_VISUALIZATION::OWNED | BASIC_VISUALIZATION::SELECTABLE | BASIC_VISUALIZATION::START_HIDDEN);}
    filename=basedir+"/%d/p_internal_energy";
    if(FILE_UTILITIES::Frame_File_Exists(filename,start_frame)){
        Add_Component(new OPENGL_COMPONENT_SCALAR_FIELD_1D<T,T,RW>(grid,filename,OPENGL_COLOR::Red(),OPENGL_COLOR::Magenta()),
            "p_internal_energy",']',BASIC_VISUALIZATION::OWNED | BASIC_VISUALIZATION::SELECTABLE | BASIC_VISUALIZATION::START_HIDDEN);}

    filename=basedir+"/%d/psi_N";
    if(FILE_UTILITIES::Frame_File_Exists(filename,start_frame)){
        OPENGL_COMPONENT_FACE_SCALAR_FIELD_1D<T,bool>* psi_N_component=new OPENGL_COMPONENT_FACE_SCALAR_FIELD_1D<T,bool>(grid,filename,OPENGL_COLOR::Cyan(),OPENGL_COLOR::Cyan());
        Add_Component(psi_N_component,"Psi_N points",'\0',BASIC_VISUALIZATION::START_HIDDEN|BASIC_VISUALIZATION::OWNED);
        opengl_world.Append_Bind_Key(OPENGL_KEY(OPENGL_KEY::F1),psi_N_component->Toggle_Draw_CB());}
    filename=basedir+"/%d/psi_D";
    if(FILE_UTILITIES::Frame_File_Exists(filename,start_frame))
        Add_Component(new OPENGL_COMPONENT_SCALAR_FIELD_1D<T,bool,RW>(grid,filename,OPENGL_COLOR::Red(),OPENGL_COLOR::Red()),
            "psi_D",'D',BASIC_VISUALIZATION::OWNED | BASIC_VISUALIZATION::SELECTABLE | BASIC_VISUALIZATION::START_HIDDEN);
    filename=basedir+"/%d/euler_psi";
    if(FILE_UTILITIES::Frame_File_Exists(filename,start_frame))
        Add_Component(new OPENGL_COMPONENT_SCALAR_FIELD_1D<T,bool,RW>(grid,filename,OPENGL_COLOR::Red(),OPENGL_COLOR::Red()),
            "euler_psi",'P',BASIC_VISUALIZATION::OWNED | BASIC_VISUALIZATION::SELECTABLE | BASIC_VISUALIZATION::START_HIDDEN);


    // add scaling controls to scalar fields
    opengl_world.Set_Key_Binding_Category("Scaling");
    for(int c=1;c<=component_list.m;c++){
        if(OPENGL_COMPONENT_SCALAR_FIELD_1D<T,T,RW>* scalar_field_component=dynamic_cast<OPENGL_COMPONENT_SCALAR_FIELD_1D<T,T,RW>*>(component_list(c))){
            opengl_world.Append_Bind_Key('>',scalar_field_component->Increase_Scale_CB());
            opengl_world.Append_Bind_Key('<',scalar_field_component->Decrease_Scale_CB());}
        else if(OPENGL_COMPONENT_FACE_SCALAR_FIELD_1D<T,T,RW>* face_scalar_field_component=dynamic_cast<OPENGL_COMPONENT_FACE_SCALAR_FIELD_1D<T,T,RW>*>(component_list(c))){
            opengl_world.Append_Bind_Key('>',face_scalar_field_component->Increase_Scale_CB());
            opengl_world.Append_Bind_Key('<',face_scalar_field_component->Decrease_Scale_CB());}}

#if COMPILE_WITH_BINTREE_SUPPORT
    if(bintree_grid_component){
        LOG::cout<<"Adding bintree grid component"<<std::endl;
        Add_Component(bintree_grid_component,"Bintree Grid",'6',BASIC_VISUALIZATION::SELECTABLE | BASIC_VISUALIZATION::OWNED);}

    filename=basedir+"/%d/dyadic_rho";
    if(has_valid_bintree_grid && FILE_UTILITIES::Frame_File_Exists(filename,start_frame)){
        OPENGL_COMPONENT_BINTREE_CELL_SCALAR_FIELD<T>* bintree_density_component=new OPENGL_COMPONENT_BINTREE_CELL_SCALAR_FIELD<T>(bintree_grid_component->opengl_grid.grid,filename,new OPENGL_CONSTANT_COLOR_MAP<T>(OPENGL_COLOR::Yellow()),false);
        Add_Component(bintree_density_component,"Bintree Density",'d',BASIC_VISUALIZATION::SELECTABLE | BASIC_VISUALIZATION::OWNED);}

    filename=basedir+"/%d/dyadic_face_velocities";
    if(has_valid_bintree_grid && FILE_UTILITIES::Frame_File_Exists(filename,start_frame)){
        OPENGL_COMPONENT_BINTREE_FACE_SCALAR_FIELD<T>* bintree_velocity_component=new OPENGL_COMPONENT_BINTREE_FACE_SCALAR_FIELD<T>(bintree_grid_component->opengl_grid.grid,filename,new OPENGL_CONSTANT_COLOR_MAP<T>(OPENGL_COLOR::Red()),true);
        Add_Component(bintree_velocity_component,"Bintree Face Velocities",'V',BASIC_VISUALIZATION::SELECTABLE | BASIC_VISUALIZATION::OWNED);}

    filename=basedir+"/%d/dyadic_velocities";
    if(has_valid_bintree_grid && FILE_UTILITIES::Frame_File_Exists(filename,start_frame)){
        OPENGL_COMPONENT_BINTREE_CELL_SCALAR_FIELD<T>* bintree_velocity_component=new OPENGL_COMPONENT_BINTREE_CELL_SCALAR_FIELD<T>(bintree_grid_component->opengl_grid.grid,filename,new OPENGL_CONSTANT_COLOR_MAP<T>(OPENGL_COLOR::Magenta()),false);
        Add_Component(bintree_velocity_component,"Bintree Velocities",'v',BASIC_VISUALIZATION::SELECTABLE | BASIC_VISUALIZATION::OWNED);}

    filename=basedir+"/%d/dyadic_psi_N";
    if(has_valid_bintree_grid && FILE_UTILITIES::Frame_File_Exists(filename,start_frame)){
        OPENGL_COMPONENT_BINTREE_FACE_SCALAR_FIELD<T,bool>* bintree_psi_N_component=new OPENGL_COMPONENT_BINTREE_FACE_SCALAR_FIELD<T,bool>(bintree_grid_component->opengl_grid.grid,filename,new OPENGL_CONSTANT_COLOR_MAP<bool>(OPENGL_COLOR::Cyan()),true);
        bintree_psi_N_component->Set_Draw(false);
        Add_Component(bintree_psi_N_component,"Bintree Psi_N points",'\0',BASIC_VISUALIZATION::SELECTABLE | BASIC_VISUALIZATION::OWNED);
        opengl_world.Append_Bind_Key(OPENGL_KEY(OPENGL_KEY::F1),bintree_psi_N_component->Toggle_Draw_CB());}

    filename=basedir+"/%d/dyadic_psi_D";
    if(has_valid_bintree_grid && FILE_UTILITIES::Frame_File_Exists(filename,start_frame)){
        OPENGL_COMPONENT_BINTREE_CELL_SCALAR_FIELD<T,bool>* bintree_psi_D_component=new OPENGL_COMPONENT_BINTREE_CELL_SCALAR_FIELD<T,bool>(bintree_grid_component->opengl_grid.grid,filename,new OPENGL_CONSTANT_COLOR_MAP<bool>(OPENGL_COLOR::Magenta()),true);
        bintree_psi_D_component->Set_Draw(false);
        Add_Component(bintree_psi_D_component,"Bintree Psi_D points",'\0',BASIC_VISUALIZATION::SELECTABLE | BASIC_VISUALIZATION::OWNED);
        opengl_world.Append_Bind_Key(OPENGL_KEY(OPENGL_KEY::F2),bintree_psi_D_component->Toggle_Draw_CB());}

#endif
}
//#####################################################################
// Function Set_Frame_Extra
//#####################################################################
template<class T,class RW> void OPENGL_1D_VISUALIZATION<T,RW>::
Set_Frame_Extra()
{
    std::string filename=STRING_UTILITIES::string_sprintf("%s/%d/frame_title",basedir.c_str(),frame);
    if(FILE_UTILITIES::File_Exists(filename)){std::ifstream input(filename.c_str());getline(input,frame_title);}
    else frame_title="";
    filename=STRING_UTILITIES::string_sprintf("%s/%d/time",basedir.c_str(),frame);
    if(FILE_UTILITIES::File_Exists(filename)){T time;FILE_UTILITIES::template Read_From_File<RW>(filename,time);frame_title=STRING_UTILITIES::string_sprintf("(%.05f) ",time)+frame_title;}
}
//#####################################################################
// Function Pre_Frame_Extra
//#####################################################################
template<class T,class RW> void OPENGL_1D_VISUALIZATION<T,RW>::
Pre_Frame_Extra()
{
#if COMPILE_WITH_BINTREE_SUPPORT
    if(bintree_grid_component) bintree_grid_component->Set_Frame(frame);
#endif
}
//#####################################################################
// Function Add_OpenGL_Initialization
//#####################################################################
template<class T,class RW> void OPENGL_1D_VISUALIZATION<T,RW>::
Add_OpenGL_Initialization()
{
    ANIMATED_VISUALIZATION::Add_OpenGL_Initialization();
    opengl_world.Set_2D_Mode(true);
}
//#####################################################################

int main(int argc,char *argv[])
{
//    Initialize_General_Particle();
    Initialize_Read_Write_Structures();
    Initialize_Geometry_Particle();
    PROCESS_UTILITIES::Set_Floating_Point_Exception_Handling(true);
    PROCESS_UTILITIES::Set_Backtrace(true);
    bool type_double=false; // float by default
    if(PARSE_ARGS::Find_And_Remove("-float",argc,argv)) type_double=false;
    if(PARSE_ARGS::Find_And_Remove("-double",argc,argv)) type_double=true;

    ANIMATED_VISUALIZATION* visualization=0;
    if(!type_double) visualization=new OPENGL_1D_VISUALIZATION<float>;
#ifndef COMPILE_WITHOUT_DOUBLE_SUPPORT
    else visualization=new OPENGL_1D_VISUALIZATION<double>;
#else
    else{std::cerr<<"Double support not enabled."<<std::endl;exit(1);}
#endif
    visualization->Initialize_And_Run(argc,argv);

    return 0;
}
