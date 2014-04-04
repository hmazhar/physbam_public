//#####################################################################
// Copyright 2003-2009, Zhaosheng Bao, Jon Gretarsson, Eran Guendelman, Geoffrey Irving, Nipun Kwatra, Michael Lentine, Frank Losasso, Nick Rasmussen, Avi Robinson-Mosher, Craig Schroeder, Andrew Selle, Andrew Selle, Tamar Shinar, Jonathan Su, Jerry Talton, Joseph Teran, Rachel Weinstein.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#include <PhysBAM_Tools/Grids_RLE/RLE_GRID_ITERATOR_CELL_2D.h>
#include <PhysBAM_Tools/Grids_Uniform_Interpolation/LINEAR_INTERPOLATION_UNIFORM.h>
#include <PhysBAM_Tools/Parsing/PARSE_ARGS.h>
#include <PhysBAM_Tools/Read_Write/Grids_Uniform/READ_WRITE_GRID.h>
#include <PhysBAM_Tools/Read_Write/Grids_Uniform_Arrays/READ_WRITE_ARRAYS.h>
#include <PhysBAM_Tools/Read_Write/Utilities/FILE_UTILITIES.h>
#include <PhysBAM_Tools/Utilities/PROCESS_UTILITIES.h>
#include <PhysBAM_Geometry/Basic_Geometry/BASIC_GEOMETRY_FORWARD.h>
#include <PhysBAM_Geometry/Geometry_Particles/REGISTER_GEOMETRY_READ_WRITE.h>
#include <PhysBAM_Geometry/Read_Write/Grids_Uniform_Level_Sets/READ_WRITE_LEVELSET_2D.h>
#include <PhysBAM_Rendering/PhysBAM_OpenGL/OpenGL/ANIMATED_VISUALIZATION.h>
#include <PhysBAM_Rendering/PhysBAM_OpenGL/OpenGL/OPENGL_COLOR_RAMP.h>
#include <PhysBAM_Rendering/PhysBAM_OpenGL/OpenGL/OPENGL_CONSTANT_COLOR_MAP.h>
#include <PhysBAM_Rendering/PhysBAM_OpenGL/OpenGL/OPENGL_GRID_2D.h>
#include <PhysBAM_Rendering/PhysBAM_OpenGL/OpenGL/OPENGL_INDEXED_COLOR_MAP.h>
#include <PhysBAM_Rendering/PhysBAM_OpenGL/OpenGL/OPENGL_LEVELSET_COLOR_MAP.h>
#include <PhysBAM_Rendering/PhysBAM_OpenGL/OpenGL/OPENGL_SEGMENTED_CURVE_2D.h>
#include <PhysBAM_Rendering/PhysBAM_OpenGL/OpenGL/OPENGL_WORLD.h>
#include <PhysBAM_Rendering/PhysBAM_OpenGL/OpenGL_Components/OPENGL_COMPONENT_ADAPTIVE_NODE_SCALAR_FIELD.h>
#include <PhysBAM_Rendering/PhysBAM_OpenGL/OpenGL_Components/OPENGL_COMPONENT_BASIC.h>
#include <PhysBAM_Rendering/PhysBAM_OpenGL/OpenGL_Components/OPENGL_COMPONENT_DEBUG_PARTICLES_2D.h>
#include <PhysBAM_Rendering/PhysBAM_OpenGL/OpenGL_Components/OPENGL_COMPONENT_DEFORMABLE_GEOMETRY_COLLECTION_2D.h>
#include <PhysBAM_Rendering/PhysBAM_OpenGL/OpenGL_Components/OPENGL_COMPONENT_FACE_SCALAR_FIELD_2D.h>
#include <PhysBAM_Rendering/PhysBAM_OpenGL/OpenGL_Components/OPENGL_COMPONENT_GRID_BASED_VECTOR_FIELD_2D.h>
#include <PhysBAM_Rendering/PhysBAM_OpenGL/OpenGL_Components/OPENGL_COMPONENT_LEVELSET_2D.h>
#include <PhysBAM_Rendering/PhysBAM_OpenGL/OpenGL_Components/OPENGL_COMPONENT_MAC_VELOCITY_FIELD_2D.h>
#include <PhysBAM_Rendering/PhysBAM_OpenGL/OpenGL_Components/OPENGL_COMPONENT_PARTICLES_2D.h>
#include <PhysBAM_Rendering/PhysBAM_OpenGL/OpenGL_Components/OPENGL_COMPONENT_PARTICLES_2D_DEFINITIONS.h>
#include <PhysBAM_Rendering/PhysBAM_OpenGL/OpenGL_Components/OPENGL_COMPONENT_QUADTREE_CELL_SCALAR_FIELD.h>
#include <PhysBAM_Rendering/PhysBAM_OpenGL/OpenGL_Components/OPENGL_COMPONENT_QUADTREE_FACE_SCALAR_FIELD.h>
#include <PhysBAM_Rendering/PhysBAM_OpenGL/OpenGL_Components/OPENGL_COMPONENT_QUADTREE_GRID.h>
#include <PhysBAM_Rendering/PhysBAM_OpenGL/OpenGL_Components/OPENGL_COMPONENT_QUADTREE_NODE_VECTOR_FIELD.h>
#include <PhysBAM_Rendering/PhysBAM_OpenGL/OpenGL_Components/OPENGL_COMPONENT_RIGID_GEOMETRY_COLLECTION_2D.h>
#include <PhysBAM_Rendering/PhysBAM_OpenGL/OpenGL_Components/OPENGL_COMPONENT_RLE_CELL_SCALAR_FIELD_2D.h>
#include <PhysBAM_Rendering/PhysBAM_OpenGL/OpenGL_Components/OPENGL_COMPONENT_RLE_FACE_SCALAR_FIELD_2D.h>
#include <PhysBAM_Rendering/PhysBAM_OpenGL/OpenGL_Components/OPENGL_COMPONENT_RLE_GRID_2D.h>
#include <PhysBAM_Rendering/PhysBAM_OpenGL/OpenGL_Components/OPENGL_COMPONENT_SCALAR_FIELD_2D.h>
#include <PhysBAM_Rendering/PhysBAM_OpenGL/OpenGL_Components/OPENGL_COMPONENT_SYMMETRIC_MATRIX_FIELD_2D.h>
#include <PhysBAM_Rendering/PhysBAM_OpenGL/OpenGL_Components/OPENGL_COMPONENT_THIN_SHELLS_DEBUGGING_2D.h>
#include <PhysBAM_Rendering/PhysBAM_OpenGL/OpenGL_Components/OPENGL_COMPONENT_TRIANGULATED_AREA.h>
#include <PhysBAM_Rendering/PhysBAM_OpenGL/OpenGL_Components/OPENGL_COMPONENT_TWO_PHASE_VELOCITY_MAGNITUDE_2D.h>
#include <fstream>
#include <sstream>
#include <string.h>
using namespace PhysBAM;
using namespace std;
//#####################################################################
// OPENGL_2D_VISUALIZATION
//#####################################################################
template<class T,class RW=T>
class OPENGL_2D_VISUALIZATION : public ANIMATED_VISUALIZATION
{
    typedef VECTOR<T,2> TV;
    typedef VECTOR<int,2> TV_INT;
public:
    OPENGL_2D_VISUALIZATION();
    ~OPENGL_2D_VISUALIZATION();

protected:
    virtual void Add_Arguments(PARSE_ARGS& parse_args);
    virtual void Parse_Arguments(PARSE_ARGS& parse_args);
    virtual void Initialize_Components_And_Key_Bindings();
    virtual void Add_OpenGL_Initialization();

private:
    void Read_Grid();
    virtual void Pre_Frame_Extra();
    virtual void Set_Frame_Extra();

    void Cycle_Draw_Extra();
    void Command_Prompt();
    void Command_Prompt_Response();
    void Toggle_2d_Mode();
    void Change_Level(int new_level);
    DEFINE_CALLBACK_CREATOR(OPENGL_2D_VISUALIZATION,Cycle_Draw_Extra);
    DEFINE_CALLBACK_CREATOR(OPENGL_2D_VISUALIZATION,Command_Prompt);
    DEFINE_CALLBACK_CREATOR(OPENGL_2D_VISUALIZATION,Command_Prompt_Response);
    DEFINE_CALLBACK_CREATOR(OPENGL_2D_VISUALIZATION,Toggle_2d_Mode);

    // Components
    OPENGL_COMPONENT_SCALAR_FIELD_2D<T> *sph_cell_weight_component;
    OPENGL_COMPONENT_SCALAR_FIELD_2D<T>* pressure_component,* coarse_pressure_component;
    OPENGL_COMPONENT_BASIC<OPENGL_GRID_2D<T> >* grid_component,*coarse_grid_component;
#ifndef COMPILE_WITHOUT_DYADIC_SUPPORT
    OPENGL_COMPONENT_QUADTREE_GRID<T>* quadtree_grid_component;
#endif
#ifndef COMPILE_WITHOUT_RLE_SUPPORT
    OPENGL_COMPONENT_RLE_GRID_2D<T>* rle_grid_component;
#endif
    // Options
    std::string basedir;

    GRID<TV> grid,mac_grid,regular_grid,coarse_grid,coarse_mac_grid,coarse_regular_grid;
    bool has_valid_grid;
#ifndef COMPILE_WITHOUT_RLE_SUPPORT
    bool has_valid_rle_grid;
#endif
#ifndef COMPILE_WITHOUT_DYADIC_SUPPORT
    bool has_valid_quadtree_grid;
#endif
    bool has_valid_coarse_grid;
    bool node_based,coarse_node_based;

    ARRAY<int> rigid_bodies_no_draw_list;
};
//#####################################################################
// OPENGL_2D_VISUALIZATION
//#####################################################################
template<class T,class RW> OPENGL_2D_VISUALIZATION<T,RW>::
OPENGL_2D_VISUALIZATION()
    :ANIMATED_VISUALIZATION(),
    sph_cell_weight_component(0),
    pressure_component(0),coarse_pressure_component(0),grid_component(0),coarse_grid_component(0)
#ifndef COMPILE_WITHOUT_DYADIC_SUPPORT
    ,quadtree_grid_component(0)
#endif
#ifndef COMPILE_WITHOUT_RLE_SUPPORT
    ,rle_grid_component(0)
#endif
{
    add_axes=false;
}
//#####################################################################
// OPENGL_2D_VISUALIZATION
//#####################################################################
template<class T,class RW> OPENGL_2D_VISUALIZATION<T,RW>::
~OPENGL_2D_VISUALIZATION()
{
    delete grid_component;
#ifndef COMPILE_WITHOUT_RLE_SUPPORT
    delete rle_grid_component;
#endif
#ifndef COMPILE_WITHOUT_DYADIC_SUPPORT
    delete quadtree_grid_component;
#endif
    delete coarse_grid_component;
}
//#####################################################################
// Add_Arguments
//#####################################################################
template<class T,class RW> void OPENGL_2D_VISUALIZATION<T,RW>::
Add_Arguments(PARSE_ARGS& parse_args)
{
    basedir=".";

    ANIMATED_VISUALIZATION::Add_Arguments(parse_args);

    parse_args.Add_Option_Argument("-float");
    parse_args.Add_Option_Argument("-double");
    parse_args.Add_Option_Argument("-no_ghost");
    parse_args.Add_String_Argument("-rigid_bodies_no_draw","");
    parse_args.Set_Extra_Arguments(-1,"[<basedir>]");
}
//#####################################################################
// Parse_Arguments
//#####################################################################
template<class T,class RW> void OPENGL_2D_VISUALIZATION<T,RW>::
Parse_Arguments(PARSE_ARGS& parse_args)
{
    ANIMATED_VISUALIZATION::Parse_Arguments(parse_args);

    if(parse_args.Num_Extra_Args()>1)
        parse_args.Print_Usage(true);
    else if(parse_args.Num_Extra_Args()==1)
        basedir=parse_args.Extra_Arg(1);

    if(parse_args.Is_Value_Set("-rigid_bodies_no_draw")){ARRAY<int> int_list;
        STRING_UTILITIES::Parse_Integer_List(parse_args.Get_String_Value("-rigid_bodies_no_draw"),int_list);
        rigid_bodies_no_draw_list.Resize(int_list.Size());
        for(int i=1;i<=rigid_bodies_no_draw_list.m;i++) rigid_bodies_no_draw_list(i)=int(int_list(i));}

    last_frame_filename=basedir+"/common/last_frame";

    if(FILE_UTILITIES::File_Exists(basedir+"/common/first_frame") && !parse_args.Is_Value_Set("-start_frame")) FILE_UTILITIES::Read_From_Text_File(basedir+"/common/first_frame",start_frame);

    // don't override camera script filename if it was already set in base class based on command line argument
    if(camera_script_filename.empty()) camera_script_filename=basedir+"/camera_script";

#ifdef __linux__
    if(!parse_args.Is_Value_Set("-window_title")) opengl_window_title="opengl_2d: " + FILE_UTILITIES::Real_Path(basedir);
#endif
}
//#####################################################################
// Read_Grid
//#####################################################################
template<class T,class RW> void OPENGL_2D_VISUALIZATION<T,RW>::
Read_Grid()
{
    has_valid_grid=false;
    has_valid_coarse_grid=false;
    std::string filename,coarse_filename;

    filename=STRING_UTILITIES::string_sprintf("%s/%d/levelset",basedir.c_str(),start_frame);
    // For backwards compatibility
    if(!FILE_UTILITIES::File_Exists(filename)) filename=STRING_UTILITIES::string_sprintf("%s/%d/levelset.phi",basedir.c_str(),start_frame);
    if(!FILE_UTILITIES::File_Exists(filename)) filename=STRING_UTILITIES::string_sprintf("%s/%d/levelset_1",basedir.c_str(),start_frame);

#ifndef COMPILE_WITHOUT_RLE_SUPPORT
    has_valid_rle_grid=false;
    std::string rle_filename=STRING_UTILITIES::string_sprintf("%s/%d/rle_grid",basedir.c_str(),start_frame);
    if(FILE_UTILITIES::File_Exists(rle_filename)){
        LOG::cout<<"Reading rle grid '"<<rle_filename<<"'"<<std::endl<<std::flush;
        rle_filename=STRING_UTILITIES::string_sprintf("%s/rle_grid.%%d",basedir.c_str());
        rle_grid_component=new OPENGL_COMPONENT_RLE_GRID_2D<T>(rle_filename);
        // Add_Component(rle_grid_component); // add later instead of here, but Pre_Frame_Extra will ensure rle grid is read in
        has_valid_rle_grid=true;}
#endif
#ifndef COMPILE_WITHOUT_DYADIC_SUPPORT
    has_valid_quadtree_grid=false;
    std::string quadtree_filename=STRING_UTILITIES::string_sprintf("%s/%d/quadtree_grid",basedir.c_str(),start_frame);
    if(FILE_UTILITIES::File_Exists(quadtree_filename)){
        LOG::cout<<"Reading quadtree '"<<quadtree_filename<<"'"<<std::endl<<std::flush;
        quadtree_filename=STRING_UTILITIES::string_sprintf("%s/quadtree_grid.%%d",basedir.c_str());
        quadtree_grid_component=new OPENGL_COMPONENT_QUADTREE_GRID<T>(quadtree_filename);
        // Add_Component(quadtree_grid_component); // add later instead of here, but Pre_Frame_Extra will ensure quadtree grid is read in
        has_valid_quadtree_grid=true;}
#endif

    if(FILE_UTILITIES::File_Exists(basedir+"/common/coarse_grid")){
        coarse_filename=basedir+"/common/coarse_grid";
        LOG::cout<<"Reading coarse grid from '"<<coarse_filename<<"'..."<<std::flush;
        FILE_UTILITIES::Read_From_File<RW>(coarse_filename,coarse_grid);
        has_valid_coarse_grid=true;}

    if(FILE_UTILITIES::File_Exists(filename)){
        LOG::cout<<"Reading grid from '"<<filename<<"'..."<<std::endl<<std::flush;
        ARRAY<T,VECTOR<int,2> > phi;
        LEVELSET_2D<GRID<TV> > levelset(grid,phi);
        FILE_UTILITIES::Read_From_File<RW>(filename,levelset);
        has_valid_grid=true;}
    else if(FILE_UTILITIES::File_Exists(basedir+"/common/grid")){
        filename=basedir+"/common/grid";
        LOG::cout<<"Reading grid from '"<<filename<<"'..."<<std::flush;
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
template<class T,class RW> void OPENGL_2D_VISUALIZATION<T,RW>::
Initialize_Components_And_Key_Bindings()
{
    ANIMATED_VISUALIZATION::Initialize_Components_And_Key_Bindings();
    //opengl_world.Unbind_Keys("abcdDeilLMnoOStTuvVZ 1!2@3#4$5%67&*890 ;'~\t=-`\\^[]");
    Read_Grid();
    std::string filename,filename2;
    opengl_world.Set_Key_Binding_Category_Priority(1);

    // Density
    filename=basedir+"/%d/density";
    if(has_valid_grid && FILE_UTILITIES::Frame_File_Exists(filename,start_frame)){
        OPENGL_COMPONENT_SCALAR_FIELD_2D<T>* density_component=new OPENGL_COMPONENT_SCALAR_FIELD_2D<T>(grid,filename,OPENGL_COLOR_RAMP<T>::Two_Color_Ramp(0,1,OPENGL_COLOR::Gray(0,0),OPENGL_COLOR::Gray(1,1)));
        density_component->opengl_scalar_field.Set_Uniform_Contour_Values(.2568,6.067,(6.067-.2568)/30.);
        density_component->opengl_scalar_field.Update();
        opengl_world.Set_Key_Binding_Category("Density");
        Add_Component(density_component,"Density",'d',BASIC_VISUALIZATION::OWNED);
        opengl_world.Append_Bind_Key('D',density_component->Toggle_Color_Map_CB());
        opengl_world.Append_Bind_Key(OPENGL_KEY(OPENGL_KEY::F8),density_component->Toggle_Draw_Mode_CB());
        opengl_world.Append_Bind_Key('`',density_component->Toggle_Smooth_CB());}
    filename=basedir+"/%d/density_gradient";
    if(has_valid_grid && FILE_UTILITIES::Frame_File_Exists(filename,start_frame)){
        OPENGL_COMPONENT_SCALAR_FIELD_2D<T>* density_gradient_component=new OPENGL_COMPONENT_SCALAR_FIELD_2D<T>(grid,filename,
            OPENGL_COLOR_RAMP<T>::Two_Color_Ramp(0,1,OPENGL_COLOR::Gray(0,0),OPENGL_COLOR::Gray(1,1)));
        density_gradient_component->opengl_scalar_field.Update();
        density_gradient_component->Toggle_Draw();
        opengl_world.Set_Key_Binding_Category("Density");
        Add_Component(density_gradient_component,"Density Gradient",'8',BASIC_VISUALIZATION::OWNED);}
    // Soot
    filename=basedir+"/%d/soot";
    if(has_valid_grid && FILE_UTILITIES::Frame_File_Exists(filename,start_frame)){
        OPENGL_COMPONENT_SCALAR_FIELD_2D<T>* soot_component=new OPENGL_COMPONENT_SCALAR_FIELD_2D<T>(grid,filename,OPENGL_COLOR_RAMP<T>::Two_Color_Ramp(0,.01,OPENGL_COLOR::Gray(0,0),OPENGL_COLOR::Gray(1,1)));
        soot_component->opengl_scalar_field.Set_Uniform_Contour_Values(0,2,.1);
        soot_component->opengl_scalar_field.Update();
        opengl_world.Set_Key_Binding_Category("Soot");
        Add_Component(soot_component,"Soot",'i',BASIC_VISUALIZATION::OWNED|BASIC_VISUALIZATION::START_HIDDEN);
        opengl_world.Append_Bind_Key('D',soot_component->Toggle_Color_Map_CB());
        opengl_world.Append_Bind_Key(OPENGL_KEY(OPENGL_KEY::F8),soot_component->Toggle_Draw_Mode_CB());
        opengl_world.Append_Bind_Key('`',soot_component->Toggle_Smooth_CB());}
#ifndef COMPILE_WITHOUT_DYADIC_SUPPORT
    filename=basedir+"/%d/quadtree_density";
    if(has_valid_quadtree_grid && FILE_UTILITIES::Frame_File_Exists(filename,start_frame)){
        OPENGL_COMPONENT_ADAPTIVE_NODE_SCALAR_FIELD<QUADTREE_GRID<T> >* quadtree_density_component=new OPENGL_COMPONENT_ADAPTIVE_NODE_SCALAR_FIELD<QUADTREE_GRID<T> >(quadtree_grid_component->opengl_grid.grid,filename,
            OPENGL_COLOR_RAMP<T>::Two_Color_Ramp(0,1,OPENGL_COLOR::Gray(0,0),OPENGL_COLOR::Gray(1,1)));
        opengl_world.Set_Key_Binding_Category("Density");
        Add_Component(quadtree_density_component,"Quadtree Density",'d',BASIC_VISUALIZATION::OWNED);
        opengl_world.Append_Bind_Key('+',quadtree_density_component->Increase_Point_Size_CB());
        opengl_world.Append_Bind_Key('_',quadtree_density_component->Decrease_Point_Size_CB());}
#endif
    filename=basedir+"/%d/internal_energy";
    if(has_valid_grid && FILE_UTILITIES::Frame_File_Exists(filename,start_frame)){
        OPENGL_COMPONENT_SCALAR_FIELD_2D<T>* internal_energy_component=new OPENGL_COMPONENT_SCALAR_FIELD_2D<T>(grid,filename,
            OPENGL_COLOR_RAMP<T>::Two_Color_Ramp(0,1,OPENGL_COLOR::Gray(0,0),OPENGL_COLOR::Yellow(1,1)));
        internal_energy_component->opengl_scalar_field.Update();
        internal_energy_component->Toggle_Draw();
        opengl_world.Set_Key_Binding_Category("Density");
        Add_Component(internal_energy_component,"Internal_Energy",'E',BASIC_VISUALIZATION::OWNED);}

    // Temperature
    filename=basedir+"/%d/temperature";
    if(has_valid_grid && FILE_UTILITIES::Frame_File_Exists(filename,start_frame)){
        ARRAY<T,VECTOR<int,2> > temp0;FILE_UTILITIES::Read_From_File<RW>(FILE_UTILITIES::Get_Frame_Filename(filename,start_frame),temp0);
        OPENGL_COMPONENT_SCALAR_FIELD_2D<T>* temperature_component=new OPENGL_COMPONENT_SCALAR_FIELD_2D<T>(grid,filename,OPENGL_COLOR_RAMP<T>::Two_Color_Ramp(0,1,OPENGL_COLOR::Gray(0,0),OPENGL_COLOR::Red(1,1)));
        temperature_component->opengl_scalar_field.Set_Scale_Range(temp0.Min(),temp0.Max());
        opengl_world.Set_Key_Binding_Category("Temperature");
        Add_Component(temperature_component,"Temperature",'t',BASIC_VISUALIZATION::OWNED|BASIC_VISUALIZATION::START_HIDDEN);
        opengl_world.Append_Bind_Key('T',temperature_component->Toggle_Color_Map_CB());
        opengl_world.Append_Bind_Key('`',temperature_component->Toggle_Smooth_CB());}
#ifndef COMPILE_WITHOUT_DYADIC_SUPPORT
    filename=basedir+"/%d/quadtree_temperature";
    if(has_valid_quadtree_grid && FILE_UTILITIES::Frame_File_Exists(filename,start_frame)){
        ARRAY<T> temp0;FILE_UTILITIES::Read_From_File<RW>(FILE_UTILITIES::Get_Frame_Filename(filename,start_frame),temp0);
        // TODO: put scaling back into quadtree
        //T min_temp=temp0.Min(),max_temp=temp0.Max();
        OPENGL_COMPONENT_ADAPTIVE_NODE_SCALAR_FIELD<QUADTREE_GRID<T> >* quadtree_temperature_component=new OPENGL_COMPONENT_ADAPTIVE_NODE_SCALAR_FIELD<QUADTREE_GRID<T> >(
            quadtree_grid_component->opengl_grid.grid,filename,OPENGL_COLOR_RAMP<T>::Two_Color_Ramp(0,1,OPENGL_COLOR::Gray(0,0),OPENGL_COLOR::Red(1,1)));
        opengl_world.Set_Key_Binding_Category("Quadtree Temperature");
        Add_Component(quadtree_temperature_component,"Quadtree Temperature",'t',BASIC_VISUALIZATION::OWNED|BASIC_VISUALIZATION::START_HIDDEN);
        opengl_world.Append_Bind_Key('+',quadtree_temperature_component->Increase_Point_Size_CB());
        opengl_world.Append_Bind_Key('_',quadtree_temperature_component->Decrease_Point_Size_CB());}
#endif

    // SPH
    filename=basedir+"/%d/sph_cell_weights";
    if(has_valid_grid && FILE_UTILITIES::Frame_File_Exists(filename,start_frame)){
        sph_cell_weight_component=new OPENGL_COMPONENT_SCALAR_FIELD_2D<T>(grid,filename,OPENGL_COLOR_RAMP<T>::Two_Color_Ramp(0,100,OPENGL_COLOR::Cyan(0,0),OPENGL_COLOR::Cyan(1)));
        sph_cell_weight_component->opengl_scalar_field.Update();
        Add_Component(sph_cell_weight_component,"SPH Cell Weights",'\0',BASIC_VISUALIZATION::OWNED);}

    // Grid vorticity and rasterized vortex particles
    const T vorticity_visualization_min=-10,vorticity_visualization_max=10;
    filename=basedir+"/%d/grid_vorticity";
    if(has_valid_grid && FILE_UTILITIES::Frame_File_Exists(filename,start_frame)){
        opengl_world.Set_Key_Binding_Category("Vorticity");
        Add_Component(new OPENGL_COMPONENT_SCALAR_FIELD_2D<T>(grid,filename,OPENGL_COLOR_RAMP<T>::Matlab_Jet(vorticity_visualization_min,vorticity_visualization_max))
            ,"Grid Vorticity",'u',BASIC_VISUALIZATION::OWNED|BASIC_VISUALIZATION::START_HIDDEN);}
    filename=basedir+"/%d/grid_vorticity_raw_particles";
    if(has_valid_grid && FILE_UTILITIES::Frame_File_Exists(filename,start_frame)){
        opengl_world.Set_Key_Binding_Category("Vorticity");
        Add_Component(new OPENGL_COMPONENT_SCALAR_FIELD_2D<T>(grid,filename,OPENGL_COLOR_RAMP<T>::Matlab_Jet(vorticity_visualization_min,vorticity_visualization_max)),
            "Grid Vorticity Raw Particles",'i',BASIC_VISUALIZATION::OWNED|BASIC_VISUALIZATION::START_HIDDEN);}
    filename=basedir+"/%d/grid_vorticity_particles";
    if(has_valid_grid && FILE_UTILITIES::Frame_File_Exists(filename,start_frame)){
        opengl_world.Set_Key_Binding_Category("Vorticity");
        Add_Component(new OPENGL_COMPONENT_SCALAR_FIELD_2D<T>(grid,filename,OPENGL_COLOR_RAMP<T>::Matlab_Jet(vorticity_visualization_min,vorticity_visualization_max)),
            "Grid Vorticity Particles",'o',BASIC_VISUALIZATION::OWNED|BASIC_VISUALIZATION::START_HIDDEN);}

    // Level sets
    filename=basedir+"/%d/levelset";
    OPENGL_COMPONENT_LEVELSET_2D<T>* levelset_component=0;
    if(!FILE_UTILITIES::Frame_File_Exists(filename,start_frame) && !FILE_UTILITIES::Frame_File_Exists(basedir+"/%d/levelset_1",start_frame)) filename=basedir+"/%d/levelset.phi"; // for backwards compatiblity
    if(FILE_UTILITIES::Frame_File_Exists(filename,start_frame) || FILE_UTILITIES::Frame_File_Exists(basedir+"/%d/levelset_1",start_frame)){
        levelset_component=new OPENGL_COMPONENT_LEVELSET_2D<T>(filename,basedir+"/%d/levelset_%d");
        if(levelset_component->opengl_levelsets.m>1) for(int j=1;j<=levelset_component->opengl_levelsets.m;j++){
            levelset_component->opengl_levelsets(j)->draw_cells=true;
            levelset_component->opengl_levelsets(j)->draw_area=false;
            levelset_component->opengl_levelsets(j)->draw_curve=true;
            levelset_component->opengl_levelsets(j)->Update();}
        else{
            levelset_component->opengl_levelsets(1)->Set_Inside_And_Outside_Colors(OPENGL_COLOR::Blue(),OPENGL_COLOR::Red(.5));
            levelset_component->opengl_levelsets(1)->draw_cells=true;
            levelset_component->opengl_levelsets(1)->draw_area=false;
            levelset_component->opengl_levelsets(1)->draw_curve=false;
            levelset_component->opengl_levelsets(1)->Update();}
        opengl_world.Set_Key_Binding_Category("Level Set");
        Add_Component(levelset_component,"Level Set",'l',BASIC_VISUALIZATION::OWNED);
        opengl_world.Append_Bind_Key(';',levelset_component->Toggle_Draw_Mode_CB());
        opengl_world.Append_Bind_Key('\'',levelset_component->Toggle_Draw_Sign_CB());
        opengl_world.Append_Bind_Key('L',levelset_component->Toggle_Color_Mode_CB());
        opengl_world.Append_Bind_Key('n',levelset_component->Toggle_Normals_CB());
        opengl_world.Append_Bind_Key('`',levelset_component->Toggle_Smooth_CB());
        opengl_world.Append_Bind_Key('>',levelset_component->Next_Set_CB());
        opengl_world.Append_Bind_Key('<',levelset_component->Previous_Set_CB());
        opengl_world.Append_Bind_Key('M',levelset_component->Toggle_Draw_Multiple_Levelsets_CB());
        opengl_world.Append_Bind_Key('^',levelset_component->Toggle_Draw_Ghost_Values_CB());}
    filename=basedir+"/%d/object_levelset";
    if(FILE_UTILITIES::Frame_File_Exists(filename,start_frame)){
        OPENGL_COMPONENT_LEVELSET_2D<T>* object_levelset_component=new OPENGL_COMPONENT_LEVELSET_2D<T>(filename);
        object_levelset_component->opengl_levelset->outside_color=OPENGL_COLOR::Red(.6);
        object_levelset_component->opengl_levelset->inside_color=OPENGL_COLOR::Blue(.6);
        object_levelset_component->opengl_levelset->draw_cells=true;
        object_levelset_component->opengl_levelset->draw_area=false;
        object_levelset_component->opengl_levelset->draw_curve=false;
        object_levelset_component->opengl_levelset->Update();
        opengl_world.Set_Key_Binding_Category("Object Level Set");
        Add_Component(object_levelset_component,"Object Level Set",'9',BASIC_VISUALIZATION::OWNED|BASIC_VISUALIZATION::START_HIDDEN);
        opengl_world.Append_Bind_Key('(',object_levelset_component->Toggle_Draw_Mode_CB());
        opengl_world.Append_Bind_Key(')',object_levelset_component->Toggle_Color_Mode_CB());}
#ifndef COMPILE_WITHOUT_RLE_SUPPORT
    filename=basedir+"/%d/rle_levelset";
    if(has_valid_rle_grid && FILE_UTILITIES::Frame_File_Exists(filename,start_frame)){
        OPENGL_COMPONENT_RLE_CELL_SCALAR_FIELD_2D<T>* rle_levelset_component=new OPENGL_COMPONENT_RLE_CELL_SCALAR_FIELD_2D<T>(rle_grid_component->opengl_grid.grid,filename,
            new OPENGL_LEVELSET_COLOR_MAP<T>(OPENGL_COLOR::Blue(),OPENGL_COLOR::Red()));
        rle_levelset_component->opengl_rle_cell_scalar_field.draw_filled=true;
        opengl_world.Set_Key_Binding_Category("RLE Level Set");
        Add_Component(rle_levelset_component,"RLE Level Set",'l',BASIC_VISUALIZATION::OWNED);}
    filename=basedir+"/%d/rle_object_levelset";
    if(has_valid_rle_grid && FILE_UTILITIES::Frame_File_Exists(filename,start_frame)){
        OPENGL_COMPONENT_RLE_CELL_SCALAR_FIELD_2D<T>* rle_object_levelset_component=new OPENGL_COMPONENT_RLE_CELL_SCALAR_FIELD_2D<T>(rle_grid_component->opengl_grid.grid,filename,
            new OPENGL_LEVELSET_COLOR_MAP<T>(OPENGL_COLOR::Blue(),OPENGL_COLOR::Red()));
        rle_object_levelset_component->opengl_rle_cell_scalar_field.draw_filled=true;
        opengl_world.Set_Key_Binding_Category("RLE Object Level Set");
        Add_Component(rle_object_levelset_component,"RLE Object Level Set",'9',BASIC_VISUALIZATION::OWNED|BASIC_VISUALIZATION::START_HIDDEN);}
#endif
#ifndef COMPILE_WITHOUT_DYADIC_SUPPORT
    filename=basedir+"/%d/quadtree_levelset";
    if(has_valid_quadtree_grid && FILE_UTILITIES::Frame_File_Exists(filename,start_frame)){
        OPENGL_COMPONENT_QUADTREE_CELL_SCALAR_FIELD<T>* quadtree_levelset_component=new OPENGL_COMPONENT_QUADTREE_CELL_SCALAR_FIELD<T>(quadtree_grid_component->opengl_grid.grid,filename,
            new OPENGL_LEVELSET_COLOR_MAP<T>(OPENGL_COLOR::Blue(),OPENGL_COLOR::Red()));
        opengl_world.Set_Key_Binding_Category("Quadtree Level Set");
        Add_Component(quadtree_levelset_component,"Quadtree Level Set",'l',BASIC_VISUALIZATION::OWNED);
        opengl_world.Append_Bind_Key('+',quadtree_levelset_component->Increase_Point_Size_CB());
        opengl_world.Append_Bind_Key('_',quadtree_levelset_component->Decrease_Point_Size_CB());}
    filename=basedir+"/%d/quadtree_levelset_nodes";
    if(has_valid_quadtree_grid && FILE_UTILITIES::Frame_File_Exists(filename,start_frame)){
        OPENGL_COMPONENT_ADAPTIVE_NODE_SCALAR_FIELD<QUADTREE_GRID<T> >* quadtree_levelset_nodes_component=new OPENGL_COMPONENT_ADAPTIVE_NODE_SCALAR_FIELD<QUADTREE_GRID<T> >(quadtree_grid_component->opengl_grid.grid,filename,
            new OPENGL_LEVELSET_COLOR_MAP<T>(OPENGL_COLOR::Blue(),OPENGL_COLOR::Red()));
        opengl_world.Set_Key_Binding_Category("Quadtree Nodal Level Set");
        Add_Component(quadtree_levelset_nodes_component,"Quadtree Levelset",'L',BASIC_VISUALIZATION::OWNED);
        opengl_world.Append_Bind_Key('+',quadtree_levelset_nodes_component->Increase_Point_Size_CB());
        opengl_world.Append_Bind_Key('_',quadtree_levelset_nodes_component->Decrease_Point_Size_CB());}
#endif

    // Particles
    opengl_world.Set_Key_Binding_Category("Particles");

    // Uniform velocities
    OPENGL_COMPONENT_MAC_VELOCITY_FIELD_2D<T>* mac_velocity_component=0;
    OPENGL_COMPONENT_GRID_BASED_VECTOR_FIELD_2D<T>* vector_velocity_component=0;
    opengl_world.Set_Key_Binding_Category("Velocity");
    filename=basedir+"/%d/velocities";
    if(has_valid_grid && FILE_UTILITIES::Frame_File_Exists(filename,start_frame)){
        if(node_based){
            vector_velocity_component=new OPENGL_COMPONENT_GRID_BASED_VECTOR_FIELD_2D<T>(regular_grid,filename);
            vector_velocity_component->opengl_grid_based_vector_field.size=.01;
            vector_velocity_component->opengl_grid_based_vector_field.vector_color=OPENGL_COLOR::Green();
            Add_Component(vector_velocity_component,"Node velocities",'\0',BASIC_VISUALIZATION::OWNED|BASIC_VISUALIZATION::START_HIDDEN);}
        else{
            mac_velocity_component=new OPENGL_COMPONENT_MAC_VELOCITY_FIELD_2D<T>(mac_grid,filename);
            mac_velocity_component->opengl_mac_velocity_field->size=.01;
            mac_velocity_component->opengl_mac_velocity_field->vector_color=OPENGL_COLOR::Green();
            Add_Component(mac_velocity_component,"MAC velocities",'\0',BASIC_VISUALIZATION::OWNED|BASIC_VISUALIZATION::START_HIDDEN);}}
    filename=basedir+"/%d/mac_velocities";
    if(has_valid_grid && FILE_UTILITIES::Frame_File_Exists(filename,start_frame)){
        mac_velocity_component=new OPENGL_COMPONENT_MAC_VELOCITY_FIELD_2D<T>(mac_grid,filename);
        for(int i=1;i<=mac_velocity_component->opengl_adaptive_mac_velocity_fields.m;i++){
            mac_velocity_component->opengl_adaptive_mac_velocity_fields(i)->size=.01;
            mac_velocity_component->opengl_adaptive_mac_velocity_fields(i)->vector_color=OPENGL_COLOR::Magenta();
            mac_velocity_component->opengl_adaptive_mac_velocity_fields(i)->Set_Velocity_Mode(OPENGL_MAC_VELOCITY_FIELD_2D<T>::FACE_CENTERED);}
        mac_velocity_component->Set_Psi_N_Psi_D_Basedir_For_Divergence(basedir);
        Add_Component(mac_velocity_component,"MAC velocities",'\0',BASIC_VISUALIZATION::OWNED|BASIC_VISUALIZATION::START_HIDDEN);}
    filename=basedir+"/%d/compressible_mac_velocities";
    if(has_valid_grid && FILE_UTILITIES::Frame_File_Exists(filename,start_frame)){
        mac_velocity_component=new OPENGL_COMPONENT_MAC_VELOCITY_FIELD_2D<T>(mac_grid,filename);
        for(int i=1;i<=mac_velocity_component->opengl_adaptive_mac_velocity_fields.m;i++){
            mac_velocity_component->opengl_adaptive_mac_velocity_fields(i)->size=.01;
            mac_velocity_component->opengl_adaptive_mac_velocity_fields(i)->vector_color=OPENGL_COLOR::Magenta();
            mac_velocity_component->opengl_adaptive_mac_velocity_fields(i)->Set_Velocity_Mode(OPENGL_MAC_VELOCITY_FIELD_2D<T>::FACE_CENTERED);}
        mac_velocity_component->Set_Psi_N_Psi_D_Basedir_For_Divergence(basedir);
        Add_Component(mac_velocity_component,"Compressible MAC velocities",'\0',BASIC_VISUALIZATION::OWNED|BASIC_VISUALIZATION::START_HIDDEN);}
    filename=basedir+"/%d/centered_velocities";
    if(has_valid_grid && FILE_UTILITIES::Frame_File_Exists(filename,start_frame)){
        vector_velocity_component=new OPENGL_COMPONENT_GRID_BASED_VECTOR_FIELD_2D<T>(mac_grid,filename);
        vector_velocity_component->opengl_grid_based_vector_field.size=.01;
        vector_velocity_component->opengl_grid_based_vector_field.vector_color=OPENGL_COLOR::Green();
        Add_Component(vector_velocity_component,"Centered velocities",'\0',BASIC_VISUALIZATION::OWNED|BASIC_VISUALIZATION::START_HIDDEN);}
    filename=basedir+"/%d/mass_fluxes";
    if(has_valid_grid && FILE_UTILITIES::Frame_File_Exists(filename,start_frame)){
        std::cout<<"Using mass fluxes (from "<<filename<<")"<<std::endl;
        OPENGL_COMPONENT_MAC_VELOCITY_FIELD_2D<T>* face_flux_component=new OPENGL_COMPONENT_MAC_VELOCITY_FIELD_2D<T>(mac_grid,filename);
        for(int i=1;i<=face_flux_component->opengl_adaptive_mac_velocity_fields.m;i++){
            face_flux_component->opengl_adaptive_mac_velocity_fields(i)->size=.01;
            face_flux_component->opengl_adaptive_mac_velocity_fields(i)->vector_color=OPENGL_COLOR::Magenta();
            face_flux_component->opengl_adaptive_mac_velocity_fields(i)->Set_Velocity_Mode(OPENGL_MAC_VELOCITY_FIELD_2D<T>::FACE_CENTERED);}
        face_flux_component->Set_Psi_N_Psi_D_Basedir_For_Divergence(basedir);
        Add_Component(face_flux_component,"Mass fluxes",'\0',BASIC_VISUALIZATION::OWNED|BASIC_VISUALIZATION::START_HIDDEN);
        opengl_world.Append_Bind_Key('M',face_flux_component->Toggle_Velocity_Mode_And_Draw_CB());
        opengl_world.Append_Bind_Key('=',face_flux_component->Increase_Vector_Size_CB());
        opengl_world.Append_Bind_Key('-',face_flux_component->Decrease_Vector_Size_CB());}
    // TODO: integrate these into creation stage...
    if(mac_velocity_component){
        opengl_world.Append_Bind_Key('C',mac_velocity_component->Toggle_Draw_CB());
        opengl_world.Append_Bind_Key('V',mac_velocity_component->Toggle_Velocity_Mode_And_Draw_CB());
        opengl_world.Append_Bind_Key("<F9>",mac_velocity_component->Toggle_Draw_Divergence_CB());
        opengl_world.Append_Bind_Key('=',mac_velocity_component->Increase_Vector_Size_CB());
        opengl_world.Append_Bind_Key('-',mac_velocity_component->Decrease_Vector_Size_CB());
        opengl_world.Append_Bind_Key('h',mac_velocity_component->Toggle_Arrowhead_CB());
        opengl_world.Append_Bind_Key('S',mac_velocity_component->Toggle_Draw_Streamlines_CB());
        opengl_world.Append_Bind_Key('/',mac_velocity_component->Toggle_Use_Streamline_Seed_CB());
        opengl_world.Append_Bind_Key(']',mac_velocity_component->Lengthen_Streamlines_CB());
        opengl_world.Append_Bind_Key('[',mac_velocity_component->Shorten_Streamlines_CB());
        opengl_world.Append_Bind_Key('v',mac_velocity_component->Toggle_Draw_Vorticity_CB());
        opengl_world.Append_Bind_Key('N',mac_velocity_component->Normalize_Vorticity_Color_Map_CB());}
    if(vector_velocity_component){
        opengl_world.Append_Bind_Key('v',vector_velocity_component->Toggle_Draw_CB());
        opengl_world.Append_Bind_Key('=',vector_velocity_component->Increase_Vector_Size_CB());
        opengl_world.Append_Bind_Key('-',vector_velocity_component->Decrease_Vector_Size_CB());
        opengl_world.Append_Bind_Key('h',vector_velocity_component->Toggle_Arrowhead_CB());}

    // Uniform coarse velocities
    OPENGL_COMPONENT_MAC_VELOCITY_FIELD_2D<T>* coarse_mac_velocity_component=0;
    opengl_world.Set_Key_Binding_Category("Coarse Velocity");
    filename=basedir+"/%d/coarse_mac_velocities";
    if(has_valid_coarse_grid && FILE_UTILITIES::Frame_File_Exists(filename,start_frame)){
        coarse_mac_velocity_component=new OPENGL_COMPONENT_MAC_VELOCITY_FIELD_2D<T>(coarse_mac_grid,filename);
        for(int i=1;i<=coarse_mac_velocity_component->opengl_adaptive_mac_velocity_fields.m;i++){
            coarse_mac_velocity_component->opengl_adaptive_mac_velocity_fields(i)->size=.01;
            coarse_mac_velocity_component->opengl_adaptive_mac_velocity_fields(i)->vector_color=OPENGL_COLOR::Blue();
            coarse_mac_velocity_component->opengl_adaptive_mac_velocity_fields(i)->Set_Velocity_Mode(OPENGL_MAC_VELOCITY_FIELD_2D<T>::FACE_CENTERED);}
        coarse_mac_velocity_component->Set_Psi_N_Psi_D_Basedir_For_Divergence(basedir);
        Add_Component(coarse_mac_velocity_component,"coarse MAC velocities",'\0',BASIC_VISUALIZATION::OWNED|BASIC_VISUALIZATION::START_HIDDEN);}
    if(coarse_mac_velocity_component){
        opengl_world.Append_Bind_Key('x',coarse_mac_velocity_component->Toggle_Draw_CB());
        opengl_world.Append_Bind_Key('V',coarse_mac_velocity_component->Toggle_Velocity_Mode_And_Draw_CB());
        opengl_world.Append_Bind_Key("<F9>",coarse_mac_velocity_component->Toggle_Draw_Divergence_CB());
        opengl_world.Append_Bind_Key('=',coarse_mac_velocity_component->Increase_Vector_Size_CB());
        opengl_world.Append_Bind_Key('-',coarse_mac_velocity_component->Decrease_Vector_Size_CB());
        opengl_world.Append_Bind_Key('h',coarse_mac_velocity_component->Toggle_Arrowhead_CB());
        opengl_world.Append_Bind_Key('S',coarse_mac_velocity_component->Toggle_Draw_Streamlines_CB());
        opengl_world.Append_Bind_Key('/',coarse_mac_velocity_component->Toggle_Use_Streamline_Seed_CB());
        opengl_world.Append_Bind_Key(']',coarse_mac_velocity_component->Lengthen_Streamlines_CB());
        opengl_world.Append_Bind_Key('[',coarse_mac_velocity_component->Shorten_Streamlines_CB());}

#ifndef COMPILE_WITHOUT_RLE_SUPPORT
    // rle velocities
    filename=basedir+"/%d/rle_face_velocities";
    if(has_valid_rle_grid && FILE_UTILITIES::Frame_File_Exists(filename,start_frame)){
         OPENGL_COMPONENT_RLE_FACE_SCALAR_FIELD_2D<T,T>* rle_face_velocity_component=new OPENGL_COMPONENT_RLE_FACE_SCALAR_FIELD_2D<T,T>(rle_grid_component->opengl_grid.grid,filename,
            new OPENGL_CONSTANT_COLOR_MAP<T>(OPENGL_COLOR::Green()),false);
        rle_face_velocity_component->opengl_rle_face_scalar_field.line_size=.01;
        Add_Component(rle_face_velocity_component,"RLE Face Velocities",'V',BASIC_VISUALIZATION::OWNED|BASIC_VISUALIZATION::START_HIDDEN);
        opengl_world.Append_Bind_Key('C',rle_face_velocity_component->Toggle_Draw_CB());
        opengl_world.Append_Bind_Key('=',rle_face_velocity_component->Increase_Vector_Size_CB());
        opengl_world.Append_Bind_Key('-',rle_face_velocity_component->Decrease_Vector_Size_CB());}
#endif
#ifndef COMPILE_WITHOUT_DYADIC_SUPPORT
    // quadtree velocities
    filename=basedir+"/%d/quadtree_face_velocities";
    if(has_valid_quadtree_grid && FILE_UTILITIES::Frame_File_Exists(filename,start_frame)){
        OPENGL_COMPONENT_QUADTREE_FACE_SCALAR_FIELD<T,T>* quadtree_face_velocity_component=new OPENGL_COMPONENT_QUADTREE_FACE_SCALAR_FIELD<T,T>(quadtree_grid_component->opengl_grid.grid,filename,
            new OPENGL_CONSTANT_COLOR_MAP<T>(OPENGL_COLOR::Magenta()),false);
        quadtree_face_velocity_component->opengl_quadtree_face_scalar_field.line_size=.01;
        Add_Component(quadtree_face_velocity_component,"Quadtree Face Velocities",'V',BASIC_VISUALIZATION::OWNED|BASIC_VISUALIZATION::START_HIDDEN);
        opengl_world.Append_Bind_Key('=',quadtree_face_velocity_component->Increase_Vector_Size_CB());
        opengl_world.Append_Bind_Key('-',quadtree_face_velocity_component->Decrease_Vector_Size_CB());}
    filename=basedir+"/%d/quadtree_face_velocities_fuel";
    if(has_valid_quadtree_grid && FILE_UTILITIES::Frame_File_Exists(filename,start_frame)){
        OPENGL_COMPONENT_QUADTREE_FACE_SCALAR_FIELD<T,T>* quadtree_face_velocity_fuel_component=new OPENGL_COMPONENT_QUADTREE_FACE_SCALAR_FIELD<T,T>(quadtree_grid_component->opengl_grid.grid,filename,
            new OPENGL_CONSTANT_COLOR_MAP<T>(OPENGL_COLOR::Magenta()),false);
        quadtree_face_velocity_fuel_component->opengl_quadtree_face_scalar_field.line_size=.01;
        Add_Component(quadtree_face_velocity_fuel_component,"Quadtree Face Velocities Fuel",'B',BASIC_VISUALIZATION::OWNED|BASIC_VISUALIZATION::START_HIDDEN);
        opengl_world.Append_Bind_Key('=',quadtree_face_velocity_fuel_component->Increase_Vector_Size_CB());
        opengl_world.Append_Bind_Key('-',quadtree_face_velocity_fuel_component->Decrease_Vector_Size_CB());}
#endif

    filename=basedir+"/%d/beta_face";
    if(has_valid_grid && FILE_UTILITIES::Frame_File_Exists(filename,start_frame)){
        OPENGL_COMPONENT_FACE_SCALAR_FIELD_2D<T,T>* beta_face_component=new OPENGL_COMPONENT_FACE_SCALAR_FIELD_2D<T,T>(mac_grid,filename,OPENGL_COLOR_RAMP<T>::Two_Color_Ramp(0,.002,OPENGL_COLOR::Gray(1),OPENGL_COLOR::Gray(0)));
        Add_Component(beta_face_component,"Beta Face",'\0',BASIC_VISUALIZATION::OWNED|BASIC_VISUALIZATION::START_HIDDEN);
        opengl_world.Append_Bind_Key(OPENGL_KEY(OPENGL_KEY::F6),beta_face_component->Toggle_Draw_CB());}

#ifndef COMPILE_WITHOUT_DYADIC_SUPPORT
    filename=basedir+"/%d/quadtree_beta_face";
    if(has_valid_quadtree_grid && FILE_UTILITIES::Frame_File_Exists(filename,start_frame)){
        OPENGL_COMPONENT_QUADTREE_FACE_SCALAR_FIELD<T,T>* quadtree_beta_face_component=new OPENGL_COMPONENT_QUADTREE_FACE_SCALAR_FIELD<T,T>(quadtree_grid_component->opengl_grid.grid,filename,
            new OPENGL_CONSTANT_COLOR_MAP<T>(OPENGL_COLOR::Magenta()),false);
        quadtree_beta_face_component->opengl_quadtree_face_scalar_field.line_size=.01;
        Add_Component(quadtree_beta_face_component,"Quadtree Beta Face",'\0',BASIC_VISUALIZATION::OWNED);}
#endif

    filename=basedir+"/%d/mac_velocities_fuel";
    if(has_valid_grid && FILE_UTILITIES::Frame_File_Exists(filename,start_frame)){
        OPENGL_COMPONENT_MAC_VELOCITY_FIELD_2D<T>* mac_velocity_fuel_component=new OPENGL_COMPONENT_MAC_VELOCITY_FIELD_2D<T>(mac_grid,filename);
        mac_velocity_fuel_component->opengl_mac_velocity_field->size=.01;
        mac_velocity_fuel_component->opengl_mac_velocity_field->vector_color=OPENGL_COLOR::Cyan();
        mac_velocity_fuel_component->opengl_mac_velocity_field->Set_Velocity_Mode(OPENGL_MAC_VELOCITY_FIELD_2D<T>::FACE_CENTERED);
        Add_Component(mac_velocity_fuel_component,"MAC velocities fuel",'B',BASIC_VISUALIZATION::OWNED|BASIC_VISUALIZATION::START_HIDDEN);
        opengl_world.Append_Bind_Key('=',mac_velocity_fuel_component->Increase_Vector_Size_CB());
        opengl_world.Append_Bind_Key('-',mac_velocity_fuel_component->Decrease_Vector_Size_CB());
        opengl_world.Append_Bind_Key('h',mac_velocity_fuel_component->Toggle_Arrowhead_CB());}

    filename=basedir+"/%d/velocities_ghost_fuel";
    OPENGL_COMPONENT_GRID_BASED_VECTOR_FIELD_2D<T>* vector_velocity_ghost_minus_component=0;
    if(has_valid_grid && FILE_UTILITIES::Frame_File_Exists(filename,start_frame)){
        if(node_based){
            vector_velocity_ghost_minus_component=new OPENGL_COMPONENT_GRID_BASED_VECTOR_FIELD_2D<T>(regular_grid,filename);
            vector_velocity_ghost_minus_component->opengl_grid_based_vector_field.size=.01;
            vector_velocity_ghost_minus_component->opengl_grid_based_vector_field.vector_color=OPENGL_COLOR::Green();
            Add_Component(vector_velocity_ghost_minus_component,"Node ghost fuel velocities",'\0',BASIC_VISUALIZATION::OWNED|BASIC_VISUALIZATION::START_HIDDEN);
            opengl_world.Append_Bind_Key('c',vector_velocity_ghost_minus_component->Toggle_Draw_CB());
            opengl_world.Append_Bind_Key('h',vector_velocity_ghost_minus_component->Toggle_Arrowhead_CB());
            opengl_world.Append_Bind_Key('=',vector_velocity_ghost_minus_component->Increase_Vector_Size_CB());
            opengl_world.Append_Bind_Key('-',vector_velocity_ghost_minus_component->Decrease_Vector_Size_CB());}
        else{
            OPENGL_COMPONENT_MAC_VELOCITY_FIELD_2D<T>* mac_velocity_ghost_minus_component=new OPENGL_COMPONENT_MAC_VELOCITY_FIELD_2D<T>(mac_grid,filename);
            mac_velocity_ghost_minus_component->opengl_mac_velocity_field->size=.01;
            mac_velocity_ghost_minus_component->opengl_mac_velocity_field->vector_color=OPENGL_COLOR::Green();
            mac_velocity_ghost_minus_component->Set_Draw(false);
            Add_Component(mac_velocity_ghost_minus_component,"MAC ghost fuel velocities",'\0',BASIC_VISUALIZATION::OWNED|BASIC_VISUALIZATION::START_HIDDEN);
            opengl_world.Append_Bind_Key('c',mac_velocity_ghost_minus_component->Toggle_Velocity_Mode_And_Draw_CB());
            opengl_world.Append_Bind_Key('=',mac_velocity_ghost_minus_component->Increase_Vector_Size_CB());
            opengl_world.Append_Bind_Key('-',mac_velocity_ghost_minus_component->Decrease_Vector_Size_CB());
            opengl_world.Append_Bind_Key('h',mac_velocity_ghost_minus_component->Toggle_Arrowhead_CB());}}



#ifndef COMPILE_WITHOUT_DYADIC_SUPPORT
    filename=basedir+"/%d/quadtree_velocities_ghost_fuel";
    if(has_valid_quadtree_grid && FILE_UTILITIES::Frame_File_Exists(filename,start_frame)){
        OPENGL_COMPONENT_QUADTREE_NODE_VECTOR_FIELD<T>* quadtree_velocity_minus_component=new OPENGL_COMPONENT_QUADTREE_NODE_VECTOR_FIELD<T>(quadtree_grid_component->opengl_grid.grid,filename);
        quadtree_velocity_minus_component->opengl_quadtree_node_vector_field.size=.01;
        quadtree_velocity_minus_component->opengl_quadtree_node_vector_field.vector_color=OPENGL_COLOR::Green();
        Add_Component(quadtree_velocity_minus_component,"Quadtree ghost fuel velocities",'c',BASIC_VISUALIZATION::OWNED|BASIC_VISUALIZATION::START_HIDDEN);
        opengl_world.Append_Bind_Key('h',quadtree_velocity_minus_component->Toggle_Arrowhead_CB());
        opengl_world.Append_Bind_Key('=',quadtree_velocity_minus_component->Increase_Vector_Size_CB());
        opengl_world.Append_Bind_Key('-',quadtree_velocity_minus_component->Decrease_Vector_Size_CB());}
#endif

    filename=basedir+"/%d/velocities_ghost";
    OPENGL_COMPONENT_GRID_BASED_VECTOR_FIELD_2D<T>* vector_velocity_ghost_plus_component=0;
    if(has_valid_grid && FILE_UTILITIES::Frame_File_Exists(filename,start_frame)){
        if(node_based){
            vector_velocity_ghost_plus_component=new OPENGL_COMPONENT_GRID_BASED_VECTOR_FIELD_2D<T>(regular_grid,filename);
            vector_velocity_ghost_plus_component->opengl_grid_based_vector_field.size=.01;
            vector_velocity_ghost_plus_component->opengl_grid_based_vector_field.vector_color=OPENGL_COLOR::Green();
            vector_velocity_ghost_plus_component->Set_Draw(false);
            Add_Component(vector_velocity_ghost_plus_component,"Node ghost velocities",'b',BASIC_VISUALIZATION::OWNED|BASIC_VISUALIZATION::START_HIDDEN);
            opengl_world.Append_Bind_Key('h',vector_velocity_ghost_plus_component->Toggle_Arrowhead_CB());
            opengl_world.Append_Bind_Key('=',vector_velocity_ghost_plus_component->Increase_Vector_Size_CB());
            opengl_world.Append_Bind_Key('-',vector_velocity_ghost_plus_component->Decrease_Vector_Size_CB());}
        else{
            OPENGL_COMPONENT_MAC_VELOCITY_FIELD_2D<T>* mac_velocity_ghost_plus_component=new OPENGL_COMPONENT_MAC_VELOCITY_FIELD_2D<T>(mac_grid,filename);
            mac_velocity_ghost_plus_component->opengl_mac_velocity_field->size=.01;
            mac_velocity_ghost_plus_component->opengl_mac_velocity_field->vector_color=OPENGL_COLOR::Green();
            mac_velocity_ghost_plus_component->Set_Draw(false);
            Add_Component(mac_velocity_ghost_plus_component,"MAC ghost velocities",'\0',BASIC_VISUALIZATION::OWNED|BASIC_VISUALIZATION::START_HIDDEN);
            opengl_world.Append_Bind_Key('b',mac_velocity_ghost_plus_component->Toggle_Velocity_Mode_And_Draw_CB());
            opengl_world.Append_Bind_Key('=',mac_velocity_ghost_plus_component->Increase_Vector_Size_CB());
            opengl_world.Append_Bind_Key('-',mac_velocity_ghost_plus_component->Decrease_Vector_Size_CB());
            opengl_world.Append_Bind_Key('h',mac_velocity_ghost_plus_component->Toggle_Arrowhead_CB());}}

#ifndef COMPILE_WITHOUT_DYADIC_SUPPORT
    filename=basedir+"/%d/quadtree_velocities_ghost";
    if(has_valid_quadtree_grid && FILE_UTILITIES::Frame_File_Exists(filename,start_frame)){
        OPENGL_COMPONENT_QUADTREE_NODE_VECTOR_FIELD<T>* quadtree_velocity_plus_component=new OPENGL_COMPONENT_QUADTREE_NODE_VECTOR_FIELD<T>(quadtree_grid_component->opengl_grid.grid,filename);
        quadtree_velocity_plus_component->opengl_quadtree_node_vector_field.size=.01;
        quadtree_velocity_plus_component->opengl_quadtree_node_vector_field.vector_color=OPENGL_COLOR::Green();
        Add_Component(quadtree_velocity_plus_component,"Quadtree ghost velocities",'b',BASIC_VISUALIZATION::START_HIDDEN|BASIC_VISUALIZATION::OWNED);
        opengl_world.Append_Bind_Key('h',quadtree_velocity_plus_component->Toggle_Arrowhead_CB());
        opengl_world.Append_Bind_Key('=',quadtree_velocity_plus_component->Increase_Vector_Size_CB());
        opengl_world.Append_Bind_Key('-',quadtree_velocity_plus_component->Decrease_Vector_Size_CB());}
#endif

    if(has_valid_grid && vector_velocity_ghost_minus_component && vector_velocity_ghost_plus_component && levelset_component){
        OPENGL_COMPONENT_TWO_PHASE_VELOCITY_MAGNITUDE_2D<T>* two_phase_velocity_magnitude_component=new OPENGL_COMPONENT_TWO_PHASE_VELOCITY_MAGNITUDE_2D<T>
            (*vector_velocity_ghost_minus_component,*vector_velocity_ghost_plus_component,*levelset_component);
        two_phase_velocity_magnitude_component->opengl_two_phase_velocity_magnitude.plus.size=.01;
        two_phase_velocity_magnitude_component->opengl_two_phase_velocity_magnitude.minus.size=.01;
        two_phase_velocity_magnitude_component->opengl_two_phase_velocity_magnitude.plus.vector_color=OPENGL_COLOR::Magenta();
        two_phase_velocity_magnitude_component->opengl_two_phase_velocity_magnitude.minus.vector_color=OPENGL_COLOR::Cyan();
        Add_Component(two_phase_velocity_magnitude_component,"Two Phase Magnitude",'2',BASIC_VISUALIZATION::OWNED|BASIC_VISUALIZATION::START_HIDDEN);
        opengl_world.Append_Bind_Key('\\',two_phase_velocity_magnitude_component->Toggle_3D_Mode_CB());
        //opengl_world.Append_Bind_Key('h',two_phase_velocity_magnitude_component->Toggle_Arrowhead_CB());
        opengl_world.Append_Bind_Key('=',two_phase_velocity_magnitude_component->Increase_Point_Size_CB());
        opengl_world.Append_Bind_Key('-',two_phase_velocity_magnitude_component->Decrease_Point_Size_CB());}

    filename=basedir+"/%d/center_velocities";
    if(has_valid_grid && FILE_UTILITIES::Frame_File_Exists(filename,start_frame)){
        OPENGL_COMPONENT_GRID_BASED_VECTOR_FIELD_2D<T>* center_velocity_component=new OPENGL_COMPONENT_GRID_BASED_VECTOR_FIELD_2D<T>(mac_grid,filename);
        center_velocity_component->opengl_grid_based_vector_field.size=.01;
        center_velocity_component->opengl_grid_based_vector_field.vector_color=OPENGL_COLOR::Green();
        Add_Component(center_velocity_component,"Centered velocities",'V',BASIC_VISUALIZATION::OWNED|BASIC_VISUALIZATION::START_HIDDEN);
        opengl_world.Append_Bind_Key('=',center_velocity_component->Increase_Vector_Size_CB());
        opengl_world.Append_Bind_Key('-',center_velocity_component->Decrease_Vector_Size_CB());
        opengl_world.Append_Bind_Key('h',center_velocity_component->Toggle_Arrowhead_CB());}

#ifndef COMPILE_WITHOUT_DYADIC_SUPPORT
    filename=basedir+"/%d/quadtree_velocities";
    if(has_valid_quadtree_grid && FILE_UTILITIES::Frame_File_Exists(filename,start_frame)){
        OPENGL_COMPONENT_QUADTREE_NODE_VECTOR_FIELD<T>* quadtree_velocity_component=new OPENGL_COMPONENT_QUADTREE_NODE_VECTOR_FIELD<T>(quadtree_grid_component->opengl_grid.grid,filename);
        quadtree_velocity_component->opengl_quadtree_node_vector_field.size=.01;
        quadtree_velocity_component->opengl_quadtree_node_vector_field.vector_color=OPENGL_COLOR::Green();
        Add_Component(quadtree_velocity_component,"Quadtree Velocities",'v',BASIC_VISUALIZATION::OWNED|BASIC_VISUALIZATION::START_HIDDEN);
        opengl_world.Append_Bind_Key('=',quadtree_velocity_component->Increase_Vector_Size_CB());
        opengl_world.Append_Bind_Key('-',quadtree_velocity_component->Decrease_Vector_Size_CB());
        opengl_world.Append_Bind_Key('h',quadtree_velocity_component->Toggle_Arrowhead_CB());}
#endif

    opengl_world.Set_Key_Binding_Category("Pressure");
    // TODO: these ramps are leaking memory
    OPENGL_COLOR_MAP<T>* pressure_for_pressure_coupling_color_map=OPENGL_COLOR_RAMP<T>::Two_Color_Ramp(-10000,10000,OPENGL_COLOR::Cyan(0,0),OPENGL_COLOR::Cyan(1));
    filename=basedir+"/%d/pressure";
    if(has_valid_grid && FILE_UTILITIES::Frame_File_Exists(filename,start_frame)){
        pressure_component=new OPENGL_COMPONENT_SCALAR_FIELD_2D<T>(mac_grid,filename,OPENGL_COLOR_RAMP<T>::Two_Color_Ramp(0,1,OPENGL_COLOR::Cyan(0,0),OPENGL_COLOR::Cyan(1)));
        pressure_component->opengl_scalar_field.Set_Scale_Range(0,30);
        pressure_component->opengl_scalar_field.Set_Uniform_Contour_Values(2,28,(T)26/(T)53);
        Add_Component(pressure_component,"Pressure",'7',BASIC_VISUALIZATION::OWNED|BASIC_VISUALIZATION::START_HIDDEN);
        opengl_world.Append_Bind_Key(OPENGL_KEY(OPENGL_KEY::F10),pressure_component->Toggle_Draw_Mode_CB());}

    opengl_world.Set_Key_Binding_Category("Coarse Pressure");
    // TODO: these ramps are leaking memory
    filename=basedir+"/%d/coarse_pressure";
    if(has_valid_coarse_grid && FILE_UTILITIES::Frame_File_Exists(filename,start_frame)){
        coarse_pressure_component=new OPENGL_COMPONENT_SCALAR_FIELD_2D<T>(coarse_mac_grid,filename,OPENGL_COLOR_RAMP<T>::Two_Color_Ramp(0,1,OPENGL_COLOR::Cyan(0,0),OPENGL_COLOR::Cyan(1)));
        coarse_pressure_component->opengl_scalar_field.Set_Scale_Range(-20,20);
        coarse_pressure_component->opengl_scalar_field.Set_Uniform_Contour_Values(0,20,1);
        Add_Component(coarse_pressure_component,"Coarse Pressure",'7',BASIC_VISUALIZATION::OWNED|BASIC_VISUALIZATION::START_HIDDEN);
        opengl_world.Append_Bind_Key(OPENGL_KEY(OPENGL_KEY::F10),coarse_pressure_component->Toggle_Draw_Mode_CB());}

    opengl_world.Set_Key_Binding_Category("Compressible_Implicit_Pressure");
    filename=basedir+"/%d/compressible_implicit_pressure";
    if(has_valid_grid && FILE_UTILITIES::Frame_File_Exists(filename,start_frame)){
        OPENGL_COMPONENT_SCALAR_FIELD_2D<T>* compressible_implicit_pressure_component=new OPENGL_COMPONENT_SCALAR_FIELD_2D<T>(mac_grid,filename,OPENGL_COLOR_RAMP<T>::Matlab_Jet(0,1));
        compressible_implicit_pressure_component->opengl_scalar_field.Set_Uniform_Contour_Values(0,20,1);
        Add_Component(compressible_implicit_pressure_component,"Compressible_Implicit_Pressure",'\0',BASIC_VISUALIZATION::OWNED|BASIC_VISUALIZATION::START_HIDDEN);
        opengl_world.Append_Bind_Key(OPENGL_KEY(OPENGL_KEY::F9),compressible_implicit_pressure_component->Toggle_Draw_CB());}
    filename=basedir+"/%d/pressure_gradient";
    if(has_valid_grid && FILE_UTILITIES::Frame_File_Exists(filename,start_frame)){
        OPENGL_COMPONENT_SCALAR_FIELD_2D<T>* pressure_gradient_component=new OPENGL_COMPONENT_SCALAR_FIELD_2D<T>(grid,filename,OPENGL_COLOR_RAMP<T>::Two_Color_Ramp(0,1,OPENGL_COLOR::Gray(0,0),OPENGL_COLOR::Magenta(1,1)));
        pressure_gradient_component->opengl_scalar_field.Update();
        Add_Component(pressure_gradient_component,"Pressure Gradient",'9',BASIC_VISUALIZATION::OWNED|BASIC_VISUALIZATION::START_HIDDEN);}
    filename=basedir+"/%d/pressure_for_pressure_coupling";
    if(has_valid_grid && FILE_UTILITIES::Frame_File_Exists(filename,start_frame)){
        OPENGL_COMPONENT_SCALAR_FIELD_2D<T>* pressure2_component=new OPENGL_COMPONENT_SCALAR_FIELD_2D<T>(mac_grid,filename,pressure_for_pressure_coupling_color_map);
        pressure2_component->opengl_scalar_field.Set_Uniform_Contour_Values(-10000,10000,100);
        Add_Component(pressure2_component,"Pressure2",'\0',BASIC_VISUALIZATION::OWNED|BASIC_VISUALIZATION::START_HIDDEN);
        opengl_world.Append_Bind_Key(OPENGL_KEY(OPENGL_KEY::F7),pressure2_component->Toggle_Draw_CB());}
#ifndef COMPILE_WITHOUT_RLE_SUPPORT
    filename=basedir+"/%d/rle_pressure";
    if(has_valid_rle_grid && FILE_UTILITIES::Frame_File_Exists(filename,start_frame)){
        OPENGL_COMPONENT_RLE_CELL_SCALAR_FIELD_2D<T>* rle_pressure_component=new OPENGL_COMPONENT_RLE_CELL_SCALAR_FIELD_2D<T>(rle_grid_component->opengl_grid.grid,filename,
            OPENGL_COLOR_RAMP<T>::Matlab_Jet(0,1));
        Add_Component(rle_pressure_component,"RLE Pressure",'7',BASIC_VISUALIZATION::OWNED|BASIC_VISUALIZATION::START_HIDDEN);}
#endif
#ifndef COMPILE_WITHOUT_DYADIC_SUPPORT
    filename=basedir+"/%d/quadtree_pressure";
    if(has_valid_quadtree_grid && FILE_UTILITIES::Frame_File_Exists(filename,start_frame)){
        OPENGL_COMPONENT_QUADTREE_CELL_SCALAR_FIELD<T>* quadtree_pressure_component=new OPENGL_COMPONENT_QUADTREE_CELL_SCALAR_FIELD<T>(quadtree_grid_component->opengl_grid.grid,filename,
            OPENGL_COLOR_RAMP<T>::Matlab_Jet(0,1));
        Add_Component(quadtree_pressure_component,"Quadtree Pressure",'7',BASIC_VISUALIZATION::OWNED|BASIC_VISUALIZATION::START_HIDDEN);}
#endif

    opengl_world.Set_Key_Binding_Category("Volume of Material");
    filename=basedir+"/%d/negative_material";
    if(has_valid_grid && FILE_UTILITIES::Frame_File_Exists(filename,start_frame)){
        OPENGL_COMPONENT_SCALAR_FIELD_2D<T>* volume_of_material_component=new OPENGL_COMPONENT_SCALAR_FIELD_2D<T>(mac_grid,filename,OPENGL_COLOR_RAMP<T>::Matlab_Jet(0,1));
        volume_of_material_component->opengl_scalar_field.Set_Uniform_Contour_Values(-2,2,(T).1);
        Add_Component(volume_of_material_component,"Volume of Material",'0',BASIC_VISUALIZATION::OWNED|BASIC_VISUALIZATION::START_HIDDEN);}
    filename=basedir+"/%d/vof_object";
    if(has_valid_grid && FILE_UTILITIES::Frame_File_Exists(filename,start_frame)){
        std::string color_file=basedir+"/%d/vof_colors";
        Add_Component(new OPENGL_COMPONENT_TRIANGULATED_AREA<T>(filename,*new std::string(color_file)),
            "VOF Refinement",'Z',BASIC_VISUALIZATION::OWNED|BASIC_VISUALIZATION::START_HIDDEN);}

    // Draw grid here so it'll be above particles and pressure
    opengl_world.Set_Key_Binding_Category("Grid");
    if(has_valid_grid){
        OPENGL_GRID_2D<T>* opengl_grid=new OPENGL_GRID_2D<T>(*(new GRID<TV>(grid)),OPENGL_COLOR::Gray(.5));
        grid_component=new OPENGL_COMPONENT_BASIC<OPENGL_GRID_2D<T> >(*opengl_grid);
        Add_Component(grid_component,"Grid",'6',BASIC_VISUALIZATION::SELECTABLE);
        opengl_world.Append_Bind_Key('^',grid_component->object.Toggle_Draw_Ghost_Values_CB());}
#ifndef COMPILE_WITHOUT_RLE_SUPPORT
    if(rle_grid_component)
        Add_Component(rle_grid_component,"RLE Grid",'6',BASIC_VISUALIZATION::SELECTABLE);
#endif
#ifndef COMPILE_WITHOUT_DYADIC_SUPPORT
    if(quadtree_grid_component) 
        Add_Component(quadtree_grid_component,"Quadtree Grid",'6',BASIC_VISUALIZATION::SELECTABLE);
#endif
    opengl_world.Set_Key_Binding_Category("Coarse Grid");
    if(has_valid_coarse_grid){
        OPENGL_GRID_2D<T>* opengl_grid=new OPENGL_GRID_2D<T>(*(new GRID<TV>(coarse_grid)),OPENGL_COLOR::Ground_Tan(.5));
        coarse_grid_component=new OPENGL_COMPONENT_BASIC<OPENGL_GRID_2D<T> >(*opengl_grid);
        Add_Component(coarse_grid_component,"Coarse Grid",'y',BASIC_VISUALIZATION::SELECTABLE);
        opengl_world.Append_Bind_Key('^',coarse_grid_component->object.Toggle_Draw_Ghost_Values_CB());}
    opengl_world.Set_Key_Binding_Category("Sub Grids");

    // deformable and rigid bodies
    OPENGL_COMPONENT_RIGID_GEOMETRY_COLLECTION_2D<T>* rigid_bodies_component=0;
    //OPENGL_COMPONENT_RIGID_BODIES_2D<T>* rigid_bodies_component=0;
    OPENGL_COMPONENT_DEFORMABLE_GEOMETRY_COLLECTION_2D<T>* deformable_objects_component=0;
    std::string deformable_object_filename=STRING_UTILITIES::string_sprintf("%s/%d/deformable_object_particles",basedir.c_str(),start_frame);
    if(FILE_UTILITIES::File_Exists(deformable_object_filename)){
        deformable_objects_component=new OPENGL_COMPONENT_DEFORMABLE_GEOMETRY_COLLECTION_2D<T>(basedir+"/",start_frame);
        deformable_objects_component->Set_Vector_Size(.01);
        deformable_objects_component->selectable=true;
        // TODO: what the hell?
        if(!FILE_UTILITIES::Frame_File_Exists(basedir+"/%d/rigid_geometry_particles",start_frame)){
            }}
    if(FILE_UTILITIES::Frame_File_Exists(basedir+"/%d/rigid_geometry_particles",start_frame)){
        //if(deformable_objects_component) rigid_bodies_component=new OPENGL_COMPONENT_RIGID_BODIES_2D<T>(deformable_objects_component->solid_body_collection.rigid_body_collection,basedir);
        //else rigid_bodies_component=new OPENGL_COMPONENT_RIGID_BODIES_2D<T>(basedir);
        rigid_bodies_component=new OPENGL_COMPONENT_RIGID_GEOMETRY_COLLECTION_2D<T>(0,basedir);
        //rigid_bodies_component=new OPENGL_COMPONENT_RIGID_BODIES_2D<T>(basedir);
        rigid_bodies_component->Set_Vector_Size(.01);
        rigid_bodies_component->selectable=true;
        if(FILE_UTILITIES::Frame_File_Exists(basedir+"/%d/colors",start_frame)) FILE_UTILITIES::template Read_From_File<RW>(STRING_UTILITIES::string_sprintf("%s/%d/colors",basedir.c_str(),start_frame),rigid_bodies_component->colors);
        for(int i=1;i<=rigid_bodies_no_draw_list.m;i++){
            LOG::cout<<"Rigid bodies: not drawing object "<<rigid_bodies_no_draw_list(i)<<std::endl;
            rigid_bodies_component->Set_Draw_Object(rigid_bodies_no_draw_list(i),false);}
        opengl_world.Set_Key_Binding_Category("Rigid Bodies");
        Add_Component(rigid_bodies_component,"Rigid Bodies",'\0',BASIC_VISUALIZATION::OWNED|BASIC_VISUALIZATION::SELECTABLE);
        opengl_world.Append_Bind_Key('5',rigid_bodies_component->Toggle_Draw_Mode_CB());
        opengl_world.Append_Bind_Key('%',rigid_bodies_component->Toggle_Velocity_Vectors_CB());
        opengl_world.Append_Bind_Key('%',rigid_bodies_component->Toggle_Show_Object_Names_CB());
        opengl_world.Append_Bind_Key('=',rigid_bodies_component->Increase_Vector_Size_CB());
        opengl_world.Append_Bind_Key('-',rigid_bodies_component->Decrease_Vector_Size_CB());}
    if(deformable_objects_component){
        opengl_world.Set_Key_Binding_Category("Deformable Objects");
        Add_Component(deformable_objects_component,"Deformable Objects",'\0',BASIC_VISUALIZATION::OWNED|BASIC_VISUALIZATION::SELECTABLE);
        opengl_world.Append_Bind_Key('8',deformable_objects_component->Toggle_Draw_CB());
        opengl_world.Append_Bind_Key('e',deformable_objects_component->Cycle_Display_Mode_CB());
        opengl_world.Append_Bind_Key('*',deformable_objects_component->Toggle_Draw_Velocities_CB());
        opengl_world.Append_Bind_Key('=',deformable_objects_component->Increase_Vector_Size_CB());
        opengl_world.Append_Bind_Key('-',deformable_objects_component->Decrease_Vector_Size_CB());}

    std::string soft_constraints_deformable_object_filename=basedir+STRING_UTILITIES::string_sprintf("/%d/soft_constraints_deformable_object_particles",start_frame);
    if(FILE_UTILITIES::File_Exists(soft_constraints_deformable_object_filename)){
        OPENGL_COMPONENT_DEFORMABLE_GEOMETRY_COLLECTION_2D<T>* soft_constraints_deformable_objects_component=new OPENGL_COMPONENT_DEFORMABLE_GEOMETRY_COLLECTION_2D<T>(basedir+"/soft_constraints_",start_frame);
        soft_constraints_deformable_objects_component->Set_Vector_Size(.01);
        Add_Component(soft_constraints_deformable_objects_component,"Soft Constraints Deformable Objects",'e',BASIC_VISUALIZATION::OWNED|BASIC_VISUALIZATION::SELECTABLE);}

    opengl_world.Set_Key_Binding_Category("Fluid Boundaries");
    if(has_valid_coarse_grid && FILE_UTILITIES::Frame_File_Exists(basedir+"/%d/coarse_psi_N",start_frame)){
        OPENGL_COMPONENT_FACE_SCALAR_FIELD_2D<T,bool>* psi_N_component=new OPENGL_COMPONENT_FACE_SCALAR_FIELD_2D<T,bool>(coarse_grid,basedir+"/%d/coarse_psi_N",
            new OPENGL_CONSTANT_COLOR_MAP<bool>(OPENGL_COLOR::Cyan()));
        Add_Component(psi_N_component,"Coarse Psi_N points",'\0',BASIC_VISUALIZATION::START_HIDDEN|BASIC_VISUALIZATION::OWNED);
        opengl_world.Append_Bind_Key(OPENGL_KEY(OPENGL_KEY::F2),psi_N_component->Toggle_Draw_CB());}
    if(has_valid_grid && FILE_UTILITIES::Frame_File_Exists(basedir+"/%d/psi_N",start_frame)){
        OPENGL_COMPONENT_FACE_SCALAR_FIELD_2D<T,bool>* psi_N_component=new OPENGL_COMPONENT_FACE_SCALAR_FIELD_2D<T,bool>(grid,basedir+"/%d/psi_N",
            new OPENGL_CONSTANT_COLOR_MAP<bool>(OPENGL_COLOR::Cyan()));
        Add_Component(psi_N_component,"Psi_N points",'\0',BASIC_VISUALIZATION::START_HIDDEN|BASIC_VISUALIZATION::OWNED);
        opengl_world.Append_Bind_Key(OPENGL_KEY(OPENGL_KEY::F1),psi_N_component->Toggle_Draw_CB());}
    filename=basedir+"/%d/debug_particles";
    if(FILE_UTILITIES::Frame_File_Exists(filename,start_frame)){
        OPENGL_COMPONENT_DEBUG_PARTICLES_2D<T>* component=new OPENGL_COMPONENT_DEBUG_PARTICLES_2D<T>(filename);
        Add_Component(component,"Debug particles",'w',BASIC_VISUALIZATION::START_HIDDEN|BASIC_VISUALIZATION::SELECTABLE|BASIC_VISUALIZATION::OWNED);
        opengl_world.Append_Bind_Key('W',component->Toggle_Draw_Velocities_CB());
        opengl_world.Append_Bind_Key('h',component->Toggle_Arrowhead_CB());
        opengl_world.Append_Bind_Key('=',component->Increase_Vector_Size_CB());
        opengl_world.Append_Bind_Key('-',component->Decrease_Vector_Size_CB());}
    filename=basedir+"/%d/residual_energy";
    if(FILE_UTILITIES::Frame_File_Exists(filename,start_frame)){
        OPENGL_COMPONENT_PARTICLES_2D<T,GEOMETRY_PARTICLES<TV> >* component=new OPENGL_COMPONENT_PARTICLES_2D<T,GEOMETRY_PARTICLES<TV> >(filename,"",false,false);
        Add_Component(component,"Residual Energy",'k',BASIC_VISUALIZATION::START_HIDDEN|BASIC_VISUALIZATION::OWNED);}
    filename=basedir+"/%d/collision_iterators";
    if(has_valid_grid && FILE_UTILITIES::Frame_File_Exists(filename,start_frame)){
        OPENGL_COMPONENT_PARTICLES_2D<T,GEOMETRY_PARTICLES<TV> >* component=new OPENGL_COMPONENT_PARTICLES_2D<T,GEOMETRY_PARTICLES<TV> >(filename,"",false,false);
        Add_Component(component,"Collision Iterators",'I',BASIC_VISUALIZATION::START_HIDDEN|BASIC_VISUALIZATION::OWNED);}
#ifndef COMPILE_WITHOUT_RLE_SUPPORT
    filename=basedir+"/%d/rle_psi_N";
    if(has_valid_rle_grid && FILE_UTILITIES::Frame_File_Exists(filename,start_frame)){
        OPENGL_COMPONENT_RLE_FACE_SCALAR_FIELD_2D<T,bool>* rle_psi_N_component=new OPENGL_COMPONENT_RLE_FACE_SCALAR_FIELD_2D<T,bool>(rle_grid_component->opengl_grid.grid,filename,new OPENGL_CONSTANT_COLOR_MAP<bool>(OPENGL_COLOR::Cyan()),true);
        rle_psi_N_component->Set_Draw(false);
        Add_Component(rle_psi_N_component,"RLE Psi_N points",'\0',BASIC_VISUALIZATION::START_HIDDEN|BASIC_VISUALIZATION::OWNED);
        opengl_world.Append_Bind_Key(OPENGL_KEY(OPENGL_KEY::F1),rle_psi_N_component->Toggle_Draw_CB());}
#endif
#ifndef COMPILE_WITHOUT_DYADIC_SUPPORT
    filename=basedir+"/%d/quadtree_psi_N";
    if(has_valid_quadtree_grid && FILE_UTILITIES::Frame_File_Exists(filename,start_frame)){
        OPENGL_COMPONENT_QUADTREE_FACE_SCALAR_FIELD<T,bool>* quadtree_psi_N_component=new OPENGL_COMPONENT_QUADTREE_FACE_SCALAR_FIELD<T,bool>(quadtree_grid_component->opengl_grid.grid,filename,new OPENGL_CONSTANT_COLOR_MAP<bool>(OPENGL_COLOR::Cyan()),true);
        quadtree_psi_N_component->Set_Draw(false);
        Add_Component(quadtree_psi_N_component,"Quadtree Psi_N points",'\0',BASIC_VISUALIZATION::START_HIDDEN|BASIC_VISUALIZATION::OWNED);
        opengl_world.Append_Bind_Key(OPENGL_KEY(OPENGL_KEY::F1),quadtree_psi_N_component->Toggle_Draw_CB());}
#endif
    filename=basedir+"/%d/coarse_psi_D";
    if(has_valid_coarse_grid && FILE_UTILITIES::Frame_File_Exists(filename,start_frame)){
        OPENGL_COMPONENT_SCALAR_FIELD_2D<T,bool>* psi_D_component=new OPENGL_COMPONENT_SCALAR_FIELD_2D<T,bool>(coarse_mac_grid,filename,new OPENGL_CONSTANT_COLOR_MAP<bool>(OPENGL_COLOR::Magenta()),OPENGL_SCALAR_FIELD_2D<T,bool>::DRAW_POINTS);
        Add_Component(psi_D_component,"Coarse Psi_D points",'\0',BASIC_VISUALIZATION::START_HIDDEN|BASIC_VISUALIZATION::OWNED);
        opengl_world.Append_Bind_Key(OPENGL_KEY(OPENGL_KEY::F2),psi_D_component->Toggle_Draw_CB());}
    filename=basedir+"/%d/psi_D";
    if(has_valid_grid && FILE_UTILITIES::Frame_File_Exists(filename,start_frame)){
        OPENGL_COMPONENT_SCALAR_FIELD_2D<T,bool>* psi_D_component=new OPENGL_COMPONENT_SCALAR_FIELD_2D<T,bool>(mac_grid,filename,new OPENGL_CONSTANT_COLOR_MAP<bool>(OPENGL_COLOR::Magenta()),OPENGL_SCALAR_FIELD_2D<T,bool>::DRAW_POINTS);
        Add_Component(psi_D_component,"Psi_D points",'\0',BASIC_VISUALIZATION::START_HIDDEN|BASIC_VISUALIZATION::OWNED);
        opengl_world.Append_Bind_Key(OPENGL_KEY(OPENGL_KEY::F1),psi_D_component->Toggle_Draw_CB());}
#ifndef COMPILE_WITHOUT_RLE_SUPPORT
    filename=basedir+"/%d/rle_psi_D";
    if(has_valid_rle_grid && FILE_UTILITIES::Frame_File_Exists(filename,start_frame)){
        OPENGL_COMPONENT_RLE_CELL_SCALAR_FIELD_2D<T,bool>* rle_psi_D_component=new OPENGL_COMPONENT_RLE_CELL_SCALAR_FIELD_2D<T,bool>(rle_grid_component->opengl_grid.grid,filename,new OPENGL_CONSTANT_COLOR_MAP<bool>(OPENGL_COLOR::Magenta()));
        rle_psi_D_component->Set_Draw(false);
        Add_Component(rle_psi_D_component,"RLE Psi_D points",'\0',BASIC_VISUALIZATION::START_HIDDEN|BASIC_VISUALIZATION::OWNED);
        opengl_world.Append_Bind_Key(OPENGL_KEY(OPENGL_KEY::F1),rle_psi_D_component->Toggle_Draw_CB());}
#endif
#ifndef COMPILE_WITHOUT_DYADIC_SUPPORT
    filename=basedir+"/%d/quadtree_psi_D";
    if(has_valid_quadtree_grid && FILE_UTILITIES::Frame_File_Exists(filename,start_frame)){
        OPENGL_COMPONENT_QUADTREE_CELL_SCALAR_FIELD<T,bool>* quadtree_psi_D_component=new OPENGL_COMPONENT_QUADTREE_CELL_SCALAR_FIELD<T,bool>(quadtree_grid_component->opengl_grid.grid,filename,new OPENGL_CONSTANT_COLOR_MAP<bool>(OPENGL_COLOR::Magenta()));
        quadtree_psi_D_component->Set_Draw(false);
        Add_Component(quadtree_psi_D_component,"Quadtree Psi_D points",'\0',BASIC_VISUALIZATION::START_HIDDEN|BASIC_VISUALIZATION::OWNED);
        opengl_world.Append_Bind_Key(OPENGL_KEY(OPENGL_KEY::F1),quadtree_psi_D_component->Toggle_Draw_CB());}        
#endif

    filename=basedir+"/%d/euler_psi";
    if(has_valid_grid && FILE_UTILITIES::Frame_File_Exists(filename,start_frame)){
        OPENGL_COMPONENT_SCALAR_FIELD_2D<T,bool>* psi_component=new OPENGL_COMPONENT_SCALAR_FIELD_2D<T,bool>(mac_grid,filename,new OPENGL_CONSTANT_COLOR_MAP<bool>(OPENGL_COLOR::Red()),OPENGL_SCALAR_FIELD_2D<T,bool>::DRAW_POINTS);
        Add_Component(psi_component,"Psi points",'\0',BASIC_VISUALIZATION::START_HIDDEN|BASIC_VISUALIZATION::OWNED);
        opengl_world.Append_Bind_Key(OPENGL_KEY(OPENGL_KEY::F3),psi_component->Toggle_Draw_CB());}

    OPENGL_INDEXED_COLOR_MAP* colors_color_map=OPENGL_INDEXED_COLOR_MAP::Basic_16_Color_Map();colors_color_map->Set_Index_Mode(OPENGL_INDEXED_COLOR_MAP::PERIODIC);
    filename=basedir+"/%d/colors";
    if(has_valid_grid && FILE_UTILITIES::Frame_File_Exists(filename,start_frame)){
        OPENGL_COMPONENT_SCALAR_FIELD_2D<T,int>* psi_colors_component=new OPENGL_COMPONENT_SCALAR_FIELD_2D<T,int>(mac_grid,filename,colors_color_map,OPENGL_SCALAR_FIELD_2D<T,int>::DRAW_POINTS);
        Add_Component(psi_colors_component,"Psi colors",'\0',BASIC_VISUALIZATION::OWNED|BASIC_VISUALIZATION::START_HIDDEN);
        opengl_world.Append_Bind_Key(OPENGL_KEY(OPENGL_KEY::F2),psi_colors_component->Toggle_Draw_CB());}

#ifndef COMPILE_WITHOUT_RLE_SUPPORT
    filename=basedir+"/%d/rle_colors";
    if(has_valid_rle_grid && FILE_UTILITIES::Frame_File_Exists(filename,start_frame)){
        OPENGL_COMPONENT_RLE_CELL_SCALAR_FIELD_2D<T,int>* rle_psi_colors_component=new OPENGL_COMPONENT_RLE_CELL_SCALAR_FIELD_2D<T,int>(rle_grid_component->opengl_grid.grid,filename,colors_color_map);
        rle_psi_colors_component->Set_Draw(false);
        Add_Component(rle_psi_colors_component,"RLE Psi colors",'\0',BASIC_VISUALIZATION::OWNED|BASIC_VISUALIZATION::START_HIDDEN);
        opengl_world.Append_Bind_Key(OPENGL_KEY(OPENGL_KEY::F2),rle_psi_colors_component->Toggle_Draw_CB());}
#endif

#ifndef COMPILE_WITHOUT_DYADIC_SUPPORT
    filename=basedir+"/%d/quadtree_colors";
    if(has_valid_quadtree_grid && FILE_UTILITIES::Frame_File_Exists(filename,start_frame)){
        OPENGL_COMPONENT_QUADTREE_CELL_SCALAR_FIELD<T,int>* quadtree_psi_colors_component=new OPENGL_COMPONENT_QUADTREE_CELL_SCALAR_FIELD<T,int>(quadtree_grid_component->opengl_grid.grid,filename,colors_color_map);
        quadtree_psi_colors_component->Set_Draw(false);
        Add_Component(quadtree_psi_colors_component,"Quadtree Psi colors",'\0',BASIC_VISUALIZATION::OWNED|BASIC_VISUALIZATION::START_HIDDEN);
        opengl_world.Append_Bind_Key(OPENGL_KEY(OPENGL_KEY::F2),quadtree_psi_colors_component->Toggle_Draw_CB());}
#endif

    filename=basedir+"/%d/pressure_jumps";
    if(has_valid_grid && FILE_UTILITIES::Frame_File_Exists(filename,start_frame)){
        OPENGL_COMPONENT_GRID_BASED_VECTOR_FIELD_2D<T>* pressure_jump_component=new OPENGL_COMPONENT_GRID_BASED_VECTOR_FIELD_2D<T>(grid,filename);
        pressure_jump_component->opengl_grid_based_vector_field.size=.01;
        pressure_jump_component->opengl_grid_based_vector_field.vector_color=OPENGL_COLOR::Magenta();
        Add_Component(pressure_jump_component,"Pressure jumps",'&',BASIC_VISUALIZATION::OWNED|BASIC_VISUALIZATION::START_HIDDEN);
        opengl_world.Append_Bind_Key('=',pressure_jump_component->Increase_Vector_Size_CB());
        opengl_world.Append_Bind_Key('-',pressure_jump_component->Decrease_Vector_Size_CB());}

#ifndef COMPILE_WITHOUT_DYADIC_SUPPORT
    filename=basedir+"/%d/quadtree_pressure_jumps";
    if(has_valid_quadtree_grid && FILE_UTILITIES::Frame_File_Exists(filename,start_frame)){
        OPENGL_COMPONENT_QUADTREE_NODE_VECTOR_FIELD<T>* quadtree_pressure_jump_component=new OPENGL_COMPONENT_QUADTREE_NODE_VECTOR_FIELD<T>(quadtree_grid_component->opengl_grid.grid,filename);
        quadtree_pressure_jump_component->opengl_quadtree_node_vector_field.size=.01;
        quadtree_pressure_jump_component->opengl_quadtree_node_vector_field.vector_color=OPENGL_COLOR::Magenta();
        Add_Component(quadtree_pressure_jump_component,"Quadtree pressure jumps",'&',BASIC_VISUALIZATION::OWNED|BASIC_VISUALIZATION::START_HIDDEN);
        opengl_world.Append_Bind_Key('=',quadtree_pressure_jump_component->Increase_Vector_Size_CB());
        opengl_world.Append_Bind_Key('-',quadtree_pressure_jump_component->Decrease_Vector_Size_CB());}
#endif

    filename=basedir+"/%d/thin_shells_grid_visibility";
    if(has_valid_grid && FILE_UTILITIES::Frame_File_Exists(filename,start_frame)){
        OPENGL_COMPONENT_THIN_SHELLS_DEBUGGING_2D<T>* thin_shells_debugging_component=new OPENGL_COMPONENT_THIN_SHELLS_DEBUGGING_2D<T>(grid,basedir);
        opengl_world.Set_Key_Binding_Category("Thin Shells");
        Add_Component(thin_shells_debugging_component,"thin shells debugging",'\0',BASIC_VISUALIZATION::OWNED);
        opengl_world.Append_Bind_Key(OPENGL_KEY(OPENGL_KEY::F4),thin_shells_debugging_component->Toggle_Draw_Grid_Visibility_CB());
        opengl_world.Append_Bind_Key(OPENGL_KEY(OPENGL_KEY::F5),thin_shells_debugging_component->Toggle_Draw_Density_Valid_Mask_CB());
        opengl_world.Append_Bind_Key(OPENGL_KEY(OPENGL_KEY::F5),thin_shells_debugging_component->Toggle_Draw_Phi_Valid_Mask_CB());}

    filename=basedir+"/%d/strain";
    if(has_valid_grid && FILE_UTILITIES::Frame_File_Exists(filename,start_frame)){
        OPENGL_COMPONENT_SYMMETRIC_MATRIX_FIELD_2D<T>* strain_component=new OPENGL_COMPONENT_SYMMETRIC_MATRIX_FIELD_2D<T>(grid,filename);
        Add_Component(strain_component,"Strain",'e',BASIC_VISUALIZATION::OWNED|BASIC_VISUALIZATION::START_HIDDEN);
        opengl_world.Append_Bind_Key('+',strain_component->Increase_Size_CB());
        opengl_world.Append_Bind_Key('_',strain_component->Decrease_Size_CB());}

    filename=basedir+"/%d/strain_1";
    if(has_valid_grid && FILE_UTILITIES::Frame_File_Exists(filename,start_frame)){
        OPENGL_COMPONENT_SYMMETRIC_MATRIX_FIELD_2D<T>* strain_component=new OPENGL_COMPONENT_SYMMETRIC_MATRIX_FIELD_2D<T>(grid,filename);
        Add_Component(strain_component,"Strain",'e',BASIC_VISUALIZATION::OWNED|BASIC_VISUALIZATION::START_HIDDEN);
        opengl_world.Append_Bind_Key('+',strain_component->Increase_Size_CB());
        opengl_world.Append_Bind_Key('_',strain_component->Decrease_Size_CB());}

    // Commands and misc
    opengl_world.Set_Key_Binding_Category("Misc.");
    opengl_world.Append_Bind_Key('~',Command_Prompt_CB("Command prompt"));
    opengl_world.Append_Bind_Key('\\',Toggle_2d_Mode_CB("Cycle 2d mode."));

    // Initialize priority ordering for selections
    Selection_Priority(OPENGL_SELECTION::DEBUG_PARTICLES_2D)=110;
    Selection_Priority(OPENGL_SELECTION::POINTS_2D)=100;
    Selection_Priority(OPENGL_SELECTION::COMPONENT_PARTICLES_2D)=100;
    Selection_Priority(OPENGL_SELECTION::ARTICULATED_RIGID_BODIES_JOINT_2D)=95;
    Selection_Priority(OPENGL_SELECTION::ARTICULATED_RIGID_BODIES_MUSCLE_2D)=95;
    Selection_Priority(OPENGL_SELECTION::COMPONENT_RIGID_BODIES_2D)=80;
    Selection_Priority(OPENGL_SELECTION::COMPONENT_DEFORMABLE_OBJECT_2D)=80;
    Selection_Priority(OPENGL_SELECTION::SEGMENTED_CURVE_VERTEX_2D)=79;
    Selection_Priority(OPENGL_SELECTION::SEGMENTED_CURVE_SEGMENT_2D)=78;
    Selection_Priority(OPENGL_SELECTION::TRIANGULATED_AREA_VERTEX)=77;
    Selection_Priority(OPENGL_SELECTION::TRIANGULATED_AREA_SEGMENT)=76;
    Selection_Priority(OPENGL_SELECTION::TRIANGULATED_AREA_TRIANGLE)=75;
    Selection_Priority(OPENGL_SELECTION::GRID_NODE_2D)=70;
#ifndef COMPILE_WITHOUT_DYADIC_SUPPORT
    Selection_Priority(OPENGL_SELECTION::QUADTREE_NODE)=70;
#endif
    Selection_Priority(OPENGL_SELECTION::GRID_CELL_2D)=60;
#ifndef COMPILE_WITHOUT_RLE_SUPPORT
    Selection_Priority(OPENGL_SELECTION::RLE_CELL_2D)=60;
#endif
#ifndef COMPILE_WITHOUT_DYADIC_SUPPORT
    Selection_Priority(OPENGL_SELECTION::QUADTREE_CELL)=60;
#endif
}
//#####################################################################
// Add_OpenGL_Initialization
//#####################################################################
template<class T,class RW> void OPENGL_2D_VISUALIZATION<T,RW>::
Add_OpenGL_Initialization()
{
    ANIMATED_VISUALIZATION::Add_OpenGL_Initialization();
    opengl_world.Set_2D_Mode(true);
}
//#####################################################################
// Pre_Frame_Extra
//#####################################################################
template<class T,class RW> void OPENGL_2D_VISUALIZATION<T,RW>::
Pre_Frame_Extra()
{
#ifndef COMPILE_WITHOUT_RLE_SUPPORT
    if(rle_grid_component) rle_grid_component->Set_Frame(frame);
#endif
#ifndef COMPILE_WITHOUT_DYADIC_SUPPORT
    if(quadtree_grid_component) quadtree_grid_component->Set_Frame(frame);
#endif
    if(grid_component) grid_component->object.Set_Frame(frame);
}
//#####################################################################
// Set_Frame_Extra
//#####################################################################
template<class T,class RW> void OPENGL_2D_VISUALIZATION<T,RW>::
Set_Frame_Extra()
{
    std::string filename=STRING_UTILITIES::string_sprintf("%s/%d/frame_title",basedir.c_str(),frame);
    if(FILE_UTILITIES::File_Exists(filename)){std::ifstream input(filename.c_str());getline(input,frame_title);}
    else frame_title="";
    filename=STRING_UTILITIES::string_sprintf("%s/%d/time",basedir.c_str(),frame);
    if(FILE_UTILITIES::File_Exists(filename)){T time;FILE_UTILITIES::template Read_From_File<RW>(filename,time);frame_title=STRING_UTILITIES::string_sprintf("(%.05f) ",time)+frame_title;}
}
//#####################################################################
// Command_Prompt_Response
//#####################################################################
template<class T,class RW> void OPENGL_2D_VISUALIZATION<T,RW>::
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
                // TODO : handle multiple particle sets
                if(selection){Set_Current_Selection(selection);Update_OpenGL_Strings();}}}
        else if(command=="p"){
            std::string next_word;
            if(sstream>>next_word){
                if(next_word=="auto"){
                    if(pressure_component){
                        T min_pressure,max_pressure;
                        min_pressure=pressure_component->opengl_scalar_field.values.Min();
                        max_pressure=pressure_component->opengl_scalar_field.values.Max();
                        LOG::cout<<"min pressure="<<min_pressure<<", max pressure="<<max_pressure<<std::endl;
                        pressure_component->opengl_scalar_field.Set_Scale_Range(min_pressure,max_pressure);
                        pressure_component->opengl_scalar_field.Update();}}
                else{
                    T min_pressure=atof(next_word.c_str()),max_pressure;
                    if(sstream>>max_pressure){
                        if(pressure_component){
                            LOG::cout<<"min pressure="<<min_pressure<<", max pressure="<<max_pressure<<std::endl;
                            pressure_component->opengl_scalar_field.Set_Scale_Range(min_pressure,max_pressure);
                            pressure_component->opengl_scalar_field.Update();}}}}}
        else if(command=="skip"){
            int skip;
            sstream>>skip;
            frame_increment=skip;}
        else if(command=="sphcw"){
            if(sph_cell_weight_component){
                std::ostringstream output_stream;int x,y;
                if(sstream>>x && sstream>>y){
                    VECTOR<int,2> index(x,y);
                    if(sph_cell_weight_component->Is_Up_To_Date(frame))
                        if(sph_cell_weight_component->opengl_scalar_field.values.Valid_Index(index))
                            output_stream<<"SPH cell weight at "<<x<<" "<<y<<": "<<sph_cell_weight_component->opengl_scalar_field.values(index)<<std::endl;
                        else output_stream<<x<<" "<<y<<" is not a valid cell index"<<std::endl;
                    else output_stream<<"frame "<<frame<<" is not up to date"<<std::endl;}
                else output_stream<<"invalid input"<<std::endl;
                opengl_world.Add_String(output_stream.str());}}}
}
//#####################################################################
// Toggle_2d_Mode
//#####################################################################
template<class T,class RW> void OPENGL_2D_VISUALIZATION<T,RW>::
Toggle_2d_Mode()
{
    opengl_world.Set_2D_Mode(!opengl_world.mode_2d);
}
//#####################################################################
// Command_Prompt
//#####################################################################
template<class T,class RW> void OPENGL_2D_VISUALIZATION<T,RW>::
Command_Prompt()
{
    opengl_world.Prompt_User("Command: ",Command_Prompt_Response_CB());
}
//#####################################################################
// main
//#####################################################################
int main(int argc,char* argv[])
{
    Initialize_Geometry_Particle();
    Initialize_Read_Write_Structures();
    PROCESS_UTILITIES::Set_Floating_Point_Exception_Handling(false);
    PROCESS_UTILITIES::Set_Backtrace(true);
    bool type_double=false; // float by default
    if(PARSE_ARGS::Find_And_Remove("-float",argc,argv)) type_double=false;
    if(PARSE_ARGS::Find_And_Remove("-double",argc,argv)) type_double=true;

    ANIMATED_VISUALIZATION* visualization=0;
    if(!type_double) visualization=new OPENGL_2D_VISUALIZATION<float>;
#ifndef COMPILE_WITHOUT_DOUBLE_SUPPORT
    else visualization=new OPENGL_2D_VISUALIZATION<double>;
#else
    else{std::cerr<<"Double support not enabled."<<std::endl;exit(1);}
#endif
    visualization->Initialize_And_Run(argc,argv);

    return 0;
}

template OPENGL_COMPONENT_PARTICLES_2D<float,GEOMETRY_PARTICLES<VECTOR<float,2> >,float>::OPENGL_COMPONENT_PARTICLES_2D(std::basic_string<char,std::char_traits<char>,std::allocator<char> > const&,std::basic_string<char,std::char_traits<char>,std::allocator<char> > const&,bool,bool);
#ifndef COMPILE_WITHOUT_DOUBLE_SUPPORT
template OPENGL_COMPONENT_PARTICLES_2D<double,GEOMETRY_PARTICLES<VECTOR<double,2> >,double>::OPENGL_COMPONENT_PARTICLES_2D(std::basic_string<char,std::char_traits<char>,std::allocator<char> > const&,std::basic_string<char,std::char_traits<char>,std::allocator<char> > const&,bool,bool);
#endif
