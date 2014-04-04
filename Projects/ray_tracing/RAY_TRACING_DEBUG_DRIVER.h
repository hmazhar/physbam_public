//#####################################################################
// Copyright 2002-2006, Ronald Fedkiw, Geoffrey Irving, Duc Nguyen, Andrew Selle.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class RAY_TRACING_DRIVER  
//##################################################################### 
#ifndef __RAY_TRACING_DEBUG_DRIVER__
#define __RAY_TRACING_DEBUG_DRIVER__
#include <PhysBAM_Rendering/PhysBAM_Ray_Tracing/Rendering/RAY_TRACER_DEBUG_DATA.h>
#include <PhysBAM_Rendering/PhysBAM_Ray_Tracing/Rendering/RENDER_WORLD.h>
#include <fstream>
#include <iostream>
#include <typeinfo>
#include "Ray_Debugger/Fl_ToggleTree.h"
#include "Ray_Debugger/IMAGE_WINDOW.h"
#include "Ray_Debugger/SCENE_WINDOW.h"
#include "RAY_TRACING_EXAMPLE.h"
#include <FL/Fl.h>
#include <FL/Fl_Browser.h>
#include <Fl/Fl_File_Chooser.H>
#include <Fl/Fl_Hold_Browser.H>
#include <FL/Fl_Menu_Bar.H>
#include <FL/Fl_Menu_Item.H>
#include <FL/Fl_Tile.h>
namespace PhysBAM{

template<class T>
class RAY_TRACING_DEBUG_DRIVER
{
public:
    RAY_TRACING_EXAMPLE<T>& example;
    RENDER_WORLD<T> world;
    Fl_Window *main_window;
    IMAGE_WINDOW<T> *image_window;
    SCENE_WINDOW<T> *scene_window;
    Fl_Browser *browser_window;
    Fl_ToggleTree *tree;
    Fl_Hold_Browser *object_list;
    int update_count,current_i,current_j,passes,current_pass;
    bool rendering;
    T one_over_gamma;
    
    RAY_TRACING_DEBUG_DRIVER(RAY_TRACING_EXAMPLE<T>& example_input)
        :example(example_input),main_window(0),image_window(0),scene_window(0),browser_window(0),tree(0),object_list(0),update_count(0),
         current_i(1),current_j(1),passes(5),current_pass(0),rendering(true),one_over_gamma((T)1/example.gamma_correction)
    {
        world.Track_Debugging_Information();
    }

    VECTOR<T,3> Gamma_Correct(VECTOR<T,3>& color)
    {return VECTOR<T,3>(pow(color.x,one_over_gamma),pow(color.y,one_over_gamma),pow(color.z,one_over_gamma));}

    void One_Time_Initialization()
    {example.Initialize_Scene(world,example.frame);

    // determine pixels to render
    RANGE<VECTOR<int,2> >& clip=example.clipping_region;
    clip.min_corner.x=max(clip.min_corner.x,1);clip.max_corner.x=min(clip.max_corner.x,world.camera.film.grid.m);
    clip.min_corner.y=max(clip.min_corner.y,1);clip.max_corner.y=min(clip.max_corner.y,world.camera.film.grid.n);}

    void Setup_Scene(int frame)
    {
    // update interface and scene window
    scene_window->Initialize_Objects();
    if(object_list){
        object_list->clear();
        for(int i=1;i<=world.standard_objects.m;i++)
            object_list->add(world.standard_objects(i)->name.c_str(),world.standard_objects(i));}}

    void Create_Windows_And_Widgets()
    {Fl::scheme("plastic");
    main_window=new Fl_Window(1280,990,"Ray Debugger");main_window->resizable(main_window);
    Fl_Menu_Item menu_items[] = {
        {"&File",       0,0,0,FL_SUBMENU},
        {"&Save SGI RGB",FL_CTRL+'s',(Fl_Callback *)RAY_TRACING_DEBUG_DRIVER<T>::Save_RGB_Callback,this,FL_MENU_DIVIDER},
        {"&Save PhysBAM Float Image",FL_CTRL+'e',(Fl_Callback *)RAY_TRACING_DEBUG_DRIVER<T>::Save_PhysBAM_Callback,this,FL_MENU_DIVIDER},
        {"E&xit",FL_CTRL+'q',(Fl_Callback *)RAY_TRACING_DEBUG_DRIVER<T>::Quit_Callback,0},
        {0},
        {"&View",       0,0,0,FL_SUBMENU},
        {"&Global Photons",FL_CTRL+'g',(Fl_Callback *)RAY_TRACING_DEBUG_DRIVER<T>::Toggle_View_Global_Photons,this,FL_MENU_TOGGLE},
        {"&Caustic Photons",FL_CTRL+'c',(Fl_Callback *)RAY_TRACING_DEBUG_DRIVER<T>::Toggle_View_Caustic_Photons,this,FL_MENU_TOGGLE},
        {"&Volume Photons",0,(Fl_Callback *)RAY_TRACING_DEBUG_DRIVER<T>::Toggle_View_Volume_Photons,this,FL_MENU_TOGGLE},
        {"&Irradiance Cache Samples",FL_CTRL+'i',(Fl_Callback *)RAY_TRACING_DEBUG_DRIVER<T>::Toggle_View_Irradiance_Cache,this,FL_MENU_TOGGLE},
        {"&Used Photons",FL_CTRL+'u',(Fl_Callback *)RAY_TRACING_DEBUG_DRIVER<T>::Toggle_View_Used_Photons,this,FL_MENU_TOGGLE},
        {"Used Irradiance Cache &Samples",0,(Fl_Callback *)RAY_TRACING_DEBUG_DRIVER<T>::Toggle_View_Used_Irradiance_Cache_Samples,this,FL_MENU_TOGGLE},
        {"Used Subsur&face Scattering Samples",0,(Fl_Callback *)RAY_TRACING_DEBUG_DRIVER<T>::Toggle_View_Subsurface_Scattering_Samples,this,FL_MENU_TOGGLE},
        {0},
        {0}};
    Fl_Menu_Bar *menu_bar=new Fl_Menu_Bar(0,0,1280,30);
    menu_bar->copy(menu_items);
    Fl_Tile* tile=new Fl_Tile(0,30,1280,960);tile->resizable(tile);
    Fl_Group* group=new Fl_Group(0,30,640,640);group->resizable(group);
    Fl_Scroll* scroll=new Fl_Scroll(0,30,640,640);scroll->resizable(scroll);
    tree=new Fl_ToggleTree(0,30,640,640);tree->callback((Fl_Callback*)Click_Tree,this);tree->add_next("foo");
    group->end();
    printf("Making image window with %d %d\n",world.camera.film.grid.m,world.camera.film.grid.n);
    image_window=new IMAGE_WINDOW<T>(640,30,640,480,world.camera.film.grid.m,world.camera.film.grid.n);image_window->resizable(image_window);image_window->box(FL_NO_BOX);
    image_window->size_range(120,60,0,0);
    image_window->Set_Callback(RAY_TRACING_DEBUG_DRIVER<T>::Click_Pixel,(void*)this);    
    scene_window=new SCENE_WINDOW<T>(640,510,640,480,world);scene_window->resizable(scene_window);scene_window->box(FL_NO_BOX);
    scene_window->size_range(120,90,0,0);
    main_window->size_range(320,270,0,0);
    object_list=new Fl_Hold_Browser(0,670,640,160);object_list->resizable(object_list);object_list->add("test");
    object_list->callback((Fl_Callback*)RAY_TRACING_DEBUG_DRIVER<T>::Click_Object);object_list->user_data(this);
    browser_window=new Fl_Browser(0,830,640,160);browser_window->resizable(browser_window);browser_window->add("test");
    main_window->show();
    Fl::add_idle(Render_Pixel,this);}

    void Execute_Main_Program()
    {printf("One time initialization...\n");
    One_Time_Initialization();
    printf("Create windows...\n");
    Create_Windows_And_Widgets();
    printf("Setup Scene...\n");
    Setup_Scene(example.frame);
    printf("done\n");
    printf("Initialize Objects...\n");
    scene_window->Initialize_Objects();
    world.Prepare_For_Forward_Ray_Trace();
    Fl::run();}

    void Populate_Ray_Tree(Fl_ToggleTree* tree,RENDERING_RAY_DEBUG<T>* ray,int depth=0)
    {char buf[1024];
    char *type=0;
    if(ray->ray.ray_type==RENDERING_RAY<T>::UNKNOWN_RAY)type="UNKNOWN";
    if(ray->ray.ray_type==RENDERING_RAY<T>::COLOR_RAY)type="COLOR";
    if(ray->ray.ray_type==RENDERING_RAY<T>::SHADOW_RAY)type="SHADOW";
    if(ray->ray.ray_type==RENDERING_RAY<T>::DUMMY_RAY)type="DUMMY";
    if(ray->ray.ray_type==RENDERING_RAY<T>::PHOTON_GATHER_RAY)type="PHOTON_GATHER";
    if(ray->ray.ray_type==RENDERING_RAY<T>::PHOTON_RAY)type="PHOTON";
    sprintf(buf,"%s hit_object=%s semi_infinite=%s aggregate_id=%d",type,ray->hit_object?"true":"false",ray->ray.ray.semi_infinite?"true":"false",ray->ray.ray.aggregate_id);
    tree->add_sub(buf,1,0,(void*)ray);
    for(int i=1;i<=ray->children.m;i++){
        Populate_Ray_Tree(tree,ray->children(i),depth+1);
    }
    tree->traverse_up();}
    static void Click_Pixel(void *data,int x,int y)
    {RAY_TRACING_DEBUG_DRIVER<T>* driver=(RAY_TRACING_DEBUG_DRIVER<T>*)data;
    driver->tree->clear();
    RAY_TRACER_DEBUG_DATA<T> *debug_data=new RAY_TRACER_DEBUG_DATA<T>;
    driver->world.Render_Pixel(VECTOR<int,2>(x,y),debug_data);
    driver->image_window->Update_Image();
    driver->Populate_Ray_Tree(driver->tree,debug_data->ray_tree);
    driver->scene_window->Set_Debug_Ray(debug_data->ray_tree);
    debug_data->ray_tree=NULL;
    delete debug_data;}

    static void Click_Tree(void *tree,void *data)
    {RAY_TRACING_DEBUG_DRIVER<T>* driver=(RAY_TRACING_DEBUG_DRIVER<T>*)data;
    RENDERING_RAY_DEBUG<T>* ray=(RENDERING_RAY_DEBUG<T>*)(driver->tree->selected()->data());
    driver->scene_window->Select_Ray(ray);
    driver->browser_window->clear();
    Fl_Browser& out=*driver->browser_window;
    char buf[1024];char *type;
    if(ray->ray.ray_type==RENDERING_RAY<T>::UNKNOWN_RAY)type="UNKNOWN";
    else if(ray->ray.ray_type==RENDERING_RAY<T>::COLOR_RAY)type="COLOR";
    else if(ray->ray.ray_type==RENDERING_RAY<T>::SHADOW_RAY)type="SHADOW";
    else if(ray->ray.ray_type==RENDERING_RAY<T>::PHOTON_GATHER_RAY)type="PHOTON_GATHER";
    else type="Unimplemented";
    sprintf(buf,"Ray type: %s",type);out.add(buf);
    sprintf(buf,"Ray Energy: %f Ray Depth: %d",ray->ray.ray_contribution,ray->ray.recursion_depth);out.add(buf);
    VECTOR<T,3> endpoint=ray->ray.ray.endpoint;sprintf(buf,"Endpoint: %f %f %f",endpoint.x,endpoint.y,endpoint.z);out.add(buf);
    VECTOR<T,3> direction=ray->ray.ray.direction;sprintf(buf,"Direction: %f %f %f",direction.x,direction.y,direction.z);out.add(buf);
    if(ray->ray.ray.semi_infinite)out.add("Extent: semi-infinite");
    else{VECTOR<T,3> pt=ray->ray.ray.Point(ray->ray.ray.t_max);sprintf(buf,"Extent: t_max=%f point=%f %f %f",ray->ray.ray.t_max,pt.x,pt.y,pt.z);out.add(buf);}
    if(ray->hit_object){sprintf(buf,"Hit Object: aggregate_id=%d",ray->ray.ray.aggregate_id);out.add(buf);}
    // print all ray comments
    out.add("");
    for(int i=1;i<=ray->comments.m;i++)out.add(ray->comments(i).c_str());}

    static void Click_Object(void *widget)
    {Fl_Hold_Browser* browser=((Fl_Hold_Browser*)widget);
    RAY_TRACING_DEBUG_DRIVER<T>* driver=(RAY_TRACING_DEBUG_DRIVER<T>*)(browser->user_data());
    driver->scene_window->Select_Object((RENDERING_OBJECT<T>*)browser->data(browser->value()));}

    static void Quit_Callback(void *data)
    {exit(0);}

    static void Save_RGB_Callback(Fl_Menu_* item,void *data)
    {RAY_TRACING_DEBUG_DRIVER<T>* driver=(RAY_TRACING_DEBUG_DRIVER<T>*)data;
    char *filename=fl_file_chooser("Choose where to save","*.rgb",0,0);
    if(filename)driver->image_window->Save_RGB(filename);}

    static void Save_PhysBAM_Callback(Fl_Menu_* item,void *data)
    {RAY_TRACING_DEBUG_DRIVER<T>* driver=(RAY_TRACING_DEBUG_DRIVER<T>*)data;
    char *filename=fl_file_chooser("Choose where to save","*.img",0,0);
    if(filename)driver->image_window->Save_PhysBAM(filename);}

    static void Toggle_View_Global_Photons(Fl_Menu_* item,void *data)
    {RAY_TRACING_DEBUG_DRIVER<T>* driver=(RAY_TRACING_DEBUG_DRIVER<T>*)(data);
    driver->scene_window->View_Global_Photon_Map(item->mvalue()->value()!=0);}

    static void Toggle_View_Caustic_Photons(Fl_Menu_* item,void *data)
    {RAY_TRACING_DEBUG_DRIVER<T>* driver=(RAY_TRACING_DEBUG_DRIVER<T>*)(data);
    driver->scene_window->View_Caustic_Photon_Map(item->mvalue()->value()!=0);}

    static void Toggle_View_Volume_Photons(Fl_Menu_* item,void *data)
    {RAY_TRACING_DEBUG_DRIVER<T>* driver=(RAY_TRACING_DEBUG_DRIVER<T>*)(data);
    driver->scene_window->View_Volume_Photon_Map(item->mvalue()->value()!=0);}

    static void Toggle_View_Irradiance_Cache(Fl_Menu_* item,void *data)
    {RAY_TRACING_DEBUG_DRIVER<T>* driver=(RAY_TRACING_DEBUG_DRIVER<T>*)(data);
    driver->scene_window->View_Irradiance_Cache(item->mvalue()->value()!=0);}

    static void Toggle_View_Used_Photons(Fl_Menu_* item,void *data)
    {RAY_TRACING_DEBUG_DRIVER<T>* driver=(RAY_TRACING_DEBUG_DRIVER<T>*)(data);
    driver->scene_window->View_Used_Photons(item->mvalue()->value()!=0);}

    static void Toggle_View_Used_Irradiance_Cache_Samples(Fl_Menu_* item,void *data)
    {RAY_TRACING_DEBUG_DRIVER<T>* driver=(RAY_TRACING_DEBUG_DRIVER<T>*)(data);
    driver->scene_window->View_Used_Irradiance_Samples(item->mvalue()->value()!=0);}

    static void Toggle_View_Subsurface_Scattering_Samples(Fl_Menu_* item,void *data)
    {RAY_TRACING_DEBUG_DRIVER<T>* driver=(RAY_TRACING_DEBUG_DRIVER<T>*)(data);
    driver->scene_window->View_Subsurface_Scattering_Samples(item->mvalue()->value()!=0);}

    void Render_Pixel()
    {if(!rendering)return;
    if(current_i>world.camera.film.grid.m){
        current_i=1;current_j++;update_count++;
        if(update_count>2){update_count=0;image_window->Update_Image();}}
    if(current_j>world.camera.film.grid.n){
        printf("Done rendering...\n");
        image_window->Update_Image();
        rendering=false;
        return;}
    world.Render_Pixel(VECTOR<int,2>(current_i,current_j),NULL);
    current_i++;}

    static void Render_Pixel(void *data)
    {RAY_TRACING_DEBUG_DRIVER<T>* driver=(RAY_TRACING_DEBUG_DRIVER<T>*)data;
    driver->Render_Pixel();}
    //#####################################################################
};
}
#endif


