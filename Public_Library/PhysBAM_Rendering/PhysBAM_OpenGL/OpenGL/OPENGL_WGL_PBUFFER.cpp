//#####################################################################
// Copyright 2004, Eran Guendelman.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#include <PhysBAM_Rendering/PhysBAM_OpenGL/OpenGL/OPENGL_WGL_PBUFFER.h>

#ifdef OPENGL_WGL_PBUFFER_SUPPORTED
#include <iostream>


using namespace PhysBAM;

// pbuffer function prototypes
static PFNWGLMAKECONTEXTCURRENTARBPROC wglMakeContextCurrentARB = 0;
static PFNWGLCHOOSEPIXELFORMATARBPROC   wglChoosePixelFormatARB = 0;
static PFNWGLCREATEPBUFFERARBPROC       wglCreatePbufferARB = 0;
static PFNWGLDESTROYPBUFFERARBPROC     wglDestroyPbufferARB = 0;
static PFNWGLGETPBUFFERDCARBPROC           wglGetPbufferDCARB = 0;
static PFNWGLRELEASEPBUFFERDCARBPROC    wglReleasePbufferDCARB = 0;
static PFNWGLQUERYPBUFFERARBPROC           wglQueryPbufferARB = 0;
static PFNWGLBINDTEXIMAGEARBPROC       wglBindTexImageARB = 0;
static PFNWGLRELEASETEXIMAGEARBPROC    wglReleaseTexImageARB = 0;
static PFNWGLSETPBUFFERATTRIBARBPROC   wglSetPbufferAttribARB = 0;
static PFNWGLGETPIXELFORMATATTRIBIVARBPROC wglGetPixelFormatAttribivARB = 0;
static PFNWGLGETPIXELFORMATATTRIBFVARBPROC wglGetPixelFormatAttribfvARB = 0;

LRESULT CALLBACK MainWndProc(          HWND hwnd,
                            UINT uMsg,
                            WPARAM wParam,
                            LPARAM lParam
                            )
{
    return DefWindowProc(hwnd,uMsg,wParam,lParam);
}



bool OPENGL_PBUFFER::Create(int width, int height)
{
    // record size
    pbufferWidth=width;
    pbufferHeight=height;
    // make a window
    WNDCLASSEX wndclass;

    // Initialize the entire structure to zero. 
    char szMainWndClass[]="stupid_offscreen_reg_class";
    memset (&wndclass, 0, sizeof(WNDCLASSEX));
    wndclass.lpszClassName = szMainWndClass;
    wndclass.cbSize = sizeof(WNDCLASSEX);
    wndclass.style = CS_HREDRAW | CS_VREDRAW;
    wndclass.lpfnWndProc = MainWndProc;
    wndclass.hInstance = GetModuleHandle(0);
    wndclass.hIcon = LoadIcon (0, IDI_APPLICATION);
    wndclass.hIconSm = LoadIcon (0, IDI_APPLICATION);
    wndclass.hCursor = LoadCursor (0, IDC_ARROW);
    wndclass.hbrBackground = (HBRUSH) GetStockObject (WHITE_BRUSH);
    RegisterClassEx (&wndclass);
    LOG::cout<<"done"<<std::endl;
    hwnd= CreateWindow (
        szMainWndClass,             // Class name 
        "Hello",                    // Caption 
        WS_OVERLAPPEDWINDOW,        // Style 
        CW_USEDEFAULT,              // Initial x (use default) 
        CW_USEDEFAULT,              // Initial y (use default) 
        CW_USEDEFAULT,              // Initial x size (use default) 
        CW_USEDEFAULT,              // Initial y size (use default) 
        0,                       // No parent window 
        0,                       // No menu 
        GetModuleHandle(0),                       //This program instance 
        0                         //Creation parameters 
        );
    HDC dc=GetDC(hwnd);
    // give our new window a gl context
    PIXELFORMATDESCRIPTOR pfd;unsigned int selected_pf;
    pfd.nSize = sizeof( PIXELFORMATDESCRIPTOR );
    pfd.nVersion = 1;
    pfd.dwFlags = PFD_DRAW_TO_WINDOW | PFD_SUPPORT_OPENGL |
        PFD_DOUBLEBUFFER | PFD_TYPE_RGBA;
    pfd.cColorBits = 0;
    pfd.cRedBits = 0;
    pfd.cGreenBits = 0;
    pfd.cBlueBits = 0;
    pfd.cRedShift = 0;
    pfd.cGreenShift = 0;
    pfd.cBlueShift = 0;
    pfd.cAlphaBits = 0;
    pfd.cAccumAlphaBits = 0;
    pfd.cAccumBits = 0;
    pfd.cAccumBlueBits = 0;
    pfd.cAccumRedBits = 0;
    pfd.cAccumGreenBits = 0;
    pfd.cDepthBits = 16;
    pfd.cStencilBits = 0;
    pfd.cAuxBuffers = 0;
    pfd.iLayerType = PFD_MAIN_PLANE;
    pfd.bReserved = 0;
    pfd.dwLayerMask = 0;
    pfd.dwVisibleMask = 0;
    pfd.dwDamageMask = 0;
    if(!(selected_pf = ChoosePixelFormat( dc, &pfd )))goto fail;
    if(!SetPixelFormat(dc,selected_pf,&pfd)) goto fail;
    hGLRC=wglCreateContext(dc);
    wglMakeCurrent(dc,hGLRC);
    //ShowWindow(hwnd, SW_SHOW);
    hGLDC=wglGetCurrentDC();

    // check for extension and then get entry points entry points
    const char *str;
    PFNWGLGETEXTENSIONSSTRINGARBPROC wglGetExtensionsStringARB =(PFNWGLGETEXTENSIONSSTRINGARBPROC)wglGetProcAddress("wglGetExtensionsStringARB");
    if(!wglGetExtensionsStringARB)goto fail;
    str=wglGetExtensionsStringARB(hGLDC);
    if(!strstr(str,"WGL_ARB_pixel_format"))goto fail;
    if(!strstr(str,"WGL_ARB_pbuffer"))goto fail;
    
#define INIT_ENTRY_POINT( funcname, type )\
    funcname = (type) wglGetProcAddress(#funcname);\
    if ( !funcname )\
    LOG::cerr<<"#funcname() not initialized"<<std::endl;

    // Initialize WGL_ARB_pbuffer entry points. 
    INIT_ENTRY_POINT(wglCreatePbufferARB, PFNWGLCREATEPBUFFERARBPROC);
    INIT_ENTRY_POINT(wglGetPbufferDCARB, PFNWGLGETPBUFFERDCARBPROC);
    INIT_ENTRY_POINT(wglReleasePbufferDCARB, PFNWGLRELEASEPBUFFERDCARBPROC);
    INIT_ENTRY_POINT(wglDestroyPbufferARB, PFNWGLDESTROYPBUFFERARBPROC);
    INIT_ENTRY_POINT(wglQueryPbufferARB, PFNWGLQUERYPBUFFERARBPROC);
    // Initialize WGL_ARB_pixel_format entry points. 
    INIT_ENTRY_POINT(wglGetPixelFormatAttribivARB,PFNWGLGETPIXELFORMATATTRIBIVARBPROC);
    INIT_ENTRY_POINT(wglGetPixelFormatAttribfvARB,PFNWGLGETPIXELFORMATATTRIBFVARBPROC);
    INIT_ENTRY_POINT(wglChoosePixelFormatARB,PFNWGLCHOOSEPIXELFORMATARBPROC);


    LOG::cout<<"Made window and render context" << std::endl;
    
    // Get the pixel format
    unsigned int count = 0;
    int pixelFormat;
    int attr[] = {
        WGL_SUPPORT_OPENGL_ARB,         TRUE,
        WGL_DRAW_TO_PBUFFER_ARB,        TRUE,
        WGL_RED_BITS_ARB,               8,
        WGL_GREEN_BITS_ARB,             8,
        WGL_BLUE_BITS_ARB,              8,
        WGL_ALPHA_BITS_ARB,             8,
        WGL_DEPTH_BITS_ARB,             24,
        WGL_DOUBLE_BUFFER_ARB,          FALSE,
        0
    };
    wglChoosePixelFormatARB(hGLDC, (const int*)attr, 0, 1, &pixelFormat, &count);
    if(count==0)goto fail;
    int pAttrib[] ={0};
    // create the pbuffer 
    if(!(hPBuffer = wglCreatePbufferARB(hGLDC, pixelFormat, width, height, pAttrib))) goto post_gl_init_fail;
    if(!(hPBufferDC=wglGetPbufferDCARB(hPBuffer))) goto post_gl_init_fail;
    if(!(hPBufferGLRC= wglCreateContext(hPBufferDC))) goto post_gl_init_fail;
    if(!wglMakeCurrent(hPBufferDC,hPBufferGLRC)) goto post_gl_init_fail;
    // At this point it worked
    return true;
    // otherwise we jumped to failure cases
post_gl_init_fail:
    wglDeleteContext(hGLRC);
fail:
    DestroyWindow(hwnd);
    return false;
}

void OPENGL_PBUFFER::Destroy()
{
    wglReleaseTexImageARB(hPBuffer, WGL_FRONT_LEFT_ARB);
    wglDestroyPbufferARB(hPBuffer);
    wglDeleteContext(hPBufferGLRC);
    wglDeleteContext(hGLRC);
    DestroyWindow(hwnd);
}
#endif
