// Copyright (c) 2011, Eftychios Sifakis.
// Distributed under the FreeBSD license (see license.txt)

#include <PhysBAM_Geometry/Geometry_Particles/REGISTER_GEOMETRY_READ_WRITE.h>
#include <PhysBAM_Geometry/Solids_Geometry/DEFORMABLE_GEOMETRY_COLLECTION.h>

#include "SIMULATION_LAYOUT.h"
#include "SIMULATION_DRIVER.h"

using namespace PhysBAM;

int main(int argc,char* argv[])
{
    typedef float T;typedef float RW;typedef VECTOR<T,3> TV;
    RW rw=RW();STREAM_TYPE stream_type(rw);

    LOG::Initialize_Logging();

    SIMULATION_LAYOUT<T> layout(stream_type);
    SIMULATION_DRIVER<T> driver(layout);
    driver.Run();

    LOG::Finish_Logging();
} 
