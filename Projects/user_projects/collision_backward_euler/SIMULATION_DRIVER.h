// Copyright (c) 2011, Eftychios Sifakis.
// Distributed under the FreeBSD license (see license.txt)

#pragma once

#include <PhysBAM_Tools/Vectors/VECTOR_3D.h>

namespace PhysBAM{

template<class T> class SIMULATION_LAYOUT;

template<class T>
class SIMULATION_DRIVER
{
public:
    typedef VECTOR<T,3> TV;

    SIMULATION_LAYOUT<T>& layout;
    T time;

    ARRAY<TV> Vp;
    
    SIMULATION_DRIVER(SIMULATION_LAYOUT<T>& layout_input)
        :layout(layout_input) {}

    void Run()
    {
        layout.Initialize();layout.Write_Output(0);time=0;

        for(int frame=1;frame<=layout.number_of_frames;frame++){
            Simulate_Frame(frame);layout.Write_Output(frame);}
    }

    void Simulate_Frame(const int frame);
    void Simulate_Time_Step(const T time,const T dt);

};

}
