WATER_NOCONA(1)

NAME
    water_nocona - make a velocity field divergence free

SYNOPSIS
    water_nocona [OPTION]...

DESCRIPTION
    Takes a velocity field V at position x specified by V(x)=1 near the left, and phi source at position x specified by P(x)=1 near the left and runs a basic water simulation. The
    output will be dumped into a directory named output.

    -scale RESOLUTION
        Specifies the resolution of the grid used in number of cells along the y axis where as the other dimensions are specified by half that number.

    -restart FRAME
        Specifies a frame to restart the simulation from. The intended usage is when the simulation is killed or needs to be changed after a certain frame this parameter should be used to avoid
        running the entire simulation again.

    -e FRAME
        Specifies the end frame of the simulation.

    -substep LEVEL
        Specifies the amount of output written. The higher the level the more output is written and the slower the simulation is.

    -3d
        Runs the simulation in three dimensions instead of 2

    -threads NUMBER
        EXPERIMENTAL This parameter specifies the number of threads to use but is still in the development stages and may not give the desired performance.
