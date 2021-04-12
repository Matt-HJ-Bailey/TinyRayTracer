# TinyRayTracer
A C++ implementation of a path tracer to render scenes of spheres

Reads in scene data from `spheres.inpt`, which contains the following data:

`// Radius, Position, Emission, Colour, Type, Texture
1e5  , 100001 40.8 81.6 , 0 0 0 , 0.75 0.25 0.25 , DIFFUSE , NONE //LEFT
1e5  , -99901 40.8 81.6 , 0 0 0 , 0.25 0.25 0.75 , DIFFUSE , NONE //RIGHT
1e5  , 50 40.8 1e5      , 0 0 0 , 0.75 0.75 0.75 , DIFFUSE , NONE //BACK
1e5  , 50 40.8 -99830   , 0 0 0 , 0 0 0          , DIFFUSE , NONE //FRONT
1e5  , 50 1e5 81.6      , 0 0 0 , 0.75 0.75 0.75 , DIFFUSE , NONE //BOTTOM
1e5  , 50 -99918.4 81.6  , 0 0 0 , 0.75 0.75 0.75 , DIFFUSE , NONE //TOP
16.5 , 27 16.5 47       , 0.0 0.0 0.0, 0.99 0.99 0.99 , SPECULAR , mars.tga // SHINY1
16.5 , 73 16.5 78       , 0.0 0.0 0.0, 1 1 1 , DIFFUSE , earth.tga // SHINY2
600  , 50 681.33 81.6   , 12  12  12 , 0 0 0 , DIFFUSE , NONE // LIGHT`
