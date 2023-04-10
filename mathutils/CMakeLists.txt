add_library(
    mathutils STATIC 
    src/functions.cpp 
    src/linalg.cpp 
    src/interpolation.cpp
)
target_include_directories(
    mathutils PUBLIC include
)
target_link_libraries(
    mathutils PRIVATE GSL::gsl
)