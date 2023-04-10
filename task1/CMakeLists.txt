set(CMAKE_AUTOUIC ON)
set(CMAKE_AUTOMOC ON)
set(CMAKE_AUTORCC ON)

find_package(Qt5 COMPONENTS Core Widgets Charts REQUIRED)

add_executable(
    task1 
    src/main.cpp
    src/solution.cpp
    src/window.cpp
)
target_include_directories(
    task1 PRIVATE include
)
target_link_libraries(          
    task1 PUBLIC mathutils Qt5::Widgets Qt5::Charts
)