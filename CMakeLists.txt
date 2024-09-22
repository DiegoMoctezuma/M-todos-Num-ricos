cmake_minimum_required(VERSION 3.30.0)
project(programaUno LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 20)

include(FetchContent)
FetchContent_Declare(fxtui
    GIT_REPOSITORY "https://github.com/ArthurSonzogni/FTXUI.git"
    GIT_TAG "main"
    OVERRIDE_FIND_PACKAGE
)
FetchContent_Declare(ftxui-grid-container
    GIT_REPOSITORY "https://github.com/mingsheng13/grid-container-ftxui"
    GIT_TAG "master"
    OVERRIDE_FIND_PACKAGE
)
FetchContent_MakeAvailable(fxtui ftxui-grid-container)
add_library(ftxui::grid-container ALIAS ftxui-grid-container)

add_executable(programaUno "main.cpp")
target_link_libraries(programaUno
    ftxui::screen ftxui::dom ftxui::component ftxui::grid-container
)
target_include_directories(programaUno PRIVATE
    "programaUno"
)