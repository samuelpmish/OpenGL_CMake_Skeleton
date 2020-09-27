add_library(imgui STATIC 
  imconfig.h
  imgui_demo.cpp
  imgui_draw.cpp
  imgui_internal.h
  imgui_widgets.cpp
  imgui.h
  imgui.cpp
  imstb_rectpack.h
  imstb_textedit.h
  imstb_truetype.h
  imgui_impl_glfw.h
  imgui_impl_glfw.cpp
  imgui_impl_opengl3.h
  imgui_impl_opengl3.cpp
)

target_include_directories(imgui PUBLIC .)
target_include_directories(imgui PRIVATE ../glew/include)
target_include_directories(imgui PRIVATE ../glfw/include)
target_compile_definitions(imgui PUBLIC IMGUI_IMPL_OPENGL_LOADER_GLEW)
