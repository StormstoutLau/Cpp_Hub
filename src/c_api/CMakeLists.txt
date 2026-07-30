# C ABI library (ADR-009).
# Compiled as a static library; conditionally built when CPPHUB_ENABLE_C_API=ON.
# The header is always available at include/cpphub/c_api/cpphub_c_api.h.
add_library(cpphub_c_api STATIC c_api.cpp)
target_link_libraries(cpphub_c_api PUBLIC cpphub_core)
target_include_directories(cpphub_c_api PUBLIC
    ${CMAKE_SOURCE_DIR}/include
)

# MSVC: export symbols explicitly (static lib doesn't need dllexport, but
# consumers linking into a DLL need the symbols to be visible).
if(MSVC)
    target_compile_options(cpphub_c_api PRIVATE /utf-8 /W3)
endif()
