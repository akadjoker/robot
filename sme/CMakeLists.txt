project(sme)
cmake_policy(SET CMP0072 NEW)


set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ")

if (WIN32)
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS}   -D_CRT_SECURE_NO_WARNINGS")
    if (MSVC)
        if(CMAKE_BUILD_TYPE MATCHES Debug)
            add_compile_options(/RTC1 /Od /Zi)
            set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} /fsanitize=address")
        endif()     
    endif()

endif()

set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_SOURCE_DIR}/bin)

add_compile_options(

)


file(GLOB SOURCES "src/*.cpp")
add_executable(sme   ${SOURCES})


target_include_directories(sme  PUBLIC  
include 
src 
raylib 
ode
ode/include
/usr/include/python3.10
)



 

if(CMAKE_BUILD_TYPE MATCHES Debug)

if (UNIX)
#target_compile_options(sme PRIVATE -fsanitize=address -fsanitize=undefined -fsanitize=leak -g  -D_DEBUG -DVERBOSE)
#target_link_options(sme PRIVATE -fsanitize=address -fsanitize=undefined -fsanitize=leak -g  -D_DEBUG) 
endif()


elseif(CMAKE_BUILD_TYPE MATCHES Release)
    target_compile_options(sme PRIVATE -O3   -DNDEBUG )
    target_link_options(sme PRIVATE -O3   -DNDEBUG )
endif()

target_link_libraries(sme raylib ODE python3.10 )

if (WIN32)
    target_link_libraries(sme Winmm.lib)
endif()


if (UNIX)
    target_link_libraries(sme  GL m )
endif()
