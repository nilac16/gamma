if (NOT MSVC)
    set(WARN_FLAGS "-W -Wall -Wextra -Werror")
    set(ARCH_FLAGS "-march=native")
    set(BUILD_FLAGS "${WARN_FLAGS} ${ARCH_FLAGS}")
else ()
    set(BUILD_FLAGS "/W3 /WX /arch:AVX2 /Zc:__cplusplus")
endif ()

set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${BUILD_FLAGS}")
