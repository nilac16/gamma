add_library(gamma)
add_library(gamma::gamma ALIAS gamma)

target_sources(gamma
    PRIVATE
        gamma.c
        distribution.c
        psearch.c
        mat.c)

target_link_libraries(gamma PUBLIC m)
