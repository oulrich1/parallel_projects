
solution "Matrix Multiply Solution"

    root = "."

    configurations {"Release", "Debug"}
    
    project "matrixmult"
        kind "ConsoleApp"
        language "C"
        files {"matrixmult.c"}
        buildoptions {"-std=c99", "-g",  "-pg", "-O3", "-funroll-loops", "-march=native"} --, "-Werror"} -Wall
        includedirs { "$HOME/opt/openmpi/include" }
        libdirs {"/usr/local/lib", "$HOME/opt/openmpi/lib"}
        location "build"
        targetname "matrixmult"
        
        configuration "Release"
            links {"pthread", "gsl", "gslcblas", "m"}

    
