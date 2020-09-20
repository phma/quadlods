# Overview
Quadlods is a library for making Richtmyer low-discrepancy sequences using quadratic irrationals, which are approximated with multiple-precision rationals. Quadlods also generates Halton low-discrepancy sequences.

Quadlods can generate sequences up to 6542 dimensions (the number of primes less than 65536), but some pairs of primes give particularly bad distributions with high discrepancy. Each prime corresponds to a quadratic irrational (quad for short). The primes are sorted by the maximum continued fraction term of their quads. Use at most 3624 dimensions to avoid the bad pair (13691,13693).

Quadlods scrambles sequences using any of these scrambling methods: third, Thue-Morse, and Gray for Richtmyer, and power, Faure, and tipwitch for Halton. You can also specify none or default.

# Compiling
To compile, if you're not developing the program:

1. Create a subdirectory build/ inside the directory where you untarred the source code.
2. `cd build`
3. `cmake ..`
4. `make`

If you are developing the program:

1. Create a directory build/quadlods outside the directory where you cloned the source code.
2. `cd build/quadlods`
3. `cmake <directory where the source code is>`
4. `make`

It makes a file called `primes.dat`, which contains the primes in order of their continued fraction terms, and a file `permute.dat`, which contains the permutations for tipwitch-scrambling Halton, and installs them in `share/quadlods/`.

GMP is required both to compile and to run the library. Boost is required to run the program, but not the library. To compile your own code that uses the library, first run "make install" as root. If you are using CMake, put the provided file FindQuadlods.cmake in your project's cmake/Modules directory, along with the GMP and GMPXX modules, and add the following lines to your CMakeLists.txt file:

set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${CMAKE_SOURCE_DIR}/cmake/Modules/")
find_package(Quadlods REQUIRED)
set(LIBS ${LIBS} ${Quadlods_LIBRARY} ${GMP_LIBRARY} ${GMPXX_LIBRARIES})
target_link_libraries(<project> ${LIBS})
include_directories(${PROJECT_BINARY_DIR} ${Quadlods_INCLUDE_DIR})

A program called "quadlods" accompanies the library. With it, you can test a set of primes for suitability. There are three tests: scatter, circle, and fill. Scatter and circle test pairs of primes; fill shows how the tuples fill s-dimensional space.

# Commands
`quadlods sortprimes` and `quadlods scramblehalton` produce the program's data files.

`quadlods test` and `quadlods longtest` test the program.

`quadlods scatter` produces a scatter plot for each pair of quads; a 10-dimensional set produces 45 pages. Recommended number of iterations is 30000, with the plot getting too dark after 100000. You can see the effect of scrambling.

`quadlods circle` plots the error of finding the area of two quarter-circles with each pair of quads. Again, the number of pages is a triangular number. Good pairs have a variety of forms and the error does not exceed 21.25 and is usually less than 5; bad pairs have a graph like one side of a trumpet, with a large error. Recommended number of iterations is 1048576.

`quadlods fill` picks s points in s-space at random without replacement and plots the statistics of vectors from each to the nearest generated point. The first graph is deprecated, but I won't remove it until the second is calibrated. The second plots the size of a ball centered at a random point and just touching the nearest low-discrepancy point. The third plots how well the directions from the random points to the nearest low-discrepancy points are distributed on the s-1-sphere. It is normalized so that it should be 1 on average. Recommended number of iterations is 1048576, unless s is high. For high s (like at least 30), the left side of the graphs is unreliable, as most random points are 1/512 or 511/512 in at least one coordinate, and no LD point has yet fallen between the random point and the nearest wall.

`quadlods flower` draws a flower plot of each component of a sequence. Recommended number of iterations is 30000.

`quadlods discrepancy` computes a lower bound of the discrepancy of a sequence. If two runs on the same sequence give the same number, and no run gives a larger number, it's probably the true discrepancy.

`quadlods textout` outputs a sequence in text.

`quadlods interact` enters interactive mode, which can be used as an Internet server using `xinetd` or by programs in any language.
