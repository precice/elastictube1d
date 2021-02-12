# CMAKE generated file: DO NOT EDIT!
# Generated by "Unix Makefiles" Generator, CMake Version 3.16

# Delete rule output on recipe failure.
.DELETE_ON_ERROR:


#=============================================================================
# Special targets provided by cmake.

# Disable implicit rules so canonical targets will work.
.SUFFIXES:


# Remove some rules from gmake that .SUFFIXES does not remove.
SUFFIXES =

.SUFFIXES: .hpux_make_needs_suffix_list


# Suppress display of executed commands.
$(VERBOSE).SILENT:


# A target that is always out of date.
cmake_force:

.PHONY : cmake_force

#=============================================================================
# Set environment variables for the build.

# The shell in which to execute make rules.
SHELL = /bin/sh

# The CMake executable.
CMAKE_COMMAND = /usr/bin/cmake

# The command to remove a file.
RM = /usr/bin/cmake -E remove -f

# Escaping for special characters.
EQUALS = =

# The top-level source directory on which CMake was run.
CMAKE_SOURCE_DIR = /home/homspons/Documents/hiwi/elastictube1d/fluid-cpp

# The top-level build directory on which CMake was run.
CMAKE_BINARY_DIR = /home/homspons/Documents/hiwi/elastictube1d/fluid-cpp/build

# Include any dependencies generated for this target.
include CMakeFiles/FluidSolver.dir/depend.make

# Include the progress variables for this target.
include CMakeFiles/FluidSolver.dir/progress.make

# Include the compile flags for this target's objects.
include CMakeFiles/FluidSolver.dir/flags.make

CMakeFiles/FluidSolver.dir/fluid_solver.cpp.o: CMakeFiles/FluidSolver.dir/flags.make
CMakeFiles/FluidSolver.dir/fluid_solver.cpp.o: ../fluid_solver.cpp
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/homspons/Documents/hiwi/elastictube1d/fluid-cpp/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_1) "Building CXX object CMakeFiles/FluidSolver.dir/fluid_solver.cpp.o"
	/usr/bin/c++  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/FluidSolver.dir/fluid_solver.cpp.o -c /home/homspons/Documents/hiwi/elastictube1d/fluid-cpp/fluid_solver.cpp

CMakeFiles/FluidSolver.dir/fluid_solver.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/FluidSolver.dir/fluid_solver.cpp.i"
	/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /home/homspons/Documents/hiwi/elastictube1d/fluid-cpp/fluid_solver.cpp > CMakeFiles/FluidSolver.dir/fluid_solver.cpp.i

CMakeFiles/FluidSolver.dir/fluid_solver.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/FluidSolver.dir/fluid_solver.cpp.s"
	/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /home/homspons/Documents/hiwi/elastictube1d/fluid-cpp/fluid_solver.cpp -o CMakeFiles/FluidSolver.dir/fluid_solver.cpp.s

CMakeFiles/FluidSolver.dir/fluid_nl.cpp.o: CMakeFiles/FluidSolver.dir/flags.make
CMakeFiles/FluidSolver.dir/fluid_nl.cpp.o: ../fluid_nl.cpp
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/homspons/Documents/hiwi/elastictube1d/fluid-cpp/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_2) "Building CXX object CMakeFiles/FluidSolver.dir/fluid_nl.cpp.o"
	/usr/bin/c++  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/FluidSolver.dir/fluid_nl.cpp.o -c /home/homspons/Documents/hiwi/elastictube1d/fluid-cpp/fluid_nl.cpp

CMakeFiles/FluidSolver.dir/fluid_nl.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/FluidSolver.dir/fluid_nl.cpp.i"
	/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /home/homspons/Documents/hiwi/elastictube1d/fluid-cpp/fluid_nl.cpp > CMakeFiles/FluidSolver.dir/fluid_nl.cpp.i

CMakeFiles/FluidSolver.dir/fluid_nl.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/FluidSolver.dir/fluid_nl.cpp.s"
	/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /home/homspons/Documents/hiwi/elastictube1d/fluid-cpp/fluid_nl.cpp -o CMakeFiles/FluidSolver.dir/fluid_nl.cpp.s

CMakeFiles/FluidSolver.dir/fluidComputeSolution.cpp.o: CMakeFiles/FluidSolver.dir/flags.make
CMakeFiles/FluidSolver.dir/fluidComputeSolution.cpp.o: ../fluidComputeSolution.cpp
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/homspons/Documents/hiwi/elastictube1d/fluid-cpp/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_3) "Building CXX object CMakeFiles/FluidSolver.dir/fluidComputeSolution.cpp.o"
	/usr/bin/c++  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/FluidSolver.dir/fluidComputeSolution.cpp.o -c /home/homspons/Documents/hiwi/elastictube1d/fluid-cpp/fluidComputeSolution.cpp

CMakeFiles/FluidSolver.dir/fluidComputeSolution.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/FluidSolver.dir/fluidComputeSolution.cpp.i"
	/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /home/homspons/Documents/hiwi/elastictube1d/fluid-cpp/fluidComputeSolution.cpp > CMakeFiles/FluidSolver.dir/fluidComputeSolution.cpp.i

CMakeFiles/FluidSolver.dir/fluidComputeSolution.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/FluidSolver.dir/fluidComputeSolution.cpp.s"
	/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /home/homspons/Documents/hiwi/elastictube1d/fluid-cpp/fluidComputeSolution.cpp -o CMakeFiles/FluidSolver.dir/fluidComputeSolution.cpp.s

# Object files for target FluidSolver
FluidSolver_OBJECTS = \
"CMakeFiles/FluidSolver.dir/fluid_solver.cpp.o" \
"CMakeFiles/FluidSolver.dir/fluid_nl.cpp.o" \
"CMakeFiles/FluidSolver.dir/fluidComputeSolution.cpp.o"

# External object files for target FluidSolver
FluidSolver_EXTERNAL_OBJECTS =

FluidSolver: CMakeFiles/FluidSolver.dir/fluid_solver.cpp.o
FluidSolver: CMakeFiles/FluidSolver.dir/fluid_nl.cpp.o
FluidSolver: CMakeFiles/FluidSolver.dir/fluidComputeSolution.cpp.o
FluidSolver: CMakeFiles/FluidSolver.dir/build.make
FluidSolver: /usr/lib/x86_64-linux-gnu/libprecice.so.2.1.0
FluidSolver: /usr/lib/x86_64-linux-gnu/liblapack.so
FluidSolver: /usr/lib/x86_64-linux-gnu/libblas.so
FluidSolver: /usr/lib/x86_64-linux-gnu/openmpi/lib/libmpi_cxx.so
FluidSolver: /usr/lib/x86_64-linux-gnu/openmpi/lib/libmpi.so
FluidSolver: CMakeFiles/FluidSolver.dir/link.txt
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --bold --progress-dir=/home/homspons/Documents/hiwi/elastictube1d/fluid-cpp/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_4) "Linking CXX executable FluidSolver"
	$(CMAKE_COMMAND) -E cmake_link_script CMakeFiles/FluidSolver.dir/link.txt --verbose=$(VERBOSE)

# Rule to build all files generated by this target.
CMakeFiles/FluidSolver.dir/build: FluidSolver

.PHONY : CMakeFiles/FluidSolver.dir/build

CMakeFiles/FluidSolver.dir/clean:
	$(CMAKE_COMMAND) -P CMakeFiles/FluidSolver.dir/cmake_clean.cmake
.PHONY : CMakeFiles/FluidSolver.dir/clean

CMakeFiles/FluidSolver.dir/depend:
	cd /home/homspons/Documents/hiwi/elastictube1d/fluid-cpp/build && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" /home/homspons/Documents/hiwi/elastictube1d/fluid-cpp /home/homspons/Documents/hiwi/elastictube1d/fluid-cpp /home/homspons/Documents/hiwi/elastictube1d/fluid-cpp/build /home/homspons/Documents/hiwi/elastictube1d/fluid-cpp/build /home/homspons/Documents/hiwi/elastictube1d/fluid-cpp/build/CMakeFiles/FluidSolver.dir/DependInfo.cmake --color=$(COLOR)
.PHONY : CMakeFiles/FluidSolver.dir/depend

