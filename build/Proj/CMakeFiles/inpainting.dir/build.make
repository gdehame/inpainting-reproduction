# CMAKE generated file: DO NOT EDIT!
# Generated by "Unix Makefiles" Generator, CMake Version 3.22

# Delete rule output on recipe failure.
.DELETE_ON_ERROR:

#=============================================================================
# Special targets provided by cmake.

# Disable implicit rules so canonical targets will work.
.SUFFIXES:

# Disable VCS-based implicit rules.
% : %,v

# Disable VCS-based implicit rules.
% : RCS/%

# Disable VCS-based implicit rules.
% : RCS/%,v

# Disable VCS-based implicit rules.
% : SCCS/s.%

# Disable VCS-based implicit rules.
% : s.%

.SUFFIXES: .hpux_make_needs_suffix_list

# Command-line flag to silence nested $(MAKE).
$(VERBOSE)MAKESILENT = -s

#Suppress display of executed commands.
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
RM = /usr/bin/cmake -E rm -f

# Escaping for special characters.
EQUALS = =

# The top-level source directory on which CMake was run.
CMAKE_SOURCE_DIR = /home/jab/M1/S2/ProjCGDI

# The top-level build directory on which CMake was run.
CMAKE_BINARY_DIR = /home/jab/M1/S2/ProjCGDI/build

# Include any dependencies generated for this target.
include Proj/CMakeFiles/inpainting.dir/depend.make
# Include any dependencies generated by the compiler for this target.
include Proj/CMakeFiles/inpainting.dir/compiler_depend.make

# Include the progress variables for this target.
include Proj/CMakeFiles/inpainting.dir/progress.make

# Include the compile flags for this target's objects.
include Proj/CMakeFiles/inpainting.dir/flags.make

Proj/CMakeFiles/inpainting.dir/inpainting.cpp.o: Proj/CMakeFiles/inpainting.dir/flags.make
Proj/CMakeFiles/inpainting.dir/inpainting.cpp.o: ../Proj/inpainting.cpp
Proj/CMakeFiles/inpainting.dir/inpainting.cpp.o: Proj/CMakeFiles/inpainting.dir/compiler_depend.ts
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/jab/M1/S2/ProjCGDI/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_1) "Building CXX object Proj/CMakeFiles/inpainting.dir/inpainting.cpp.o"
	cd /home/jab/M1/S2/ProjCGDI/build/Proj && /usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -MD -MT Proj/CMakeFiles/inpainting.dir/inpainting.cpp.o -MF CMakeFiles/inpainting.dir/inpainting.cpp.o.d -o CMakeFiles/inpainting.dir/inpainting.cpp.o -c /home/jab/M1/S2/ProjCGDI/Proj/inpainting.cpp

Proj/CMakeFiles/inpainting.dir/inpainting.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/inpainting.dir/inpainting.cpp.i"
	cd /home/jab/M1/S2/ProjCGDI/build/Proj && /usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /home/jab/M1/S2/ProjCGDI/Proj/inpainting.cpp > CMakeFiles/inpainting.dir/inpainting.cpp.i

Proj/CMakeFiles/inpainting.dir/inpainting.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/inpainting.dir/inpainting.cpp.s"
	cd /home/jab/M1/S2/ProjCGDI/build/Proj && /usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /home/jab/M1/S2/ProjCGDI/Proj/inpainting.cpp -o CMakeFiles/inpainting.dir/inpainting.cpp.s

# Object files for target inpainting
inpainting_OBJECTS = \
"CMakeFiles/inpainting.dir/inpainting.cpp.o"

# External object files for target inpainting
inpainting_EXTERNAL_OBJECTS =

Proj/inpainting: Proj/CMakeFiles/inpainting.dir/inpainting.cpp.o
Proj/inpainting: Proj/CMakeFiles/inpainting.dir/build.make
Proj/inpainting: Proj/CMakeFiles/inpainting.dir/link.txt
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --bold --progress-dir=/home/jab/M1/S2/ProjCGDI/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_2) "Linking CXX executable inpainting"
	cd /home/jab/M1/S2/ProjCGDI/build/Proj && $(CMAKE_COMMAND) -E cmake_link_script CMakeFiles/inpainting.dir/link.txt --verbose=$(VERBOSE)

# Rule to build all files generated by this target.
Proj/CMakeFiles/inpainting.dir/build: Proj/inpainting
.PHONY : Proj/CMakeFiles/inpainting.dir/build

Proj/CMakeFiles/inpainting.dir/clean:
	cd /home/jab/M1/S2/ProjCGDI/build/Proj && $(CMAKE_COMMAND) -P CMakeFiles/inpainting.dir/cmake_clean.cmake
.PHONY : Proj/CMakeFiles/inpainting.dir/clean

Proj/CMakeFiles/inpainting.dir/depend:
	cd /home/jab/M1/S2/ProjCGDI/build && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" /home/jab/M1/S2/ProjCGDI /home/jab/M1/S2/ProjCGDI/Proj /home/jab/M1/S2/ProjCGDI/build /home/jab/M1/S2/ProjCGDI/build/Proj /home/jab/M1/S2/ProjCGDI/build/Proj/CMakeFiles/inpainting.dir/DependInfo.cmake --color=$(COLOR)
.PHONY : Proj/CMakeFiles/inpainting.dir/depend
