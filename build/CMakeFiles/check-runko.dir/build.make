# CMAKE generated file: DO NOT EDIT!
# Generated by "Unix Makefiles" Generator, CMake Version 3.13

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
CMAKE_SOURCE_DIR = /home/krissmedt/Code/runko

# The top-level build directory on which CMake was run.
CMAKE_BINARY_DIR = /home/krissmedt/Code/runko/build

# Utility rule file for check-runko.

# Include the progress variables for this target.
include CMakeFiles/check-runko.dir/progress.make

CMakeFiles/check-runko: ../corgi/lib/pycorgi.cpython-37m-x86_64-linux-gnu.so
CMakeFiles/check-runko: ../lib/pyrunko.cpython-37m-x86_64-linux-gnu.so
	cd /home/krissmedt/Code/runko/lib && /usr/bin/python3.7 -m unittest discover -s ../tests/ -v

check-runko: CMakeFiles/check-runko
check-runko: CMakeFiles/check-runko.dir/build.make

.PHONY : check-runko

# Rule to build all files generated by this target.
CMakeFiles/check-runko.dir/build: check-runko

.PHONY : CMakeFiles/check-runko.dir/build

CMakeFiles/check-runko.dir/clean:
	$(CMAKE_COMMAND) -P CMakeFiles/check-runko.dir/cmake_clean.cmake
.PHONY : CMakeFiles/check-runko.dir/clean

CMakeFiles/check-runko.dir/depend:
	cd /home/krissmedt/Code/runko/build && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" /home/krissmedt/Code/runko /home/krissmedt/Code/runko /home/krissmedt/Code/runko/build /home/krissmedt/Code/runko/build /home/krissmedt/Code/runko/build/CMakeFiles/check-runko.dir/DependInfo.cmake --color=$(COLOR)
.PHONY : CMakeFiles/check-runko.dir/depend

