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

# Include any dependencies generated for this target.
include corgi/mpi4cpp/test/CMakeFiles/send_recv.dir/depend.make

# Include the progress variables for this target.
include corgi/mpi4cpp/test/CMakeFiles/send_recv.dir/progress.make

# Include the compile flags for this target's objects.
include corgi/mpi4cpp/test/CMakeFiles/send_recv.dir/flags.make

corgi/mpi4cpp/test/CMakeFiles/send_recv.dir/test_send_recv.c++.o: corgi/mpi4cpp/test/CMakeFiles/send_recv.dir/flags.make
corgi/mpi4cpp/test/CMakeFiles/send_recv.dir/test_send_recv.c++.o: ../corgi/mpi4cpp/test/test_send_recv.c++
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/krissmedt/Code/runko/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_1) "Building CXX object corgi/mpi4cpp/test/CMakeFiles/send_recv.dir/test_send_recv.c++.o"
	cd /home/krissmedt/Code/runko/build/corgi/mpi4cpp/test && /usr/bin/mpic++  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/send_recv.dir/test_send_recv.c++.o -c /home/krissmedt/Code/runko/corgi/mpi4cpp/test/test_send_recv.c++

corgi/mpi4cpp/test/CMakeFiles/send_recv.dir/test_send_recv.c++.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/send_recv.dir/test_send_recv.c++.i"
	cd /home/krissmedt/Code/runko/build/corgi/mpi4cpp/test && /usr/bin/mpic++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /home/krissmedt/Code/runko/corgi/mpi4cpp/test/test_send_recv.c++ > CMakeFiles/send_recv.dir/test_send_recv.c++.i

corgi/mpi4cpp/test/CMakeFiles/send_recv.dir/test_send_recv.c++.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/send_recv.dir/test_send_recv.c++.s"
	cd /home/krissmedt/Code/runko/build/corgi/mpi4cpp/test && /usr/bin/mpic++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /home/krissmedt/Code/runko/corgi/mpi4cpp/test/test_send_recv.c++ -o CMakeFiles/send_recv.dir/test_send_recv.c++.s

# Object files for target send_recv
send_recv_OBJECTS = \
"CMakeFiles/send_recv.dir/test_send_recv.c++.o"

# External object files for target send_recv
send_recv_EXTERNAL_OBJECTS =

corgi/mpi4cpp/test/send_recv: corgi/mpi4cpp/test/CMakeFiles/send_recv.dir/test_send_recv.c++.o
corgi/mpi4cpp/test/send_recv: corgi/mpi4cpp/test/CMakeFiles/send_recv.dir/build.make
corgi/mpi4cpp/test/send_recv: corgi/mpi4cpp/test/CMakeFiles/send_recv.dir/link.txt
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --bold --progress-dir=/home/krissmedt/Code/runko/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_2) "Linking CXX executable send_recv"
	cd /home/krissmedt/Code/runko/build/corgi/mpi4cpp/test && $(CMAKE_COMMAND) -E cmake_link_script CMakeFiles/send_recv.dir/link.txt --verbose=$(VERBOSE)

# Rule to build all files generated by this target.
corgi/mpi4cpp/test/CMakeFiles/send_recv.dir/build: corgi/mpi4cpp/test/send_recv

.PHONY : corgi/mpi4cpp/test/CMakeFiles/send_recv.dir/build

corgi/mpi4cpp/test/CMakeFiles/send_recv.dir/clean:
	cd /home/krissmedt/Code/runko/build/corgi/mpi4cpp/test && $(CMAKE_COMMAND) -P CMakeFiles/send_recv.dir/cmake_clean.cmake
.PHONY : corgi/mpi4cpp/test/CMakeFiles/send_recv.dir/clean

corgi/mpi4cpp/test/CMakeFiles/send_recv.dir/depend:
	cd /home/krissmedt/Code/runko/build && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" /home/krissmedt/Code/runko /home/krissmedt/Code/runko/corgi/mpi4cpp/test /home/krissmedt/Code/runko/build /home/krissmedt/Code/runko/build/corgi/mpi4cpp/test /home/krissmedt/Code/runko/build/corgi/mpi4cpp/test/CMakeFiles/send_recv.dir/DependInfo.cmake --color=$(COLOR)
.PHONY : corgi/mpi4cpp/test/CMakeFiles/send_recv.dir/depend

