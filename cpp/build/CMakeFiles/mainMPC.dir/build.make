# CMAKE generated file: DO NOT EDIT!
# Generated by "Unix Makefiles" Generator, CMake Version 2.8

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

# The program to use to edit the cache.
CMAKE_EDIT_COMMAND = /usr/bin/ccmake

# The top-level source directory on which CMake was run.
CMAKE_SOURCE_DIR = /home/gkharish/softdev/DDP/cpp

# The top-level build directory on which CMake was run.
CMAKE_BINARY_DIR = /home/gkharish/softdev/DDP/cpp/build

# Include any dependencies generated for this target.
include CMakeFiles/mainMPC.dir/depend.make

# Include the progress variables for this target.
include CMakeFiles/mainMPC.dir/progress.make

# Include the compile flags for this target's objects.
include CMakeFiles/mainMPC.dir/flags.make

CMakeFiles/mainMPC.dir/test/mainMPC.cpp.o: CMakeFiles/mainMPC.dir/flags.make
CMakeFiles/mainMPC.dir/test/mainMPC.cpp.o: ../test/mainMPC.cpp
	$(CMAKE_COMMAND) -E cmake_progress_report /home/gkharish/softdev/DDP/cpp/build/CMakeFiles $(CMAKE_PROGRESS_1)
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Building CXX object CMakeFiles/mainMPC.dir/test/mainMPC.cpp.o"
	/usr/bin/c++   $(CXX_DEFINES) $(CXX_FLAGS) -o CMakeFiles/mainMPC.dir/test/mainMPC.cpp.o -c /home/gkharish/softdev/DDP/cpp/test/mainMPC.cpp

CMakeFiles/mainMPC.dir/test/mainMPC.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/mainMPC.dir/test/mainMPC.cpp.i"
	/usr/bin/c++  $(CXX_DEFINES) $(CXX_FLAGS) -E /home/gkharish/softdev/DDP/cpp/test/mainMPC.cpp > CMakeFiles/mainMPC.dir/test/mainMPC.cpp.i

CMakeFiles/mainMPC.dir/test/mainMPC.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/mainMPC.dir/test/mainMPC.cpp.s"
	/usr/bin/c++  $(CXX_DEFINES) $(CXX_FLAGS) -S /home/gkharish/softdev/DDP/cpp/test/mainMPC.cpp -o CMakeFiles/mainMPC.dir/test/mainMPC.cpp.s

CMakeFiles/mainMPC.dir/test/mainMPC.cpp.o.requires:
.PHONY : CMakeFiles/mainMPC.dir/test/mainMPC.cpp.o.requires

CMakeFiles/mainMPC.dir/test/mainMPC.cpp.o.provides: CMakeFiles/mainMPC.dir/test/mainMPC.cpp.o.requires
	$(MAKE) -f CMakeFiles/mainMPC.dir/build.make CMakeFiles/mainMPC.dir/test/mainMPC.cpp.o.provides.build
.PHONY : CMakeFiles/mainMPC.dir/test/mainMPC.cpp.o.provides

CMakeFiles/mainMPC.dir/test/mainMPC.cpp.o.provides.build: CMakeFiles/mainMPC.dir/test/mainMPC.cpp.o

CMakeFiles/mainMPC.dir/src/ilqrsolver.cpp.o: CMakeFiles/mainMPC.dir/flags.make
CMakeFiles/mainMPC.dir/src/ilqrsolver.cpp.o: ../src/ilqrsolver.cpp
	$(CMAKE_COMMAND) -E cmake_progress_report /home/gkharish/softdev/DDP/cpp/build/CMakeFiles $(CMAKE_PROGRESS_2)
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Building CXX object CMakeFiles/mainMPC.dir/src/ilqrsolver.cpp.o"
	/usr/bin/c++   $(CXX_DEFINES) $(CXX_FLAGS) -o CMakeFiles/mainMPC.dir/src/ilqrsolver.cpp.o -c /home/gkharish/softdev/DDP/cpp/src/ilqrsolver.cpp

CMakeFiles/mainMPC.dir/src/ilqrsolver.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/mainMPC.dir/src/ilqrsolver.cpp.i"
	/usr/bin/c++  $(CXX_DEFINES) $(CXX_FLAGS) -E /home/gkharish/softdev/DDP/cpp/src/ilqrsolver.cpp > CMakeFiles/mainMPC.dir/src/ilqrsolver.cpp.i

CMakeFiles/mainMPC.dir/src/ilqrsolver.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/mainMPC.dir/src/ilqrsolver.cpp.s"
	/usr/bin/c++  $(CXX_DEFINES) $(CXX_FLAGS) -S /home/gkharish/softdev/DDP/cpp/src/ilqrsolver.cpp -o CMakeFiles/mainMPC.dir/src/ilqrsolver.cpp.s

CMakeFiles/mainMPC.dir/src/ilqrsolver.cpp.o.requires:
.PHONY : CMakeFiles/mainMPC.dir/src/ilqrsolver.cpp.o.requires

CMakeFiles/mainMPC.dir/src/ilqrsolver.cpp.o.provides: CMakeFiles/mainMPC.dir/src/ilqrsolver.cpp.o.requires
	$(MAKE) -f CMakeFiles/mainMPC.dir/build.make CMakeFiles/mainMPC.dir/src/ilqrsolver.cpp.o.provides.build
.PHONY : CMakeFiles/mainMPC.dir/src/ilqrsolver.cpp.o.provides

CMakeFiles/mainMPC.dir/src/ilqrsolver.cpp.o.provides.build: CMakeFiles/mainMPC.dir/src/ilqrsolver.cpp.o

CMakeFiles/mainMPC.dir/src/pneumaticarmnonlinearmodel.cpp.o: CMakeFiles/mainMPC.dir/flags.make
CMakeFiles/mainMPC.dir/src/pneumaticarmnonlinearmodel.cpp.o: ../src/pneumaticarmnonlinearmodel.cpp
	$(CMAKE_COMMAND) -E cmake_progress_report /home/gkharish/softdev/DDP/cpp/build/CMakeFiles $(CMAKE_PROGRESS_3)
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Building CXX object CMakeFiles/mainMPC.dir/src/pneumaticarmnonlinearmodel.cpp.o"
	/usr/bin/c++   $(CXX_DEFINES) $(CXX_FLAGS) -o CMakeFiles/mainMPC.dir/src/pneumaticarmnonlinearmodel.cpp.o -c /home/gkharish/softdev/DDP/cpp/src/pneumaticarmnonlinearmodel.cpp

CMakeFiles/mainMPC.dir/src/pneumaticarmnonlinearmodel.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/mainMPC.dir/src/pneumaticarmnonlinearmodel.cpp.i"
	/usr/bin/c++  $(CXX_DEFINES) $(CXX_FLAGS) -E /home/gkharish/softdev/DDP/cpp/src/pneumaticarmnonlinearmodel.cpp > CMakeFiles/mainMPC.dir/src/pneumaticarmnonlinearmodel.cpp.i

CMakeFiles/mainMPC.dir/src/pneumaticarmnonlinearmodel.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/mainMPC.dir/src/pneumaticarmnonlinearmodel.cpp.s"
	/usr/bin/c++  $(CXX_DEFINES) $(CXX_FLAGS) -S /home/gkharish/softdev/DDP/cpp/src/pneumaticarmnonlinearmodel.cpp -o CMakeFiles/mainMPC.dir/src/pneumaticarmnonlinearmodel.cpp.s

CMakeFiles/mainMPC.dir/src/pneumaticarmnonlinearmodel.cpp.o.requires:
.PHONY : CMakeFiles/mainMPC.dir/src/pneumaticarmnonlinearmodel.cpp.o.requires

CMakeFiles/mainMPC.dir/src/pneumaticarmnonlinearmodel.cpp.o.provides: CMakeFiles/mainMPC.dir/src/pneumaticarmnonlinearmodel.cpp.o.requires
	$(MAKE) -f CMakeFiles/mainMPC.dir/build.make CMakeFiles/mainMPC.dir/src/pneumaticarmnonlinearmodel.cpp.o.provides.build
.PHONY : CMakeFiles/mainMPC.dir/src/pneumaticarmnonlinearmodel.cpp.o.provides

CMakeFiles/mainMPC.dir/src/pneumaticarmnonlinearmodel.cpp.o.provides.build: CMakeFiles/mainMPC.dir/src/pneumaticarmnonlinearmodel.cpp.o

CMakeFiles/mainMPC.dir/src/pneumaticarmelbow2ndordermodel.cpp.o: CMakeFiles/mainMPC.dir/flags.make
CMakeFiles/mainMPC.dir/src/pneumaticarmelbow2ndordermodel.cpp.o: ../src/pneumaticarmelbow2ndordermodel.cpp
	$(CMAKE_COMMAND) -E cmake_progress_report /home/gkharish/softdev/DDP/cpp/build/CMakeFiles $(CMAKE_PROGRESS_4)
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Building CXX object CMakeFiles/mainMPC.dir/src/pneumaticarmelbow2ndordermodel.cpp.o"
	/usr/bin/c++   $(CXX_DEFINES) $(CXX_FLAGS) -o CMakeFiles/mainMPC.dir/src/pneumaticarmelbow2ndordermodel.cpp.o -c /home/gkharish/softdev/DDP/cpp/src/pneumaticarmelbow2ndordermodel.cpp

CMakeFiles/mainMPC.dir/src/pneumaticarmelbow2ndordermodel.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/mainMPC.dir/src/pneumaticarmelbow2ndordermodel.cpp.i"
	/usr/bin/c++  $(CXX_DEFINES) $(CXX_FLAGS) -E /home/gkharish/softdev/DDP/cpp/src/pneumaticarmelbow2ndordermodel.cpp > CMakeFiles/mainMPC.dir/src/pneumaticarmelbow2ndordermodel.cpp.i

CMakeFiles/mainMPC.dir/src/pneumaticarmelbow2ndordermodel.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/mainMPC.dir/src/pneumaticarmelbow2ndordermodel.cpp.s"
	/usr/bin/c++  $(CXX_DEFINES) $(CXX_FLAGS) -S /home/gkharish/softdev/DDP/cpp/src/pneumaticarmelbow2ndordermodel.cpp -o CMakeFiles/mainMPC.dir/src/pneumaticarmelbow2ndordermodel.cpp.s

CMakeFiles/mainMPC.dir/src/pneumaticarmelbow2ndordermodel.cpp.o.requires:
.PHONY : CMakeFiles/mainMPC.dir/src/pneumaticarmelbow2ndordermodel.cpp.o.requires

CMakeFiles/mainMPC.dir/src/pneumaticarmelbow2ndordermodel.cpp.o.provides: CMakeFiles/mainMPC.dir/src/pneumaticarmelbow2ndordermodel.cpp.o.requires
	$(MAKE) -f CMakeFiles/mainMPC.dir/build.make CMakeFiles/mainMPC.dir/src/pneumaticarmelbow2ndordermodel.cpp.o.provides.build
.PHONY : CMakeFiles/mainMPC.dir/src/pneumaticarmelbow2ndordermodel.cpp.o.provides

CMakeFiles/mainMPC.dir/src/pneumaticarmelbow2ndordermodel.cpp.o.provides.build: CMakeFiles/mainMPC.dir/src/pneumaticarmelbow2ndordermodel.cpp.o

CMakeFiles/mainMPC.dir/src/costfunctionpneumaticarmelbow.cpp.o: CMakeFiles/mainMPC.dir/flags.make
CMakeFiles/mainMPC.dir/src/costfunctionpneumaticarmelbow.cpp.o: ../src/costfunctionpneumaticarmelbow.cpp
	$(CMAKE_COMMAND) -E cmake_progress_report /home/gkharish/softdev/DDP/cpp/build/CMakeFiles $(CMAKE_PROGRESS_5)
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Building CXX object CMakeFiles/mainMPC.dir/src/costfunctionpneumaticarmelbow.cpp.o"
	/usr/bin/c++   $(CXX_DEFINES) $(CXX_FLAGS) -o CMakeFiles/mainMPC.dir/src/costfunctionpneumaticarmelbow.cpp.o -c /home/gkharish/softdev/DDP/cpp/src/costfunctionpneumaticarmelbow.cpp

CMakeFiles/mainMPC.dir/src/costfunctionpneumaticarmelbow.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/mainMPC.dir/src/costfunctionpneumaticarmelbow.cpp.i"
	/usr/bin/c++  $(CXX_DEFINES) $(CXX_FLAGS) -E /home/gkharish/softdev/DDP/cpp/src/costfunctionpneumaticarmelbow.cpp > CMakeFiles/mainMPC.dir/src/costfunctionpneumaticarmelbow.cpp.i

CMakeFiles/mainMPC.dir/src/costfunctionpneumaticarmelbow.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/mainMPC.dir/src/costfunctionpneumaticarmelbow.cpp.s"
	/usr/bin/c++  $(CXX_DEFINES) $(CXX_FLAGS) -S /home/gkharish/softdev/DDP/cpp/src/costfunctionpneumaticarmelbow.cpp -o CMakeFiles/mainMPC.dir/src/costfunctionpneumaticarmelbow.cpp.s

CMakeFiles/mainMPC.dir/src/costfunctionpneumaticarmelbow.cpp.o.requires:
.PHONY : CMakeFiles/mainMPC.dir/src/costfunctionpneumaticarmelbow.cpp.o.requires

CMakeFiles/mainMPC.dir/src/costfunctionpneumaticarmelbow.cpp.o.provides: CMakeFiles/mainMPC.dir/src/costfunctionpneumaticarmelbow.cpp.o.requires
	$(MAKE) -f CMakeFiles/mainMPC.dir/build.make CMakeFiles/mainMPC.dir/src/costfunctionpneumaticarmelbow.cpp.o.provides.build
.PHONY : CMakeFiles/mainMPC.dir/src/costfunctionpneumaticarmelbow.cpp.o.provides

CMakeFiles/mainMPC.dir/src/costfunctionpneumaticarmelbow.cpp.o.provides.build: CMakeFiles/mainMPC.dir/src/costfunctionpneumaticarmelbow.cpp.o

CMakeFiles/mainMPC.dir/src/costfunctionromeoactuator.cpp.o: CMakeFiles/mainMPC.dir/flags.make
CMakeFiles/mainMPC.dir/src/costfunctionromeoactuator.cpp.o: ../src/costfunctionromeoactuator.cpp
	$(CMAKE_COMMAND) -E cmake_progress_report /home/gkharish/softdev/DDP/cpp/build/CMakeFiles $(CMAKE_PROGRESS_6)
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Building CXX object CMakeFiles/mainMPC.dir/src/costfunctionromeoactuator.cpp.o"
	/usr/bin/c++   $(CXX_DEFINES) $(CXX_FLAGS) -o CMakeFiles/mainMPC.dir/src/costfunctionromeoactuator.cpp.o -c /home/gkharish/softdev/DDP/cpp/src/costfunctionromeoactuator.cpp

CMakeFiles/mainMPC.dir/src/costfunctionromeoactuator.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/mainMPC.dir/src/costfunctionromeoactuator.cpp.i"
	/usr/bin/c++  $(CXX_DEFINES) $(CXX_FLAGS) -E /home/gkharish/softdev/DDP/cpp/src/costfunctionromeoactuator.cpp > CMakeFiles/mainMPC.dir/src/costfunctionromeoactuator.cpp.i

CMakeFiles/mainMPC.dir/src/costfunctionromeoactuator.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/mainMPC.dir/src/costfunctionromeoactuator.cpp.s"
	/usr/bin/c++  $(CXX_DEFINES) $(CXX_FLAGS) -S /home/gkharish/softdev/DDP/cpp/src/costfunctionromeoactuator.cpp -o CMakeFiles/mainMPC.dir/src/costfunctionromeoactuator.cpp.s

CMakeFiles/mainMPC.dir/src/costfunctionromeoactuator.cpp.o.requires:
.PHONY : CMakeFiles/mainMPC.dir/src/costfunctionromeoactuator.cpp.o.requires

CMakeFiles/mainMPC.dir/src/costfunctionromeoactuator.cpp.o.provides: CMakeFiles/mainMPC.dir/src/costfunctionromeoactuator.cpp.o.requires
	$(MAKE) -f CMakeFiles/mainMPC.dir/build.make CMakeFiles/mainMPC.dir/src/costfunctionromeoactuator.cpp.o.provides.build
.PHONY : CMakeFiles/mainMPC.dir/src/costfunctionromeoactuator.cpp.o.provides

CMakeFiles/mainMPC.dir/src/costfunctionromeoactuator.cpp.o.provides.build: CMakeFiles/mainMPC.dir/src/costfunctionromeoactuator.cpp.o

CMakeFiles/mainMPC.dir/src/costfunction.cpp.o: CMakeFiles/mainMPC.dir/flags.make
CMakeFiles/mainMPC.dir/src/costfunction.cpp.o: ../src/costfunction.cpp
	$(CMAKE_COMMAND) -E cmake_progress_report /home/gkharish/softdev/DDP/cpp/build/CMakeFiles $(CMAKE_PROGRESS_7)
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Building CXX object CMakeFiles/mainMPC.dir/src/costfunction.cpp.o"
	/usr/bin/c++   $(CXX_DEFINES) $(CXX_FLAGS) -o CMakeFiles/mainMPC.dir/src/costfunction.cpp.o -c /home/gkharish/softdev/DDP/cpp/src/costfunction.cpp

CMakeFiles/mainMPC.dir/src/costfunction.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/mainMPC.dir/src/costfunction.cpp.i"
	/usr/bin/c++  $(CXX_DEFINES) $(CXX_FLAGS) -E /home/gkharish/softdev/DDP/cpp/src/costfunction.cpp > CMakeFiles/mainMPC.dir/src/costfunction.cpp.i

CMakeFiles/mainMPC.dir/src/costfunction.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/mainMPC.dir/src/costfunction.cpp.s"
	/usr/bin/c++  $(CXX_DEFINES) $(CXX_FLAGS) -S /home/gkharish/softdev/DDP/cpp/src/costfunction.cpp -o CMakeFiles/mainMPC.dir/src/costfunction.cpp.s

CMakeFiles/mainMPC.dir/src/costfunction.cpp.o.requires:
.PHONY : CMakeFiles/mainMPC.dir/src/costfunction.cpp.o.requires

CMakeFiles/mainMPC.dir/src/costfunction.cpp.o.provides: CMakeFiles/mainMPC.dir/src/costfunction.cpp.o.requires
	$(MAKE) -f CMakeFiles/mainMPC.dir/build.make CMakeFiles/mainMPC.dir/src/costfunction.cpp.o.provides.build
.PHONY : CMakeFiles/mainMPC.dir/src/costfunction.cpp.o.provides

CMakeFiles/mainMPC.dir/src/costfunction.cpp.o.provides.build: CMakeFiles/mainMPC.dir/src/costfunction.cpp.o

CMakeFiles/mainMPC.dir/src/dynamicmodel.cpp.o: CMakeFiles/mainMPC.dir/flags.make
CMakeFiles/mainMPC.dir/src/dynamicmodel.cpp.o: ../src/dynamicmodel.cpp
	$(CMAKE_COMMAND) -E cmake_progress_report /home/gkharish/softdev/DDP/cpp/build/CMakeFiles $(CMAKE_PROGRESS_8)
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Building CXX object CMakeFiles/mainMPC.dir/src/dynamicmodel.cpp.o"
	/usr/bin/c++   $(CXX_DEFINES) $(CXX_FLAGS) -o CMakeFiles/mainMPC.dir/src/dynamicmodel.cpp.o -c /home/gkharish/softdev/DDP/cpp/src/dynamicmodel.cpp

CMakeFiles/mainMPC.dir/src/dynamicmodel.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/mainMPC.dir/src/dynamicmodel.cpp.i"
	/usr/bin/c++  $(CXX_DEFINES) $(CXX_FLAGS) -E /home/gkharish/softdev/DDP/cpp/src/dynamicmodel.cpp > CMakeFiles/mainMPC.dir/src/dynamicmodel.cpp.i

CMakeFiles/mainMPC.dir/src/dynamicmodel.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/mainMPC.dir/src/dynamicmodel.cpp.s"
	/usr/bin/c++  $(CXX_DEFINES) $(CXX_FLAGS) -S /home/gkharish/softdev/DDP/cpp/src/dynamicmodel.cpp -o CMakeFiles/mainMPC.dir/src/dynamicmodel.cpp.s

CMakeFiles/mainMPC.dir/src/dynamicmodel.cpp.o.requires:
.PHONY : CMakeFiles/mainMPC.dir/src/dynamicmodel.cpp.o.requires

CMakeFiles/mainMPC.dir/src/dynamicmodel.cpp.o.provides: CMakeFiles/mainMPC.dir/src/dynamicmodel.cpp.o.requires
	$(MAKE) -f CMakeFiles/mainMPC.dir/build.make CMakeFiles/mainMPC.dir/src/dynamicmodel.cpp.o.provides.build
.PHONY : CMakeFiles/mainMPC.dir/src/dynamicmodel.cpp.o.provides

CMakeFiles/mainMPC.dir/src/dynamicmodel.cpp.o.provides.build: CMakeFiles/mainMPC.dir/src/dynamicmodel.cpp.o

CMakeFiles/mainMPC.dir/src/romeosimpleactuator.cpp.o: CMakeFiles/mainMPC.dir/flags.make
CMakeFiles/mainMPC.dir/src/romeosimpleactuator.cpp.o: ../src/romeosimpleactuator.cpp
	$(CMAKE_COMMAND) -E cmake_progress_report /home/gkharish/softdev/DDP/cpp/build/CMakeFiles $(CMAKE_PROGRESS_9)
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Building CXX object CMakeFiles/mainMPC.dir/src/romeosimpleactuator.cpp.o"
	/usr/bin/c++   $(CXX_DEFINES) $(CXX_FLAGS) -o CMakeFiles/mainMPC.dir/src/romeosimpleactuator.cpp.o -c /home/gkharish/softdev/DDP/cpp/src/romeosimpleactuator.cpp

CMakeFiles/mainMPC.dir/src/romeosimpleactuator.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/mainMPC.dir/src/romeosimpleactuator.cpp.i"
	/usr/bin/c++  $(CXX_DEFINES) $(CXX_FLAGS) -E /home/gkharish/softdev/DDP/cpp/src/romeosimpleactuator.cpp > CMakeFiles/mainMPC.dir/src/romeosimpleactuator.cpp.i

CMakeFiles/mainMPC.dir/src/romeosimpleactuator.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/mainMPC.dir/src/romeosimpleactuator.cpp.s"
	/usr/bin/c++  $(CXX_DEFINES) $(CXX_FLAGS) -S /home/gkharish/softdev/DDP/cpp/src/romeosimpleactuator.cpp -o CMakeFiles/mainMPC.dir/src/romeosimpleactuator.cpp.s

CMakeFiles/mainMPC.dir/src/romeosimpleactuator.cpp.o.requires:
.PHONY : CMakeFiles/mainMPC.dir/src/romeosimpleactuator.cpp.o.requires

CMakeFiles/mainMPC.dir/src/romeosimpleactuator.cpp.o.provides: CMakeFiles/mainMPC.dir/src/romeosimpleactuator.cpp.o.requires
	$(MAKE) -f CMakeFiles/mainMPC.dir/build.make CMakeFiles/mainMPC.dir/src/romeosimpleactuator.cpp.o.provides.build
.PHONY : CMakeFiles/mainMPC.dir/src/romeosimpleactuator.cpp.o.provides

CMakeFiles/mainMPC.dir/src/romeosimpleactuator.cpp.o.provides.build: CMakeFiles/mainMPC.dir/src/romeosimpleactuator.cpp.o

CMakeFiles/mainMPC.dir/src/romeolinearactuator.cpp.o: CMakeFiles/mainMPC.dir/flags.make
CMakeFiles/mainMPC.dir/src/romeolinearactuator.cpp.o: ../src/romeolinearactuator.cpp
	$(CMAKE_COMMAND) -E cmake_progress_report /home/gkharish/softdev/DDP/cpp/build/CMakeFiles $(CMAKE_PROGRESS_10)
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Building CXX object CMakeFiles/mainMPC.dir/src/romeolinearactuator.cpp.o"
	/usr/bin/c++   $(CXX_DEFINES) $(CXX_FLAGS) -o CMakeFiles/mainMPC.dir/src/romeolinearactuator.cpp.o -c /home/gkharish/softdev/DDP/cpp/src/romeolinearactuator.cpp

CMakeFiles/mainMPC.dir/src/romeolinearactuator.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/mainMPC.dir/src/romeolinearactuator.cpp.i"
	/usr/bin/c++  $(CXX_DEFINES) $(CXX_FLAGS) -E /home/gkharish/softdev/DDP/cpp/src/romeolinearactuator.cpp > CMakeFiles/mainMPC.dir/src/romeolinearactuator.cpp.i

CMakeFiles/mainMPC.dir/src/romeolinearactuator.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/mainMPC.dir/src/romeolinearactuator.cpp.s"
	/usr/bin/c++  $(CXX_DEFINES) $(CXX_FLAGS) -S /home/gkharish/softdev/DDP/cpp/src/romeolinearactuator.cpp -o CMakeFiles/mainMPC.dir/src/romeolinearactuator.cpp.s

CMakeFiles/mainMPC.dir/src/romeolinearactuator.cpp.o.requires:
.PHONY : CMakeFiles/mainMPC.dir/src/romeolinearactuator.cpp.o.requires

CMakeFiles/mainMPC.dir/src/romeolinearactuator.cpp.o.provides: CMakeFiles/mainMPC.dir/src/romeolinearactuator.cpp.o.requires
	$(MAKE) -f CMakeFiles/mainMPC.dir/build.make CMakeFiles/mainMPC.dir/src/romeolinearactuator.cpp.o.provides.build
.PHONY : CMakeFiles/mainMPC.dir/src/romeolinearactuator.cpp.o.provides

CMakeFiles/mainMPC.dir/src/romeolinearactuator.cpp.o.provides.build: CMakeFiles/mainMPC.dir/src/romeolinearactuator.cpp.o

CMakeFiles/mainMPC.dir/src/pneumaticarmelbowpiecelinear.cpp.o: CMakeFiles/mainMPC.dir/flags.make
CMakeFiles/mainMPC.dir/src/pneumaticarmelbowpiecelinear.cpp.o: ../src/pneumaticarmelbowpiecelinear.cpp
	$(CMAKE_COMMAND) -E cmake_progress_report /home/gkharish/softdev/DDP/cpp/build/CMakeFiles $(CMAKE_PROGRESS_11)
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Building CXX object CMakeFiles/mainMPC.dir/src/pneumaticarmelbowpiecelinear.cpp.o"
	/usr/bin/c++   $(CXX_DEFINES) $(CXX_FLAGS) -o CMakeFiles/mainMPC.dir/src/pneumaticarmelbowpiecelinear.cpp.o -c /home/gkharish/softdev/DDP/cpp/src/pneumaticarmelbowpiecelinear.cpp

CMakeFiles/mainMPC.dir/src/pneumaticarmelbowpiecelinear.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/mainMPC.dir/src/pneumaticarmelbowpiecelinear.cpp.i"
	/usr/bin/c++  $(CXX_DEFINES) $(CXX_FLAGS) -E /home/gkharish/softdev/DDP/cpp/src/pneumaticarmelbowpiecelinear.cpp > CMakeFiles/mainMPC.dir/src/pneumaticarmelbowpiecelinear.cpp.i

CMakeFiles/mainMPC.dir/src/pneumaticarmelbowpiecelinear.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/mainMPC.dir/src/pneumaticarmelbowpiecelinear.cpp.s"
	/usr/bin/c++  $(CXX_DEFINES) $(CXX_FLAGS) -S /home/gkharish/softdev/DDP/cpp/src/pneumaticarmelbowpiecelinear.cpp -o CMakeFiles/mainMPC.dir/src/pneumaticarmelbowpiecelinear.cpp.s

CMakeFiles/mainMPC.dir/src/pneumaticarmelbowpiecelinear.cpp.o.requires:
.PHONY : CMakeFiles/mainMPC.dir/src/pneumaticarmelbowpiecelinear.cpp.o.requires

CMakeFiles/mainMPC.dir/src/pneumaticarmelbowpiecelinear.cpp.o.provides: CMakeFiles/mainMPC.dir/src/pneumaticarmelbowpiecelinear.cpp.o.requires
	$(MAKE) -f CMakeFiles/mainMPC.dir/build.make CMakeFiles/mainMPC.dir/src/pneumaticarmelbowpiecelinear.cpp.o.provides.build
.PHONY : CMakeFiles/mainMPC.dir/src/pneumaticarmelbowpiecelinear.cpp.o.provides

CMakeFiles/mainMPC.dir/src/pneumaticarmelbowpiecelinear.cpp.o.provides.build: CMakeFiles/mainMPC.dir/src/pneumaticarmelbowpiecelinear.cpp.o

CMakeFiles/mainMPC.dir/src/pneumaticarmelbowlinear.cpp.o: CMakeFiles/mainMPC.dir/flags.make
CMakeFiles/mainMPC.dir/src/pneumaticarmelbowlinear.cpp.o: ../src/pneumaticarmelbowlinear.cpp
	$(CMAKE_COMMAND) -E cmake_progress_report /home/gkharish/softdev/DDP/cpp/build/CMakeFiles $(CMAKE_PROGRESS_12)
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Building CXX object CMakeFiles/mainMPC.dir/src/pneumaticarmelbowlinear.cpp.o"
	/usr/bin/c++   $(CXX_DEFINES) $(CXX_FLAGS) -o CMakeFiles/mainMPC.dir/src/pneumaticarmelbowlinear.cpp.o -c /home/gkharish/softdev/DDP/cpp/src/pneumaticarmelbowlinear.cpp

CMakeFiles/mainMPC.dir/src/pneumaticarmelbowlinear.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/mainMPC.dir/src/pneumaticarmelbowlinear.cpp.i"
	/usr/bin/c++  $(CXX_DEFINES) $(CXX_FLAGS) -E /home/gkharish/softdev/DDP/cpp/src/pneumaticarmelbowlinear.cpp > CMakeFiles/mainMPC.dir/src/pneumaticarmelbowlinear.cpp.i

CMakeFiles/mainMPC.dir/src/pneumaticarmelbowlinear.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/mainMPC.dir/src/pneumaticarmelbowlinear.cpp.s"
	/usr/bin/c++  $(CXX_DEFINES) $(CXX_FLAGS) -S /home/gkharish/softdev/DDP/cpp/src/pneumaticarmelbowlinear.cpp -o CMakeFiles/mainMPC.dir/src/pneumaticarmelbowlinear.cpp.s

CMakeFiles/mainMPC.dir/src/pneumaticarmelbowlinear.cpp.o.requires:
.PHONY : CMakeFiles/mainMPC.dir/src/pneumaticarmelbowlinear.cpp.o.requires

CMakeFiles/mainMPC.dir/src/pneumaticarmelbowlinear.cpp.o.provides: CMakeFiles/mainMPC.dir/src/pneumaticarmelbowlinear.cpp.o.requires
	$(MAKE) -f CMakeFiles/mainMPC.dir/build.make CMakeFiles/mainMPC.dir/src/pneumaticarmelbowlinear.cpp.o.provides.build
.PHONY : CMakeFiles/mainMPC.dir/src/pneumaticarmelbowlinear.cpp.o.provides

CMakeFiles/mainMPC.dir/src/pneumaticarmelbowlinear.cpp.o.provides.build: CMakeFiles/mainMPC.dir/src/pneumaticarmelbowlinear.cpp.o

# Object files for target mainMPC
mainMPC_OBJECTS = \
"CMakeFiles/mainMPC.dir/test/mainMPC.cpp.o" \
"CMakeFiles/mainMPC.dir/src/ilqrsolver.cpp.o" \
"CMakeFiles/mainMPC.dir/src/pneumaticarmnonlinearmodel.cpp.o" \
"CMakeFiles/mainMPC.dir/src/pneumaticarmelbow2ndordermodel.cpp.o" \
"CMakeFiles/mainMPC.dir/src/costfunctionpneumaticarmelbow.cpp.o" \
"CMakeFiles/mainMPC.dir/src/costfunctionromeoactuator.cpp.o" \
"CMakeFiles/mainMPC.dir/src/costfunction.cpp.o" \
"CMakeFiles/mainMPC.dir/src/dynamicmodel.cpp.o" \
"CMakeFiles/mainMPC.dir/src/romeosimpleactuator.cpp.o" \
"CMakeFiles/mainMPC.dir/src/romeolinearactuator.cpp.o" \
"CMakeFiles/mainMPC.dir/src/pneumaticarmelbowpiecelinear.cpp.o" \
"CMakeFiles/mainMPC.dir/src/pneumaticarmelbowlinear.cpp.o"

# External object files for target mainMPC
mainMPC_EXTERNAL_OBJECTS =

mainMPC: CMakeFiles/mainMPC.dir/test/mainMPC.cpp.o
mainMPC: CMakeFiles/mainMPC.dir/src/ilqrsolver.cpp.o
mainMPC: CMakeFiles/mainMPC.dir/src/pneumaticarmnonlinearmodel.cpp.o
mainMPC: CMakeFiles/mainMPC.dir/src/pneumaticarmelbow2ndordermodel.cpp.o
mainMPC: CMakeFiles/mainMPC.dir/src/costfunctionpneumaticarmelbow.cpp.o
mainMPC: CMakeFiles/mainMPC.dir/src/costfunctionromeoactuator.cpp.o
mainMPC: CMakeFiles/mainMPC.dir/src/costfunction.cpp.o
mainMPC: CMakeFiles/mainMPC.dir/src/dynamicmodel.cpp.o
mainMPC: CMakeFiles/mainMPC.dir/src/romeosimpleactuator.cpp.o
mainMPC: CMakeFiles/mainMPC.dir/src/romeolinearactuator.cpp.o
mainMPC: CMakeFiles/mainMPC.dir/src/pneumaticarmelbowpiecelinear.cpp.o
mainMPC: CMakeFiles/mainMPC.dir/src/pneumaticarmelbowlinear.cpp.o
mainMPC: CMakeFiles/mainMPC.dir/build.make
mainMPC: CMakeFiles/mainMPC.dir/link.txt
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --red --bold "Linking CXX executable mainMPC"
	$(CMAKE_COMMAND) -E cmake_link_script CMakeFiles/mainMPC.dir/link.txt --verbose=$(VERBOSE)

# Rule to build all files generated by this target.
CMakeFiles/mainMPC.dir/build: mainMPC
.PHONY : CMakeFiles/mainMPC.dir/build

CMakeFiles/mainMPC.dir/requires: CMakeFiles/mainMPC.dir/test/mainMPC.cpp.o.requires
CMakeFiles/mainMPC.dir/requires: CMakeFiles/mainMPC.dir/src/ilqrsolver.cpp.o.requires
CMakeFiles/mainMPC.dir/requires: CMakeFiles/mainMPC.dir/src/pneumaticarmnonlinearmodel.cpp.o.requires
CMakeFiles/mainMPC.dir/requires: CMakeFiles/mainMPC.dir/src/pneumaticarmelbow2ndordermodel.cpp.o.requires
CMakeFiles/mainMPC.dir/requires: CMakeFiles/mainMPC.dir/src/costfunctionpneumaticarmelbow.cpp.o.requires
CMakeFiles/mainMPC.dir/requires: CMakeFiles/mainMPC.dir/src/costfunctionromeoactuator.cpp.o.requires
CMakeFiles/mainMPC.dir/requires: CMakeFiles/mainMPC.dir/src/costfunction.cpp.o.requires
CMakeFiles/mainMPC.dir/requires: CMakeFiles/mainMPC.dir/src/dynamicmodel.cpp.o.requires
CMakeFiles/mainMPC.dir/requires: CMakeFiles/mainMPC.dir/src/romeosimpleactuator.cpp.o.requires
CMakeFiles/mainMPC.dir/requires: CMakeFiles/mainMPC.dir/src/romeolinearactuator.cpp.o.requires
CMakeFiles/mainMPC.dir/requires: CMakeFiles/mainMPC.dir/src/pneumaticarmelbowpiecelinear.cpp.o.requires
CMakeFiles/mainMPC.dir/requires: CMakeFiles/mainMPC.dir/src/pneumaticarmelbowlinear.cpp.o.requires
.PHONY : CMakeFiles/mainMPC.dir/requires

CMakeFiles/mainMPC.dir/clean:
	$(CMAKE_COMMAND) -P CMakeFiles/mainMPC.dir/cmake_clean.cmake
.PHONY : CMakeFiles/mainMPC.dir/clean

CMakeFiles/mainMPC.dir/depend:
	cd /home/gkharish/softdev/DDP/cpp/build && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" /home/gkharish/softdev/DDP/cpp /home/gkharish/softdev/DDP/cpp /home/gkharish/softdev/DDP/cpp/build /home/gkharish/softdev/DDP/cpp/build /home/gkharish/softdev/DDP/cpp/build/CMakeFiles/mainMPC.dir/DependInfo.cmake --color=$(COLOR)
.PHONY : CMakeFiles/mainMPC.dir/depend

