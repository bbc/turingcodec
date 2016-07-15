# List to '.PHONY' all fake targets, those that are neither files nor folders.
# "all" and "clean" are good candidates.
.PHONY: all, clean, git-describe.h

# Define the final program name
PROGNAME := pluto

# Pre-processor flags to be used for includes (-I) and defines (-D) 
CPPFLAGS := -I.. 

# common options 
OPT=-O3 -flto -march=corei7 

# CXXFLAGS is used for C++ compilation options.
CXXFLAGS += -std=c++11 -w $(OPT)

# LDFLAGS is used for linker (-g enables debug symbols)
LDFLAGS  += $(OPT) -fwhole-program

# Which Boost modules to use
BOOST_MODULES = 	\
  chrono \
  filesystem \
  program_options \
  system \
  timer \

# Boost libraries' type (a suffix)
BOOST_MODULES_TYPE := .a

# Define library names with their type
BOOST_MODULES_LIBS := $(addsuffix $(BOOST_MODULES_TYPE),$(BOOST_MODULES))

# Define the linker argument to use the Boost libraries.
BOOST_LDFLAGS := $(addprefix -l:libboost_,$(BOOST_MODULES_LIBS))

# Feed compiler/linker flags with Boost's
CPPFLAGS += $(BOOST_CPPFLAGS)
LDFLAGS += $(BOOST_LDFLAGS)

# List the project' sources to compile or let the Makefile recognize
# them for you using 'wildcard' function.
#
SOURCES_CPP = $(wildcard *.cpp) 
SOURCES_C = $(wildcard *.c) 

# List the project' headers or let the Makefile recognize
# them for you using 'wildcard' function.
#
HEADERS = $(wildcard %.h) $(wildcard %.hpp)

# Construct the list of object files based on source files using
# simple extension substitution.
OBJECTS = $(SOURCES_CPP:%.cpp=%.o) $(SOURCES_C:%.c=%.o)
 
# Now declare the dependencies rules and targets
#
# Starting with 'all' make it  becomes the default target when none 
# is specified on 'make' command line.
all : $(PROGNAME)

# Declare that the final program depends on all objects
$(PROGNAME) : $(OBJECTS)
	g++ -o $@ -L/usr/local/lib -I../assembly $(OBJECTS) $(LDFLAGS) -lpthread ../hevcasm/src/lib/libhevcasm.a

# Don't use implicit rules
%.o: %.cpp $(HEADERS)
	g++ $(CXXFLAGS) $(CPPFLAGS) -c $(OUTPUT_OPTION) $<

%.o: %.c $(HEADERS)
	g++ $(CFLAGS) $(CPPFLAGS) -c $(OUTPUT_OPTION) $<

main.o : git-describe.h

git-describe.h :
	echo "#define GIT_DESCRIBE \\" > $@
	git describe >> $@

# Simple clean-up target
# notes:
# - the '@' before 'echo' informs make to hide command invocation.
# - the '-' before 'rm' command to informs make to ignore errors.
clean :
	@echo "Clean."
	-rm -f git-describe.h *.o $(PROGNAME)

