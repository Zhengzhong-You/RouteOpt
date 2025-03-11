###############################################################################
# Matthew Galati - cvrpsep.a - 10/24/03
###############################################################################

# Makefile for Lysgaards CVRPSEP

###############################################################################
#-----  UNAME           define the shell name
#-----  OS              define the OS name (short form)
#-----  OPTFLAG         define the optimization/debug level
#-----	COM_TARGETDIR	define the directory to compile the objects into
#-----	COM_DEPDIR	define the directory for the dependency files
#-----	COM_LIBSRC	define the source files to be compiled
#-----	COM_LIBOBJ	define the object files to be put in the library
#-----	COM_LIBDEP	define the dependency files for these object files
#-----  TARGET_LIB      define the name of the library
#-----  CXXFLAGS	define the compiler flags (not including include directories)
#-----  DEPFLAGS	define the dependency flags (mainly for include directories)
#-----  CXX             define the compile (dependent on OS)

UNAME = $(shell uname)
OS = $(shell uname -s)
OPTFLAG = -O3
#OPTFLAG = -g
COM_TARGETDIR = obj
COM_DEPDIR = dep
COM_LIBSRC = $(filter-out %unitTest.cpp, $(shell /bin/ls *.cpp))
COM_LIBOBJ = $(addprefix $(COM_TARGETDIR)/, $(COM_LIBSRC:.cpp=.o))
COM_LIBDEP = $(addprefix $(COM_DEPDIR)/, $(COM_LIBSRC:.cpp=.d))
TARGET_LIB = libcvrpsep
CXXFLAGS = $(OPTFLAG)
DEPFLAGS += -I. -D$(OS)

#LinuxCXX=g++
#CXX=$($(OS)CXX)
CXX=g++

###############################################################################
# Create the dependency information
$(COM_TARGETDIR)/%.o : %.cpp ${COM_DEPDIR}/%.d Makefile
	@echo ""
	@if test ! -e ${COM_DEPDIR}/$*.d ; then \
	    echo ; \
	    echo "   ${COM_DEPDIR}/$*.d is missing."; \
	    echo "   Probably a header file was not found when make examined";\
	    echo "   $*.cpp in an attempt to create that dependency."; \
	    echo ; \
	    exit 1; \
	fi
	@mkdir -p $(COM_TARGETDIR)
	$(CXX) $(DEPFLAGS) $(CXXFLAGS) -c $< -o $@

${COM_DEPDIR}/%.d : %.cpp
	@echo Creating dependency $*.d
	@mkdir -p ${COM_DEPDIR}
	@rm -f $*.d $*.dd
	g++ -MM $(DEPFLAGS) $< > $*.dd
	@sed -e "s|$*.o|$(COM_DEPDIR)/$*.d $(COM_TARGETDIR)/$*.o|g" $*.dd > $@
	@rm -f $*.dd

###############################################################################
# Create the targets

.PHONY: default all ${TARGET_LIB} clean

default: ${TARGET_LIB}

all : ${TARGET_LIB}

doc :
	${MAKE} -f Makefile.doc doc

clean :
	@rm -rf $(COM_DEPDIR)
	@rm -rf $(COM_TARGETDIR)
	@rm -rf $(UNAME)-*


${TARGET_LIB} : $(COM_TARGETDIR)/${TARGET_LIB}.a

$(COM_TARGETDIR)/${TARGET_LIB}.a: $(COM_LIBOBJ) Makefile
	@echo ""
	@echo Creating library $(notdir $@)
	@mkdir -p $(COM_TARGETDIR)
	@rm -f $@
	$(AR) -q $@ $(COM_LIBOBJ)

-include $(COM_LIBDEP) $(TESTDEP)

