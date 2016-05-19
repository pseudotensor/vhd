#
# Set USEICC=1 and USEGCC=0 to use Intel's icc compiler
# Set USEGCC=1 and USEICC=0 to use gcc compiler
USEICC=1
USEGCC=0


#General env. vars
RM = rm -f
SUFF = c
#
# Define a cpp source file directory, a fortran file subdirectory,
# an object file subdirectory and an executable file subdirectory
#
SRCD = .
OBJD = .
OBJD2 = ./Obj
BIND = ./bin

CMD = $(BIND)/twod

# Define preprocessor and compile flags, and linked libraries
#

ifeq ($(USEICC),1)
CFLAGS = -O3 -tpp7 -axiMKW -ipo -unroll -w1 -wd=175
LDFLAGS = -lm
MYCC = icc
endif

ifeq ($(USEGCC),1)
CFLAGS = -Wall -mcpu=pentiumpro -march=pentiumpro -O4 -ffast-math -finline-functions 
LDFLAGS = -lm
MYCC=gcc  
endif

CC = $(MYCC) $(CFLAGS)


#
# Define source files
#
SRCS = \
$(SRCD)/analsol.$(SUFF) \
$(SRCD)/bound.$(SUFF) \
$(SRCD)/diag.$(SUFF) \
$(SRCD)/init.$(SUFF) \
$(SRCD)/main.$(SUFF) \
$(SRCD)/numerics.$(SUFF) \
$(SRCD)/ranc.$(SUFF) \
$(SRCD)/step.$(SUFF) \
$(SRCD)/sweep.$(SUFF) \
$(SRCD)/timestep.$(SUFF) \
$(SRCD)/utilfun.$(SUFF)

#
# Define object files
#                                         
OBJS = \
$(OBJD)/analsol.o \
$(OBJD)/bound.o \
$(OBJD)/diag.o \
$(OBJD)/init.o \
$(OBJD)/main.o \
$(OBJD)/numerics.o \
$(OBJD)/ranc.o \
$(OBJD)/step.o \
$(OBJD)/sweep.o \
$(OBJD)/timestep.o \
$(OBJD)/utilfun.o




#
all: 	$(OBJD) $(BIND) $(CMD)
#
$(OBJD):
	mkdir $(OBJD)

$(BIND):
	mkdir $(BIND)

$(CMD): $(OBJS) $(SRCS) makefile
	$(CC) $(CFLAGS) -o $(CMD) $(OBJS) $(LDFLAGS)

# dependencies
$(OBJS)	: $(SRCD)/defs.h $(SRCD)/global.h $(SRCD)/makefile

$(OBJD)/timestep.o : $(SRCD)/timestep1.h $(SRCD)/timestep2.h


cleandat:
	$(RM) $(BIND)/*.dat $(BIND)/*.dat.ras
 
clean:
	$(RM) *.o $(CMD)

