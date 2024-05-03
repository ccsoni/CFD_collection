CC=gcc
CFLAGS=-std=c99 -ffast-math -funroll-loops -O3
#CFLAGS+=-D__FVS_MUSCL__
#CFLAGS+=-D__FVS_MUSCL2__
#CFLAGS+=-D__FVS_MP5__
#CFLAGS+=-D__FVS_MP5_2__
#CFLAGS+=-D__FVS_MP5SL__
#CFLAGS+=-D__AUSMP__
#CFLAGS+=-D__AUSMP_MUSCL__
#CFLAGS+=-D__AUSMP_MP5__
#CFLAGS+=-D__HLL__
#CFLAGS+=-D__HLLC__
#CFLAGS+=-D__HLLC_MUSCL__
CFLAGS+=-D__HLLC_MP5__
#CFLAGS+=-D__GODUNOV__
CFLAGS+=-D__TVD_RK3__
INIT=-D__SOD_INIT__
#INIT=-D__LAX_INIT__
#INIT=-D__SHU_INIT__
#INIT=-D__LE_BLANC_INIT__
NMESH=-D__NMESH__=128
CFLAGS += $(FLUX) $(INIT) $(NMESH)
CFLAGS += -D__ENTROPY_INTEGRATION__
CFLAGS += -D__POSITIVITY_LIMITER__

FLUID_1D_OBJ = fluid.o main.o muscl.o mp5.o semi_lagrange.o ausm.o
FLUID_1D_DEP := $(FLUID_1D_OBJ)

all: fluid_1d

fluid_1d: $(FLUID_1D_DEP)
	$(CC) $(CFLAGS) -o $@ $(FLUID_1D_OBJ) -lm

clean:
	\rm -rf fluid_1d
	\rm -rf *.o

distclean: clean
	\rm -rf *~
