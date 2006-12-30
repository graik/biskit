# own makefile since delivered makefile does not work

PYTHON = /usr
VERSION = 2.4
CC = gcc
CFLAGS = -Wall -w1 -O3 -tpp7 -vec_report0 
##CFLAGS = -Wall -O3 

## SuSE puts libgpvm3.a into a folder that is not searched for libraries.. added to -L
LFLAGS = -L. -L$(PVM_ROOT)/lib/$(PVM_ARCH) \
	-L$(PVM_ROOT)/pvmgs/$(PVM_ARCH) \
	-lm -lpvm3 -lgpvm3

IFLAGS = -I$(PYTHON)/local/include/python$(VERSION) \
	 -I./include -I$(PVM_ROOT)/include 
	

all: 		modules install

modules:	pypvm_coremodule.so

clean: 
	rm -R -f *.o pypvm_coremodule.so

pypvm_coremodule.o:	pypvm_coremodule.c
			$(CC) -c pypvm_coremodule.c \
			$(CFLAGS) $(IFLAGS)

pypvm_coremodule.so:	pypvm_coremodule.o
			$(CC) pypvm_coremodule.o \
			-o pypvm_coremodule.so \
			-shared -$(LFLAGS)

install:
##	cp pypvm_coremodule.so ../../py/modules
