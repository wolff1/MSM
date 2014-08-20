# Use icc for Intel Composer XE 2013
#CC=icc
# Compiler flags
#CFLAGS=-Iinclude -c -Wall

#all: msm

#msm: main.o gamma.o output.o phi.o stdafx.o utility.o
#	$(CC) main.o gamma.o output.o phi.o stdafx.o utility.o -o msm

#main.o: main.cpp
#	$(CC) $(CFLAGS) main.cpp

#gamma.o: gamma.cpp
#    $(CC) $(CFLAGS) gamma.cpp

#output.o: output.cpp
#    $(CC) $(CFLAGS) output.cpp

#phi.o: phi.cpp
#    $(CC) $(CFLAGS) phi.cpp

#stdafx.o: stdafx.cpp
#    $(CC) $(CFLAGS) stdafx.cpp

#utility.o: utility.cpp
#    $(CC) $(CFLAGS) utility.cpp

all: src/main.cpp src/gamma.cpp src/output.cpp src/phi.cpp src/stdafx.cpp src/utility.cpp src/operator.cpp src/phiC1.cpp
	icc -Iinclude -mkl src/main.cpp src/gamma.cpp src/output.cpp src/phi.cpp src/stdafx.cpp src/utility.cpp src/operator.cpp src/phiC1.cpp
	mv a.out bin/msm
clean:
	rm -rf *o bin/msm
