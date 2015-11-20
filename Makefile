
all: openmp rtm

rtm : rtm.c
	icc -O2 -o rtm rtm.c

openmp: rtm_openmp.cpp
	icpc -O2 -o rtm_openmp -openmp rtm_openmp.cpp

clean:
	rm -f rtm.o rtm rtm_openmp forward_wavefield.bin Wavelet.bin

