OBJS = common.o dio.o matlabfunctions.o fft.o old_matlabfunctions.o  platinum.o star.o synthesis.o macresmain.o old_fft.o utils.o wavread.o
CXX = g++
CFLAGS = -g -O -msse3
INCLUDES = -I.
LIBS = 

macres : $(OBJS)  mmsystem.h windows.h wavread.h world.h
	$(CXX) -o $@ $(OBJS) $(LIBS) $(CFLAGS)

.cpp.o:
	$(CXX) $(CFLAGS) $(INCLUDES) -c $<

clean:
	rm macres $(OBJS)
