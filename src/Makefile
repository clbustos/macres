OBJS = dio.o matlabfunctions.o platinum.o star.o synthesis.o macresmain.o fft.o wavread.o
CXX = g++
CFLAGS = -g -O -msse3
INCLUDES = -I.
LIBS = 

macres : $(OBJS) fft.h mmsystem.h windows.h wavread.h world.h
	$(CXX) -o $@ $(OBJS) $(LIBS) $(CFLAGS)

.cpp.o:
	$(CXX) $(CFLAGS) $(INCLUDES) -c $<

clean:
	rm macres $(OBJS)
