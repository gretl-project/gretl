CXX ?= g++
CFLAGS = -Wall -Wconversion -g -O3 -fPIC -fopenmp

svm.o: svm.cpp svm.h
	$(CXX) $(CFLAGS) -c svm.cpp

clean:
	rm -f svm.o
