
CXX=g++ -O3 -Wall -std=c++11 -pthread -march=native -lntl -lgmp -lm
ITPP=`itpp-config --cflags` `itpp-config --libs`
#ITPP=`itpp-config --cflags` `itpp-config`
#make func.o
#make main.o
bike:

	$(CXX) $(ITPP) -o bike.out main.o func.o

clean:
	rm *.o
	rm *.out