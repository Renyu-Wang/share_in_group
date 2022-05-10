INC_DIR=weilei_lib

CXX=g++ -O3 -Wall -std=c++11 -pthread -march=native -lntl -lgmp -lm
ITPP=`itpp-config --cflags` `itpp-config --libs`
object_files=$(INC_DIR)/mm_read.o $(INC_DIR)/mmio.o $(INC_DIR)/mm_write.o $(INC_DIR)/lib.o $(INC_DIR)/dist.o $(INC_DIR)/product_lib.o $(INC_DIR)/bp.o 
header_files=$(INC_DIR)/mm_read.h $(INC_DIR)/mmio.h $(INC_DIR)/mm_write.h $(INC_DIR)/lib.h $(INC_DIR)/dist.h $(INC_DIR)/product_lib.h $(INC_DIR)/bp.h

#make object file for target file
%.o:%.cpp $(header_files)
	$(CXX) $(START) $(END) -c $<
#compile object files for lib files
lib:
	cd weilei_lib && make all
#now link object files
%.out:%.o $(object_files)
	$(CXX) $(ITPP) -o $@ $< $(object_files)

# cmd=make lib && make $@.o && make $@.out

bike:
	make lib
	make func.o
	make main2.o

	$(CXX) $(ITPP) -o bike2.out main2.o func.o $(object_files)

bicycle:func.cpp main.cpp
	#g++ -std=c++11 `itpp-config --cflags` -pthread -march=native -o bicycle.out func.cpp main2.cpp  `itpp-config --libs` -lntl -lgmp -lm 

clean:
	rm *.o
	rm *.out