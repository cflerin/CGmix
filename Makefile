CPP = clang++ # g++
EXECUTABLE = CGmix

CFLAGS = -Wall -O3 -m64
CPPFLAGS = -O3 -D_FILE_OFFSET_BITS=64 -std=c++11
LIB = -lz 

BIN = bin
SOURCE = src

OBJS = src/CGmix.o src/readFiles.o src/hmm.o src/parameters.o

$(BIN)/CGmix: $(OBJS)
	$(CPP) $(CPPFLAGS) $(OBJS) -o $@ $(LIB)

# pull in dependency info for *existing* .o files
-include $(OBJS:.o=.d)

%.o: %.cpp
	$(CPP) -c $(CPPFLAGS) $*.cpp -o $*.o
	$(CPP) -MM $(CPPFLAGS) $*.cpp > $*.d

# remove compilation products
clean:
	@rm -f $(BIN)/$(EXECUTABLE) $(SOURCE)/*.o $(SOURCE)/*.d

