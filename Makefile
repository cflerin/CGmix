CPP = g++
EXECUTABLE = CGmix

CFLAGS = -Wall -g # -O2 -m64 
CPPFLAGS = -O2 -D_FILE_OFFSET_BITS=64 
LIB = -lz 

OBJS = CGmix.o readFiles.o hmm.o

CGmix: $(OBJS)
	$(CPP) $(CPPFLAGS) $(OBJS) -o $@ $(LIB)

# pull in dependency info for *existing* .o files
-include $(OBJS:.o=.d)

%.o: %.cpp
	$(CPP) -c $(CPPFLAGS) $*.cpp -o $*.o
	$(CPP) -MM $(CPPFLAGS) $*.cpp > $*.d

# remove compilation products
clean:
	@rm -f CGmix *.o *.d
