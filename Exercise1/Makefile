CC=g++

INCDIR=include
OBJDIR=src

CFLAGS = -c -I$(INCDIR)
OBJECTS:= main.o $(patsubst %.cpp,%.o,$(wildcard src/*.cpp))

output.x: $(OBJECTS)
	$(CC) $(OBJECTS) -o $@
main.o: main.cpp
	$(CC) $(CFLAGS) main.cpp
$(OBJDIR)/%.o: $(OBJDIR)/%.cpp
	$(CC) $(CFLAGS) $^ -o $@
clean:
	rm -rf *.o *.x
	rm $(OBJDIR)/*.o
