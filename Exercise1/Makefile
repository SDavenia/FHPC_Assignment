RN=srun
CC=mpicc -fopenmp 

INCDIR=src
OBJDIR=src

CFLAGS = -c -I$(INCDIR)
OBJECTS:= main_parallel.o $(patsubst %.c,%.o,$(wildcard src/*.c))

main_parallel.exe: $(OBJECTS)
	$(RN) $(CC) $(OBJECTS) -o $@
main_parallel.o: main_parallel.c
	$(RN) $(CC) $(CFLAGS) main_parallel.c
$(OBJDIR)/%.o: $(OBJDIR)/%.c
	$(RN) $(CC) $(CFLAGS) $^ -o $@
clean:
	rm -rf *.o *.exe
	rm $(OBJDIR)/*.o
clean_o:
	rm $(OBJDIR)/*.o
	rm -rf *.o