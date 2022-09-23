CC=g++
LDIR=./lib
CDIR=./src
ODIR=./src/obj

DEPS_=vector.h matrix.h qmclib.h
DEPS = $(patsubst %,$(LDIR)/%,$(DEPS_))

OBJ_=main.o vector.o matrix.o qmclib.o
OBJ = $(patsubst %,$(ODIR)/%,$(OBJ_))

FLAGS=-I$(LDIR) -Wall


$(ODIR)/%.o: $(CDIR)/%.cpp $(DEPS)
	$(CC) $(FLAGS) -c -o $@ $<

qdotvmc: $(OBJ)
	$(CC) $(FLAGS) -o $@ $^

.PHONY: clean

clean:
	rm -f $(ODIR)/*.o
