EXEC := PDFS

# Directories
SRC := src
OBJ := obj
BIN := bin

OBJS := $(patsubst $(SRC)/%.c,$(OBJ)/%.o,$(wildcard $(SRC)/*.c))

# Compiler options
CC := pgcc -fast -Minfo=accel -ta=tesla:managed
CCFLAGS := -g -lm
ICFLAGS := -L/cm/shared/apps/openmpi/openacc-15.7/64/slurm/2.0.1/lib -lmpi -L$(NETCDFDIR)/lib -lnetcdf -I/cm/shared/apps/hdf5_18-openmpi/1.8.17/include -L/cm/shared/apps/hdf5_18-openmpi/1.8.17/lib -lhdf5 -I/cm/shared/apps/szip/1.0/include -L/cm/shared/apps/szip/1.0/lib -lsz
ALFLAGS := $(CCFLAGS) $(ICFLAGS)

TARGET  := $(BIN)/$(EXEC)

default: mkdirs $(OBJS)
	$(CC) $(ALFLAGS) -o $(TARGET) $(OBJS)

$(OBJ)/%.o: $(SRC)/%.c
	$(CC) $(ALFLAGS) -o $@ -c $<

mkdirs:
	mkdir -p $(OBJ) $(BIN)

clean:
	rm -f $(TARGET)

clobber:
	rm -rf $(BIN) $(OBJ) GIN3D-PDFP_?.log

cmisc: clobber
	rm -rf dump* bol_* 
