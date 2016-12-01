#
# make [log=NO|ERROR|INFO|DEBUG]
#
EXEC := PDFS

# Directories
SRC := src
OBJ := obj
BIN := bin

# Modify the location of these directories
MPI := $(OPENMPIDIR)
HDF := $(HDF5DIR)
SZP := $(SZIPDIR)
NCF := $(NETCDFDIR)

OBJS := $(patsubst $(SRC)/%.c,$(OBJ)/%.o,$(wildcard $(SRC)/*.c))

# Default log level
LOG := INFO_LEVEL

# Check args
ifeq ($(log),NO)
	LOG := NO_LOG
endif
ifeq ($(log),ERROR)
	LOG := ERROR_LEVEL
endif
ifeq ($(log),DEBUG)
	LOG := DEBUG_LEVEL
endif

# Compiler options
CC := pgcc -fast -Minfo=accel -ta=tesla:managed
CCFLAGS := -g -lm 
LDFLAGS := -L$(MPI)/lib -lmpi -L$(HDF)/lib -lhdf5 -L$(NCF)/lib -lnetcdf -L$(SZP)/lib -lsz
ALFLAGS := $(CCFLAGS) $(LDFLAGS)

TARGET  := $(BIN)/$(EXEC)

default: mkdirs $(OBJS)
	$(CC) $(ALFLAGS) -DLOG_LEVEL=$(LOG) -o $(TARGET) $(OBJS)

$(OBJ)/%.o: $(SRC)/%.c
	$(CC) $(ALFLAGS) -DLOG_LEVEL=$(LOG) -o $@ -c $<

mkdirs:
	mkdir -p $(OBJ) $(BIN)

clean:
	rm -f $(TARGET)

clobber:
	rm -rf $(BIN) $(OBJ) pdfs_*.log
