EXEC   := PDFS
OBJDIR := obj
BINDIR := bin
SRCDIR := src

UOBJ   := $(patsubst $(UTLDIR)/%.c,$(OBJDIR)/%.o,$(wildcard $(UTLDIR)/*.c))
SOBJ   := $(patsubst $(SRCDIR)/%.c,$(OBJDIR)/%.o,$(wildcard $(SRCDIR)/*.c))
COBJ   := $(patsubst $(SRCDIR)/%.cu,$(OBJDIR)/%.o,$(wildcard $(SRCDIR)/*.cu)) # cuda files

TARGET := $(BINDIR)/$(EXEC)

# Compiler options
CC = gcc -std=c99
CCFLAGS =  -Wall -g -lm

NVCC = /usr/local/cuda65/toolkit/6.5.14/bin/nvcc
NVFLAGS = -arch=sm_20

all: makedirectories $(UOBJ) $(SOBJ) $(COBJ)
	$(NVCC) $(INFLAGS) $(NVFLAGS) -o $(TARGET) $(COBJ) $(UOBJ) $(SOBJ)

$(OBJDIR)/%.o: $(UTLDIR)/%.c
	$(CC) $(CCFLAGS) -o $@ -c $<

$(OBJDIR)/%.o: $(SRCDIR)/%.c
	$(CC) $(CCFLAGS) -o $@ -c $<

$(OBJDIR)/%.o: $(SRCDIR)/%.cu
	$(NVCC) $(NVFLAGS) -o $@ -c $< 
	
makedirectories:
	mkdir -p $(OBJDIR) $(BINDIR)

clean:
	rm -f $(OBJDIR)/*.o $(TARGET)

clobber:
	rm -rf $(OBJDIR) $(BINDIR) *.vti
