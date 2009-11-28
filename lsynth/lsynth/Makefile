CC=gcc

CFLAGS=-ggdb

AR = ar
RANLIB = ranlib

SRCS = lsynthcp.c band.c hose.c curve.c mathlib.c
OBJS = $(SRCS:.c=.o)

all	: lsynthcp

lsynthcp:   $(OBJS) $(LIBS)
	$(CC) $(CFLAGS) $(OBJS) -lm -o lsynthcp

# The following nonsense is for adding an icon to the windows executable.
lsynthico.o: lsynth.rc
	windres -i lsynth.rc -o lsynthico.o

lsynth4win:   $(OBJS) lsynthico.o $(LIBS)
	$(CC) $(CFLAGS) $(OBJS) lsynthico.o -lm -o lsynthcp
