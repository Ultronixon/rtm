CC = gcc
CFLAGS = -c -Wall -O2 -DNO_BLAS
LDFLAGS = librsf.a -lm
SOURCES =  Mrtm.c step.c wavefun.c
OBJECTS = $(SOURCES:.c=.o)
EXECUTABLE = rtm

all : $(SOURCES) $(EXECUTABLE)

$(EXECUTABLE) : $(OBJECTS)
	$(CC) -o $@ $(OBJECTS) $(LDFLAGS)

.c.o:
	$(CC) $(CFLAGS) $< -o $@

clean:
	rm -f $(EXECUTABLE) $(OBJECTS)
