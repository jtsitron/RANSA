# Makefile for RANSA

CFLAGS	= -g
LIBS	= /usr/lib64/libm.so /usr/lib64/libgsl.so /usr/lib64/libgslcblas.so
OBJS	= options.o input.o mobject.o nested.o
CC      = g++

########################################################

# executables
sample: sample.c $(OBJS)
	$(CC) $(CFLAGS) -o $@ $^ $(LIBS)
	
# object files
$(OBJS): %.o: %.c
	$(CC) $(CFLAGS) -c -o $@ $< 

.PHONY: clean
clean:
	-rm *.o 
	-rm sample






