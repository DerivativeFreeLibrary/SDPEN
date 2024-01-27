# makefile for compiling seqpen source code
CC=gcc
LD=gcc
RM=rm -f
CFLAGS=-Wall -ansi -pedantic -g
# -g for debugging, -O for optimization

#*******************************************************
# Seqpen code

EXEC = \
       seqpen

OBJS = \
       seqpenmain.o \
       seqpen.o \
       problem.o


SRCS = \
       seqpenmain.c \
       seqpen.c \
       problem.c
#*******************************************************
# Seqpen header

HDRS = \
           seqpen.h \
           basic.h
#*******************************************************

all:	seqpen

#*******************************************************

seqpen:	$(OBJS)
	$(LD) $(LDFLAGS) $(OBJS) -o seqpen -lm

#*******************************************************

seqpenmain.o:	seqpenmain.c $(HDRS)
	$(CC) $(CFLAGS) -c seqpenmain.c

seqpen.o:	seqpen.c $(HDRS)
	$(CC) $(CFLAGS) -c seqpen.c

problem.o:	problem.c $(HDRS) 
	$(CC) $(CFLAGS) -c problem.c

#*******************************************************

