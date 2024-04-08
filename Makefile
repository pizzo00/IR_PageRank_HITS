OBJS	= main.o Sorter.o
SOURCE	= main.cpp Sorter.cpp
HEADER	= Sorter.h
OUT	= IR_Project
FLAGS	= -g -c -Wall -O3 -std=c++20
LFLAGS	=
CC	= g++

all:	IR_Project

IR_Project: $(OBJS)
	$(CC) -g $(OBJS) -o $(OUT) $(LFLAGS)

main.o: main.cpp
	$(CC) $(FLAGS) main.cpp

Sorter.o: Sorter.cpp
	$(CC) $(FLAGS) Sorter.cpp

clean:
	rm -f $(OBJS) $(OUT)
