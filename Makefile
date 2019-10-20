# makefile for motif-search 

CC = g++
FLG = -g -Wall -Wextra -Werror -std=c++17  
OBJ = main.o 
PRG = motifsearch 
LIB = -lstdc++fs 

all: ${PRG}

opt: FLG = -O3 -std=c++17  
opt: all

$(PRG): $(OBJ)
	$(CC) $(FLG) $(OBJ) $(LIB) -o $@

main.o: 	
	$(CC) $(FLG) -c main.cpp $(LIB)

obj:
	rm $(OBJ)

clean:
	rm $(OBJ) $(PRG)

#.PHONY: all obj
