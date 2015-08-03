all:kupyna-reference
kupyna-reference: kupyna.c kupyna.h main.c makefile tables.c tables.h
	gcc kupyna.c main.c tables.c -o kupyna-reference

run:kupyna-reference
	./kupyna-reference
