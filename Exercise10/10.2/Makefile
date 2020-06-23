CC = mpicxx
CFLAGS = -Wall -O3 --std=c++11

main.exe : main.o random.o classi.o
	$(CC) random.o main.o classi.o -o main.exe
main.o : main.cpp main.h
	$(CC) -c main.cpp -o main.o $(CFLAGS)
random.o : random.cpp random.h
	$(CC) -c random.cpp -o random.o $(CFLAGS)
classi.o : classi.cpp classi.h
	$(CC) -c classi.cpp -o classi.o $(CFLAGS)
clean :
	rm *.o main.exe
	rm -rf output_*
	rm -rf seed.out
