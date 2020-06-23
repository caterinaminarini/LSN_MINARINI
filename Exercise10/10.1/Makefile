main.exe: main.cpp ../../ParallelRandomNumberGenerator/random.o classi.o
	g++ -Wall -o $@ $^
%.o: %.cpp %.h
	g++ -Wall -c -o $@ $<
clean:
	rm *.o
	rm main.exe
	rm seed.out
