MolDyn_NVE.exe: MolDyn_NVE.cpp ../../ParallelRandomNumberGenerator/random.o
	g++ -Wall -o $@ $^
%.o: %.cpp %.h
	g++ -Wall -c -o $@ $<
clean:
	rm *.o
	rm main.exe
	rm seed.out
