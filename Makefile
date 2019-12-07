make: 
	g++ main.cpp -gdwarf-3 -lm -lgsl -lblas

clean: 
	rm a.out *.dat

run: 
	./a.out
