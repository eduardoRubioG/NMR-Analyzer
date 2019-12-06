make: 
	c++ main.cpp -lm -lgsl -lblas

clean: 
	rm a.out *.dat

run: 
	./a.out
