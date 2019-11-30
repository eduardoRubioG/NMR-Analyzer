make: 
	c++ main.cpp -lm -lGSL -lBLAS

clean: 
	rm a.out *.dat

run: 
	./a.out
