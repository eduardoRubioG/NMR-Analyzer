make: 
	c++ main.cpp -lm

clean: 
	rm a.out analysis.txt filter.txt roots.dat spline.dat

run: 
	./a.out dat1.txt
