COMP = g++ 
LDLIBS = -lm 
OBJS =
HEADERS = 

main : main.o $(OBS)
	$(COMP) main.o -o main $(OBJS) $(LDLIBS)

main.o : main.cpp $(HEADERS)
	$(COMP) $(CCFLAG) -c main.cpp

clean: 
	rm *.o
