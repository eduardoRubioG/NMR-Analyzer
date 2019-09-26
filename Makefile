COMP = g++ 
LDLIBS = -lm 
OBJS = hello.o
HEADERS = 

main : main.o $(OBS)
	$(COMP) main.o -o main $(OBJS) $(LDLIBS)

main.o : main.cpp $(HEADERS)
	$(COMP) $(CCFLAG) -c main.cpp

hello.o : hello.cpp $(HEADERS)
	$(COMP) -c hello.cpp

clean: 
	rm *.o
