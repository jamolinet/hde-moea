all: hde clean

main.o: main.c global.h rand.h
	gcc -O3 -c  main.c

hungarian.o: hungarian.c global.h rand.h
	gcc -O3 -c  hungarian.c

weights.o: weights.c global.h rand.h
	gcc -O3 -c  weights.c

rand.o: rand.c rand.h
	gcc -O3 -c  rand.c

functions.o: functions.cpp global.h rand.h wfg/ExampleProblems.h wfg/TransFunctions.h
	g++ -O3 -c  functions.cpp

ExampleProblems.o: wfg/ExampleProblems.cpp wfg/ExampleProblems.h wfg/ExampleTransitions.h wfg/ExampleShapes.h
	g++ -O3 -c  wfg/ExampleProblems.cpp

ExampleShapes.o: wfg/ExampleShapes.cpp wfg/ExampleShapes.h wfg/ShapeFunctions.h wfg/FrameworkFunctions.h wfg/Misc.h
	g++ -O3 -c  wfg/ExampleShapes.cpp

ExampleTransitions.o: wfg/ExampleTransitions.cpp wfg/ExampleTransitions.h wfg/Misc.h wfg/TransFunctions.h
	g++ -O3 -c  wfg/ExampleTransitions.cpp

FrameworkFunctions.o: wfg/FrameworkFunctions.cpp wfg/FrameworkFunctions.h wfg/Misc.h
	g++ -O3 -c  wfg/FrameworkFunctions.cpp

Misc.o: wfg/Misc.cpp wfg/Misc.h
	g++ -O3 -c  wfg/Misc.cpp

ShapeFunctions.o: wfg/ShapeFunctions.cpp wfg/ShapeFunctions.h wfg/Misc.h
	g++ -O3 -c  wfg/ShapeFunctions.cpp
 
TransFunctions.o: wfg/TransFunctions.cpp wfg/TransFunctions.h wfg/Misc.h
	g++ -O3 -c  wfg/TransFunctions.cpp

hde: main.o hungarian.o weights.o rand.o functions.o ExampleProblems.o ExampleShapes.o ExampleTransitions.o FrameworkFunctions.o Misc.o ShapeFunctions.o TransFunctions.o
	g++ -O3 main.o hungarian.o weights.o rand.o functions.o ExampleProblems.o ExampleShapes.o ExampleTransitions.o FrameworkFunctions.o Misc.o ShapeFunctions.o TransFunctions.o -lm -o hde

clean:
	rm -rf *.o

 
