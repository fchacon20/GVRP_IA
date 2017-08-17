OBJS = main.o
CC = g++
CFLAGS = -Wall -c
LFLAGS = -Wall -std=c++11

SA : $(OBJS)
	$(CC) $(LFLAGS) $(OBJS) -o SA

main.o : city.h haversine.h solutionUtilities.h main.cpp
	$(CC) $(CFLAGS) -std=c++11 main.cpp

haversine.o : haversine.h haversine.cpp
	$(CC) $(CFLAGS) haversine.cpp

solutionUtilities.o : city.h solutionUtilities.h solutionUtilities.cpp
	$(CC) $(CFLAGS) solutionUtilities.cpp

clean:
	\rm main.o SA
