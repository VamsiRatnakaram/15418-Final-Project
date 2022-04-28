APP_NAME=seq_aStar
OBJS += seq_aStar.o

APP_NAME_MP=openMP_aStar
OBJS_MP += openMP_aStar.o

CXX = g++ -m64 -std=c++11
CXXFLAGS = -I../common -O3 -Wall -Wextra

CXX_MP = g++ -m64 -std=c++11
CXXFLAGS_MP = -I. -O3 -Wall -fopenmp -Wno-unknown-pragmas

default: $(APP_NAME)

$(APP_NAME): $(OBJS)
	$(CXX) $(CXXFLAGS) -o $@ $(OBJS)

openMP:  $(APP_NAME_MP)

$(APP_NAME_MP): $(OBJS_MP)
	$(CXX_MP) $(CXXFLAGS_MP) -o $@ $(OBJS_MP)


seq_aStar.o: seq_aStar.cpp
	$(CXX) $< $(CXXFLAGS) -c -o $@

openMP_aStar.o: openMP_aStar.cpp
	$(CXX_MP) $< $(CXXFLAGS_MP) -c -o $@

clean:
	/bin/rm -rf *~ *.o $(APP_NAME) *.class
	/bin/rm -rf *~ *.o $(APP_NAME_MP) *.class
