APP_NAME=seq_aStar
OBJS += seq_aStar.o

CXX = g++ -m64 -std=c++11
CXXFLAGS = -O3 -Wall -Wextra

APP_NAME_MP1=openMP_aStar
OBJS_MP1 += openMP_aStar.o

APP_NAME_MP2=openMP_aStar_v2
OBJS_MP2 += openMP_aStar_v2.o

CXX_MP = g++ -m64 -std=c++11
CXXFLAGS_MP = -I. -g -O3 -Wall -fopenmp -Wno-unknown-pragmas

CXX_MP2 = g++ -m64 -std=c++17 
CXXFLAGS_MP2 = -I. -g -O3 -Wall -fopenmp -Wno-unknown-pragmas -lpthread

APP_NAME_MPI1=openMPI_aStarv2
OBJS_MPI1 += openMPI_aStarv2.o

APP_NAME_MPI2=openMPI_aStarv3
OBJS_MPI2 += openMPI_aStarv3.o

APP_NAME_MPI3=openMPI_aStarv4
OBJS_MPI3 += openMPI_aStarv4.o

CXX_MPI = mpic++ -std=c++11
CXXFLAGS_MPI = -I. -O3 #-Wall -Wextra

default: $(APP_NAME)

$(APP_NAME): $(OBJS)
	$(CXX) $(CXXFLAGS) -o $@ $(OBJS)

openMP1: $(APP_NAME_MP1)

$(APP_NAME_MP1): $(OBJS_MP1)
	$(CXX_MP) $(CXXFLAGS_MP) -o $@ $(OBJS_MP1)

openMP2: $(APP_NAME_MP2)

$(APP_NAME_MP2): $(OBJS_MP2)
	$(CXX_MP2) $(CXXFLAGS_MP2) -o $@ $(OBJS_MP2)

openMPI1: $(APP_NAME_MPI1)

$(APP_NAME_MPI1): $(OBJS_MPI1)
	$(CXX_MPI) $(CXXFLAGS_MPI) -o $@ $(OBJS_MPI1)

openMPI2: $(APP_NAME_MPI2)

$(APP_NAME_MPI2): $(OBJS_MPI2)
	$(CXX_MPI) $(CXXFLAGS_MPI) -o $@ $(OBJS_MPI2)

openMPI3: $(APP_NAME_MPI3)

$(APP_NAME_MPI3): $(OBJS_MPI3)
	$(CXX_MPI) $(CXXFLAGS_MPI) -o $@ $(OBJS_MPI3)

seq_aStar.o: seq_aStar.cpp
	$(CXX) $< $(CXXFLAGS) -c -o $@

openMP_aStar.o: openMP_aStar.cpp
	$(CXX_MP) $< $(CXXFLAGS_MP) -c -o $@

openMP_aStar_v2.o: openMP_aStar_v2.cpp
	$(CXX_MP2) $< $(CXXFLAGS_MP2) -c -o $@

openMPI_aStarv2.o: openMPI_aStarv2.cpp
	$(CXX_MPI) $< $(CXXFLAGS_MPI) -c -o $@

openMPI_aStarv3.o: openMPI_aStarv3.cpp
	$(CXX_MPI) $< $(CXXFLAGS_MPI) -c -o $@

openMPI_aStarv4.o: openMPI_aStarv4.cpp
	$(CXX_MPI) $< $(CXXFLAGS_MPI) -c -o $@

clean:
	/bin/rm -rf *~ *.o $(APP_NAME) *.class
	/bin/rm -rf *~ *.o $(APP_NAME_MP1) *.class
	/bin/rm -rf *~ *.o $(APP_NAME_MP2) *.class
	/bin/rm -rf *~ *.o $(APP_NAME_MPI1) *.class
	/bin/rm -rf *~ *.o $(APP_NAME_MPI2) *.class