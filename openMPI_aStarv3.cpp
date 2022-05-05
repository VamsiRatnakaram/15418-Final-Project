#include <assert.h>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <libgen.h>
#include <math.h>
#include <queue>
#include <vector>
#include <algorithm>
#include <iostream>
#include <sstream>
#include <linux/limits.h>
#include <string>
#include <unistd.h>
#include <stack>
#include <bits/stdc++.h>
#include <chrono>
#include <time.h>
#include <sys/time.h>
#include <mpi.h>

#define BILLION  1E9;

#define TAG_VERT 7
#define TAG_BOOL 6
#define TAG_PAIR 5
#define TAG_REQ  4
#define TAG_CELL 3
#define TAG_DATA 2
#define TAG_DONE 1

// #include "CycleTimer.h"
using namespace std::chrono;

// A C++ Program to implement A* Search Algorithm
// Note: We adpated an existing sequential implementation from www.geeksforgeeks.org
using namespace std;

// Creating a shortcut for int, int pair type
typedef pair<int, int> Pair;

// Creating a shortcut for pair<int, pair<int, int>> type
typedef pair<double, pair<int, int> > pPair;

// A structure to hold the necessary parameters
struct cell {
	// Row and Column index of its parent
	// Note that 0 <= i <= ROW-1 & 0 <= j <= COL-1
	int parent_i, parent_j;
	// f = g + h
	double f, g, h;
};

struct compare{
    bool operator() (const pPair& p1,const pPair& p2 ){
         return p1.first>p2.first;
    }
};

// A Utility Function to check whether given cell (row, col)
// is a valid cell or not.
bool isValid(int row, int col, int dim_x, int dim_y)
{
	// Returns true if row number and column number
	// is in range
	return (row >= 0) && (row < dim_y) && (col >= 0) && (col < dim_x);
}

// A Utility Function to check whether the given cell is
// blocked or not
bool isUnBlocked(int *map, int row, int col, int dim_y)
{
	// Returns true if the cell is not blocked else false
	return (map[row*dim_y + col] == 0);
}

// A Utility Function to check whether destination cell has
// been reached or not
bool isDestination(int row, int col, Pair dest)
{
	return (row == dest.first && col == dest.second);
}

// A Utility Function to calculate the 'h' heuristics.
double calculateHValue(int row, int col, Pair dest)
{
	// Return Manhatten distance
	return abs(dest.second - col) + abs(dest.first - row);
}

void hashedSendFunction(pPair pairToSend, int dim_x, int nproc, int procID, std::priority_queue<pPair, vector<pPair>, compare> *openList, cell cellToSend, cell *cellDetails) {
	int x = pairToSend.second.first;
	int y = pairToSend.second.second;
	int index = x*dim_x + y;

	int indexToSend = index % nproc;

	MPI_Request request = MPI_REQUEST_NULL;
	if (indexToSend == procID) {
        if (cellDetails[x*dim_x + y].f == INT_MAX || (cellDetails[x*dim_x + y].f > pairToSend.first)) {
            openList->push(pairToSend);
            cellDetails[x*dim_x + y] = cellToSend;
        }
	}
	else {
		MPI_Isend((void*)(&pairToSend), 16, MPI_BYTE, indexToSend, TAG_DATA, MPI_COMM_WORLD, &request);
        MPI_Wait(&request, MPI_STATUS_IGNORE);
		MPI_Isend((void*)(&cellToSend), 32, MPI_BYTE, indexToSend, TAG_CELL, MPI_COMM_WORLD, &request);
        MPI_Wait(&request, MPI_STATUS_IGNORE);
	}
}

// A Function to find the shortest path between
// a given source cell to a destination cell according
// to A* Search Algorithm
double aStarSearch(int *map, Pair src, Pair dest, int dim_x, int dim_y, int procID, int nproc, const char* input_filename)
{
    (void) input_filename;
    // StartTime after intialization
    double startTime = MPI_Wtime();
	// If the source is out of range
	if (isValid(src.first, src.second, dim_x, dim_y) == false) {
		printf("Source is invalid\n");
		return -1;
	}

	// If the destination is out of range
	if (isValid(dest.first, dest.second, dim_x, dim_y) == false) {
		printf("Destination is invalid\n");
		return -1;
	}

	// Either the source or the destination is blocked
	if (isUnBlocked(map, src.first, src.second, dim_y) == false
		|| isUnBlocked(map, dest.first, dest.second, dim_y)
			== false) {
		printf("Source or the destination is blocked\n");
		return -1;
	}

	// If the destination cell is the same as source cell
	if (isDestination(src.first, src.second, dest)
		== true) {
		printf("We are already at the destination\n");
		return -1;
	}

	// Declare a 2D array of structure to hold the details
	// of that cell
	// cell cellDetails[dim_x][dim_y];
    cell *cellDetails = (cell*)calloc(dim_x*dim_y, sizeof(cell));

	int i, j;

	for (i = 0; i < dim_y; i++) {
		for (j = 0; j < dim_x; j++) {
			cellDetails[i*dim_y+j].f = INT_MAX;
			cellDetails[i*dim_y+j].g = INT_MAX;
			cellDetails[i*dim_y+j].h = INT_MAX;
			cellDetails[i*dim_y+j].parent_i = -1;
			cellDetails[i*dim_y+j].parent_j = -1;
		}
	}

	// Initialising the parameters of the starting node
	i = src.first, j = src.second;
	cellDetails[i*dim_y+j].f = 0.0;
	cellDetails[i*dim_y+j].g = 0.0;
	cellDetails[i*dim_y+j].h = 0.0;
	cellDetails[i*dim_y+j].parent_i = i;
	cellDetails[i*dim_y+j].parent_j = j;

	/*
	Create an priority queue having information as-
	<f, <i, j>>
	where f = g + h,
	and i, j are the row and column index of that cell
	Note that 0 <= i <= ROW-1 & 0 <= j <= COL-1
	This open list is implemented as a set of pair of
	pair.*/
    std::priority_queue<pPair, vector<pPair>, compare> openList;

	// Put the starting cell on the open list and set its
	// 'f' as 0
    if (procID == 0) {
        openList.push(make_pair(0.0, make_pair(i, j)));
    }

	// We set this boolean value as false as initially
	// the destination is not reached.
    MPI_Request request = MPI_REQUEST_NULL;
    bool localFoundDest = false;    
    int  localDestCost= INT_MAX;
    pPair localDestVertex;
    int conditionFlag = true;
    bool emptyBool[nproc];
	while (conditionFlag) {
        for (int node = 0; node < nproc; node++)
        {
            pPair tmpPair;
            cell tmpCell;
            int flag_data = 0;
            int flag_done = 0;
            int flag_cell = 0;
            if (node != procID)
            {	
                MPI_Iprobe(node, TAG_DATA, MPI_COMM_WORLD, &flag_data, MPI_STATUS_IGNORE);
                MPI_Iprobe(node, TAG_CELL, MPI_COMM_WORLD, &flag_cell, MPI_STATUS_IGNORE);
                while (flag_data && flag_cell)
                {
                    MPI_Irecv((void*)(&tmpPair), 16, MPI_BYTE, node, TAG_DATA, MPI_COMM_WORLD, &request);
                    MPI_Wait(&request, MPI_STATUS_IGNORE);
                    MPI_Irecv((void*)(&tmpCell), 32, MPI_BYTE, node, TAG_CELL, MPI_COMM_WORLD, &request);
                    MPI_Wait(&request, MPI_STATUS_IGNORE);

                    int x = tmpPair.second.first;
                    int y = tmpPair.second.second;
                    int index = x*dim_x + y;

                    int indexToSend = index % nproc;

                    if (indexToSend == procID) {
                        if (cellDetails[x*dim_x + y].f == INT_MAX || (cellDetails[x*dim_x + y].f > tmpPair.first)) {
                            openList.push(tmpPair);
                            cellDetails[x*dim_y+y] = tmpCell;
                        }
                    }
                    MPI_Iprobe(node, TAG_DATA, MPI_COMM_WORLD, &flag_data, MPI_STATUS_IGNORE);
                    MPI_Iprobe(node, TAG_CELL, MPI_COMM_WORLD, &flag_cell, MPI_STATUS_IGNORE);
                }
                MPI_Iprobe(node, TAG_DONE, MPI_COMM_WORLD, &flag_done, MPI_STATUS_IGNORE);
                pair<bool, int> tmp;
                if (flag_done){
                    MPI_Irecv((void*)(&tmp), 8, MPI_BYTE, node, TAG_DONE, MPI_COMM_WORLD, &request);
                    MPI_Wait(&request, MPI_STATUS_IGNORE);
                    conditionFlag=false;
                    break;
                }
            }
        }
        if(!conditionFlag){
            break;
        }

        if (openList.size() != 0) {
            pPair p = openList.top();
            openList.pop();

            // Add this vertex to the closed list
            i = p.second.first;
            j = p.second.second;

            // To store the 'g', 'h' and 'f' of the 4 successors
            double gNew, hNew, fNew;
            for (int z = 0; z < 4; z++) {
                // Only process this cell if this is a valid one
                int x, y;
                if (z==0) {
                    x = i-1;
                    y = j;
                }
                else if (z==1) {
                    x = i+1;
                    y = j;
                }
                else if (z==2) {
                    x = i;
                    y = j+1;
                }
                else {
                    x = i;
                    y = j-1;
                }
                
                if (isValid(x, y, dim_x, dim_y) == true) {
                    // If the destination cell is the same as the
                    // current successor
                    if (isDestination(x, y, dest) == true) {
                        localFoundDest = true;
                        //update if you have better path
                        int fnew=cellDetails[i*dim_y+j].g+1.0;
                        localDestCost=fnew;
                        localDestVertex=make_pair(fnew, make_pair(i, j));
                        conditionFlag=false;
                        for (int node = 0; node < nproc; node++) {
                                if (node != procID) {
                                    pair<bool, int> tmp = make_pair(localFoundDest, localDestCost);
                                    MPI_Isend((void*)(&tmp), 8, MPI_BYTE, node, TAG_DONE, MPI_COMM_WORLD, &request);
                                    MPI_Wait(&request, MPI_STATUS_IGNORE);
                                }
                        }
                        break;
                    }
                    else if (isUnBlocked(map, x, y, dim_y) == true) {
                        gNew = cellDetails[i*dim_y+j].g + 1.0;
                        hNew = calculateHValue(x, y, dest);
                        fNew = gNew + hNew;

                        pPair pairToSend = make_pair(fNew, make_pair(x, y));
                        cell cellToSend;
                        // Update the details of this cell
                        cellToSend.f = fNew;
                        cellToSend.g = gNew;
                        cellToSend.h = hNew;
                        cellToSend.parent_i = i;
                        cellToSend.parent_j = j;
                        hashedSendFunction(pairToSend, dim_x, nproc, procID, &openList, cellToSend, cellDetails);        
                    }
                }
            }

        }
    }

    MPI_Barrier(MPI_COMM_WORLD);
    // EndTime before I/O
    double endTime = MPI_Wtime();
	double computeTime = endTime - startTime;

    MPI_Status status;

    pPair localDestVertexArray[nproc];
    bool localFoundDestArray[nproc];
    MPI_Gather(&localDestVertex, 16, MPI_BYTE, localDestVertexArray, 16, MPI_BYTE, 0, MPI_COMM_WORLD);
    MPI_Gather(&localFoundDest, 1, MPI_BYTE, localFoundDestArray, 1, MPI_BYTE, 0, MPI_COMM_WORLD);

    pPair globalDestVertex = make_pair(INT_MAX, make_pair(0,0));

    for (int node = 0; node < nproc; node++) {
        if(localFoundDestArray[node]) {
            if (localDestVertexArray[node].first < globalDestVertex.first || globalDestVertex.first == INT_MAX) {
                globalDestVertex.first = localDestVertexArray[node].first;
                globalDestVertex.second.first = localDestVertexArray[node].second.first;
                globalDestVertex.second.second = localDestVertexArray[node].second.second;
            }
        }
    }

    MPI_Barrier(MPI_COMM_WORLD);

    bool traceDone = false;

    if (procID == 0) {
        printf("The destination was Found\n");
        stack<Pair> Path;
        Pair tmpPair = make_pair(dest.first, dest.second);
        Path.push(tmpPair);
        tmpPair = make_pair(globalDestVertex.second.first, globalDestVertex.second.second);
	    
        while(tmpPair.first != src.first || tmpPair.second != src.second) {
            Path.push(tmpPair);

            int ownerNode = (tmpPair.first*dim_y + tmpPair.second) % nproc;
            if (ownerNode != procID) {
                MPI_Send((void*)(&tmpPair), 8, MPI_BYTE, ownerNode, TAG_REQ, MPI_COMM_WORLD);
                MPI_Recv((void*)(&tmpPair), 8, MPI_BYTE, ownerNode, TAG_PAIR, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            }
            else {
                tmpPair = make_pair(cellDetails[tmpPair.first*dim_y + tmpPair.second].parent_i, cellDetails[tmpPair.first*dim_y + tmpPair.second].parent_j);
            }
        }
        traceDone=true;
        Path.push(tmpPair);

        for (int node = 0; node < nproc; node++) {
            if (node != procID) {
                MPI_Isend((void*)(&traceDone), 1, MPI_BYTE, node, TAG_BOOL, MPI_COMM_WORLD, &request);
            }
        }

        printf("\nThe Path is ");
        while (!Path.empty()) {
		    pair<int, int> p = Path.top();
		    Path.pop();
		    printf("-> (%d,%d) ", p.first, p.second);
	    }
        printf("\n");
    }

    while (!traceDone) {
        int flag;
        int flag_done;
        Pair tmpPair;
        Pair parentPair;
        MPI_Iprobe(0, TAG_BOOL, MPI_COMM_WORLD, &flag_done, MPI_STATUS_IGNORE);
        if (flag_done) {
            MPI_Irecv((void*)(&traceDone), 1, MPI_BYTE, 0, TAG_BOOL, MPI_COMM_WORLD, &request);
        }

        MPI_Iprobe(0, TAG_REQ, MPI_COMM_WORLD, &flag, MPI_STATUS_IGNORE);
		if (flag) {
			MPI_Recv((void*)(&tmpPair), 8, MPI_BYTE, 0, TAG_REQ, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            parentPair = make_pair(cellDetails[tmpPair.first*dim_y + tmpPair.second].parent_i, cellDetails[tmpPair.first*dim_y + tmpPair.second].parent_j);
            MPI_Send((void*)(&parentPair), 8, MPI_BYTE, 0, TAG_PAIR, MPI_COMM_WORLD);
        }
    }

    MPI_Barrier(MPI_COMM_WORLD);

    for (int node = 0; node < nproc; node++)
    {
        pPair tmpPair;
        cell tmpCell;
        pair<bool, int> tmp;
        int flag_data = 0;
        int flag_done = 0;
        int flag_cell = 0;
        if (node != procID)
        {	
            MPI_Iprobe(node, TAG_DATA, MPI_COMM_WORLD, &flag_data, MPI_STATUS_IGNORE);
            MPI_Iprobe(node, TAG_CELL, MPI_COMM_WORLD, &flag_cell, MPI_STATUS_IGNORE);
            while (flag_data && flag_cell)
            {
                MPI_Irecv((void*)(&tmpPair), 16, MPI_BYTE, node, TAG_DATA, MPI_COMM_WORLD, &request);
                MPI_Wait(&request, MPI_STATUS_IGNORE);
                MPI_Irecv((void*)(&tmpCell), 32, MPI_BYTE, node, TAG_CELL, MPI_COMM_WORLD, &request);
                MPI_Wait(&request, MPI_STATUS_IGNORE);

                MPI_Iprobe(node, TAG_DATA, MPI_COMM_WORLD, &flag_data, MPI_STATUS_IGNORE);
                MPI_Iprobe(node, TAG_CELL, MPI_COMM_WORLD, &flag_cell, MPI_STATUS_IGNORE);
            }
            
            MPI_Iprobe(node, TAG_DONE, MPI_COMM_WORLD, &flag_done, MPI_STATUS_IGNORE);
            if (flag_done){
                MPI_Irecv((void*)(&tmp), 8, MPI_BYTE, node, TAG_DONE, MPI_COMM_WORLD, &request);
                MPI_Wait(&request, MPI_STATUS_IGNORE);
            }
        }
    }

	return computeTime;
}


int main(int argc, char *argv[]) {
    char *inputFilename = NULL;
    int numProc;
    int opt = 0;

    MPI_Init(&argc, &argv);      // initialize MPI environment
    
    // Read command line arguments
    do {
        opt = getopt(argc, argv, "f:");
        switch (opt) {
        case 'f':
            inputFilename = optarg;
            break;

        case -1:
            break;

        default:
            break;
        }
    } while (opt != -1);

    if (inputFilename == NULL) {
        printf("Usage: %s -f <filename>\n", argv[0]);
        MPI_Finalize();
        return -1;
    }
    
    // I/O Read in Map Traces
    FILE *input;

    input = fopen(inputFilename, "r");

    if (!input) {
        printf("Unable to open file: %s.\n", inputFilename);
        return -1;
    }

    int dim_x, dim_y;
    fscanf(input, "%d %d\n", &dim_x, &dim_y);

    int *map = (int*)calloc(dim_x*dim_y, sizeof(int));
    for (int i = 0; i < dim_x; i++) {
        for (int j = 0; j < dim_y; j++) {
            fscanf(input, "%d ", &map[i*dim_y + j]);
		}
    }

    // Source is the left-most top-most corner
    Pair src = make_pair(0, 0);
    // Destination is the right-most bottom-most corner
    Pair dest = make_pair(dim_x - 1, dim_y - 1);

	int procID, nproc;

	// Get process rank
    MPI_Comm_rank(MPI_COMM_WORLD, &procID);

    // Get total number of processes specificed at start of run
    MPI_Comm_size(MPI_COMM_WORLD, &nproc);

    // A* Sequential Algorithm
    double computeTime;
	for (int i = 0; i < 1; i++) {
	    computeTime = aStarSearch(map, src, dest, dim_x, dim_y, procID, nproc, inputFilename);
	}

    // Cleanup
    MPI_Finalize();

	printf("Elapsed time for proc %d: %f\n", procID, computeTime);
 
    return (0);
}


