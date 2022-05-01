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
#include "mpi.h"

#define BILLION  1E9

#define TAG_CELL 3
#define TAG_DATA 2
#define TAG_DONE 1

// #include "CycleTimer.h"
using namespace std::chrono;

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

// A Utility Function to trace the path from the source
// to destination
void tracePath(cell *cellDetails, Pair dest, int dim_y)
{
	printf("\nThe Path is ");
	int row = dest.first;
	int col = dest.second;

	stack<Pair> Path;

	while (!(cellDetails[row*dim_y + col].parent_i == row
			&& cellDetails[row*dim_y + col].parent_j == col)) {
		Path.push(make_pair(row, col));
		int temp_row = cellDetails[row*dim_y + col].parent_i;
		int temp_col = cellDetails[row*dim_y + col].parent_j;
		row = temp_row;
		col = temp_col;
	}

	Path.push(make_pair(row, col));
	while (!Path.empty()) {
		pair<int, int> p = Path.top();
		Path.pop();
		printf("-> (%d,%d) ", p.first, p.second);
	}

	printf("\n");

	return;
}

void hashedSendFunction(pPair pairToSend, int dim_x, int nproc, int procID, std::priority_queue<pPair> *openList, cell cellToSend, cell *cellDetails) {
	int x = pairToSend.second.first;
	int y = pairToSend.second.second;
	int index = x*dim_x + y;

	int indexToSend = index % nproc;

	 MPI_Request request = MPI_REQUEST_NULL;

	if (indexToSend == procID) {
		openList->push(pairToSend);
		// Update the details of this cell
		cellDetails[x*dim_x+y] = cellToSend;
	}
	else {
		MPI_Isend((void*)(&pairToSend), 16, MPI_BYTE, indexToSend, TAG_DATA, MPI_COMM_WORLD, &request);
		MPI_Isend((void*)(&cellToSend), 32, MPI_BYTE, indexToSend, TAG_CELL, MPI_COMM_WORLD, &request);
		if (procID != 0) {
			MPI_Isend((void*)(&pairToSend), 16, MPI_BYTE, 0, TAG_DATA, MPI_COMM_WORLD, &request);
			MPI_Isend((void*)(&cellToSend), 32, MPI_BYTE, 0, TAG_CELL, MPI_COMM_WORLD, &request);
		}
	}
}

// A Function to find the shortest path between
// a given source cell to a destination cell according
// to A* Search Algorithm
void aStarSearch(int *map, Pair src, Pair dest, int dim_x, int dim_y, int procID, int nproc)
{
	// If the source is out of range
	if (isValid(src.first, src.second, dim_x, dim_y) == false) {
		printf("Source is invalid\n");
		return;
	}

	// If the destination is out of range
	if (isValid(dest.first, dest.second, dim_x, dim_y) == false) {
		printf("Destination is invalid\n");
		return;
	}

	// Either the source or the destination is blocked
	if (isUnBlocked(map, src.first, src.second, dim_y) == false
		|| isUnBlocked(map, dest.first, dest.second, dim_y)
			== false) {
		printf("Source or the destination is blocked\n");
		return;
	}

	// If the destination cell is the same as source cell
	if (isDestination(src.first, src.second, dest)
		== true) {
		printf("We are already at the destination\n");
		return;
	}

	// Create a closed list and initialise it to false which
	// means that no cell has been included yet This closed
	// list is implemented as a boolean 2D array
	bool closedList[dim_x][dim_y];
	memset(closedList, false, sizeof(closedList));

	// Declare a 2D array of structure to hold the details
	// of that cell
	cell cellDetails[dim_x][dim_y];

	int i, j;

	for (i = 0; i < dim_y; i++) {
		for (j = 0; j < dim_x; j++) {
			cellDetails[i][j].f = INT_MAX;
			cellDetails[i][j].g = INT_MAX;
			cellDetails[i][j].h = INT_MAX;
			cellDetails[i][j].parent_i = -1;
			cellDetails[i][j].parent_j = -1;
		}
	}

	/*
	Create an priority queue having information as-
	<f, <i, j>>
	where f = g + h,
	and i, j are the row and column index of that cell
	Note that 0 <= i <= ROW-1 & 0 <= j <= COL-1
	This open list is implemented as a set of pair of
	pair.*/
    std::priority_queue<pPair> openList;

	// Put the starting cell on the open list and set its
	// 'f' as 0
	if (procID == 0) {
		// Initialising the parameters of the starting node
		i = src.first, j = src.second;
		cellDetails[i][j].f = 0.0;
		cellDetails[i][j].g = 0.0;
		cellDetails[i][j].h = 0.0;
		cellDetails[i][j].parent_i = i;
		cellDetails[i][j].parent_j = j;
		openList.push(make_pair(0.0, make_pair(i, j)));
	}

	// We set this boolean value as false as initially
	// the destination is not reached.
	bool foundDest = false;
	bool localFoundDest = false;

	MPI_Request request = MPI_REQUEST_NULL;

	while (true) {

		//During each iteration, check for data messages from other nodes
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
					int tmpx = tmpPair.second.first;
					int tmpy = tmpPair.second.second;
					int tmpindex = tmpx*dim_x + tmpy;

					int tmpindexToSend = tmpindex % nproc;
					if (tmpindexToSend == procID && !closedList[tmpx][tmpy]) {
						openList.push(tmpPair);
					}
					MPI_Irecv((void*)(&tmpCell), 32, MPI_BYTE, node, TAG_CELL, MPI_COMM_WORLD, &request);
					if (!closedList[tmpx][tmpy]) {
						cellDetails[tmpPair.second.first][tmpPair.second.second] = tmpCell;
					}
					// printf("data being recieved");
					MPI_Iprobe(node, TAG_DATA, MPI_COMM_WORLD, &flag_data, MPI_STATUS_IGNORE);
					MPI_Iprobe(node, TAG_CELL, MPI_COMM_WORLD, &flag_cell, MPI_STATUS_IGNORE);
				}
				MPI_Iprobe(node, TAG_DONE, MPI_COMM_WORLD, &flag_done, MPI_STATUS_IGNORE);
				while (flag_done)
				{
					MPI_Irecv((void*)(&foundDest), 1, MPI_BYTE, node, TAG_DONE, MPI_COMM_WORLD, &request);
					// printf("recieved\n");
					MPI_Iprobe(node, TAG_DONE, MPI_COMM_WORLD, &flag_done, MPI_STATUS_IGNORE);
				}
			}
		}

		if (openList.size() == 0 && foundDest) {
			break;
		}

		if (openList.size() == 0) {
			continue;
		}

		pPair p = openList.top();
        openList.pop();

		i = p.second.first;
		j = p.second.second;
		if (closedList[i][j]) {
			continue;
		}

		// // Add this vertex to the closed list
		// i = p.second.first;
		// j = p.second.second;
		// closedList[i][j] = true;

		// printf("This is thread %d, I'm working on node %d, %d\n", procID, i, j);

		/*
		Generating all the 4 successor of this cell
		Cell-->Popped Cell (i, j)
		N --> North	 (i-1, j)
		S --> South	 (i+1, j)
		E --> East	 (i, j+1)
		W --> West   (i, j-1)*/

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
                    // Set the Parent of the destination cell
                    cellDetails[x][y].parent_i = i;
                    cellDetails[x][y].parent_j = j;
                    printf("The destination cell is found by thread %d\n", procID);
                    // tracePath((cell*)((void*)&cellDetails), dest, dim_x);
                    foundDest = true;
					localFoundDest = true;
					printf("past Trace\n");
					for (int node = 0; node < nproc; node++) {
						if (node != procID) {
							MPI_Isend((void*)(&foundDest), 1, MPI_BYTE, node, TAG_DONE, MPI_COMM_WORLD, &request);
						}
					}
                    continue;
                }
                else if (isUnBlocked(map, i, j, dim_x) == true) {
                    gNew = cellDetails[i][j].g + 1.0;
                    hNew = calculateHValue(x, y, dest);
                    fNew = gNew + hNew;

                    if (cellDetails[x][y].f == INT_MAX
                        || cellDetails[x][y].f > fNew) {

						pPair pairToSend = make_pair(fNew, make_pair(x, y));
						cell cellToSend = {i, j, fNew, gNew, hNew};
						// cellDetails[x][y] = cellToSend;
			
						hashedSendFunction(pairToSend, dim_y, nproc, procID, &openList, cellToSend, (cell*)(&cellDetails));
                    }
                }
            }

			// Add this vertex to the closed list
			closedList[i][j] = true;
        }
	}

	printf("hi thread %d is at the barrier\n", procID);


	// /* Create a datatype to receive into. */
	// MPI_Type_vector( NUM_LOCAL_ELE, /* # of blocks */
	// 1, /* # of datatypes in a block (one for this
	// array) */
	// gsize, /* Stride between successive blocks */
	// MPI_CHAR, /* Type of each block */
	// &old_type);
	// MPI_Type_commit( &old_type);

	// /* Resize the type to allow interleaving,
	// * so make it only one MPI_CHAR wide
	// */
	// MPI_Type_create_resized(old_type,
	// 0, /* Lower Bound */
	// 1, /* Uppoer Bound change to one block */
	// &new_type);
	// MPI_Type_commit( &new_type);

	MPI_Barrier(MPI_COMM_WORLD);

	// //Node 0 should check for leftover messages
	// if (procID == 0) {
	// 	for (int node = 0; node < nproc; node++)
	// 	{
	// 		pPair tmpPair;
	// 		cell tmpCell;
	// 		int flag_data = 0;
	// 		int flag_done = 0;
	// 		int flag_cell = 0;
	// 		if (node != procID)
	// 		{	
	// 			MPI_Iprobe(node, TAG_DATA, MPI_COMM_WORLD, &flag_data, MPI_STATUS_IGNORE);
	// 			MPI_Iprobe(node, TAG_CELL, MPI_COMM_WORLD, &flag_cell, MPI_STATUS_IGNORE);
	// 			while (flag_data && flag_cell)
	// 			{
	// 				MPI_Irecv((void*)(&tmpPair), 16, MPI_BYTE, node, TAG_DATA, MPI_COMM_WORLD, &request);
	// 				MPI_Irecv((void*)(&tmpCell), 32, MPI_BYTE, node, TAG_CELL, MPI_COMM_WORLD, &request);
	// 				cellDetails[tmpPair.second.first][tmpPair.second.second] = tmpCell;
	// 				printf("data being recieved");
	// 				MPI_Iprobe(node, TAG_DATA, MPI_COMM_WORLD, &flag_data, MPI_STATUS_IGNORE);
	// 				MPI_Iprobe(node, TAG_CELL, MPI_COMM_WORLD, &flag_cell, MPI_STATUS_IGNORE);
	// 			}
	// 		}
	// 	}
	// }

	if (localFoundDest) {
		MPI_Isend((void*)(&cellDetails[dest.first][dest.second]), 32, MPI_BYTE, 0, TAG_CELL, MPI_COMM_WORLD, &request);
	}

	MPI_Barrier(MPI_COMM_WORLD);

	// Put the starting cell on the open list and set its
	// 'f' as 0
	// if (procID == 0) {
	// 	// Initialising the parameters of the starting node
	// 	i = src.first, j = src.second;
	// 	cellDetails[i][j].parent_i = i;
	// 	cellDetails[i][j].parent_j = j;
	// }

	if (procID == 0) {
		for (int node = 0; node < nproc; node++) {
			int flag = 0;
			MPI_Iprobe(node, TAG_CELL, MPI_COMM_WORLD, &flag, MPI_STATUS_IGNORE);
			if (flag) {
				MPI_Irecv((void*)(&cellDetails[dest.first][dest.second]), 32, MPI_BYTE, node, TAG_CELL, MPI_COMM_WORLD, &request);
			}
		}
		// cellDetails[tmpPair.second.first][tmpPair.second.second] = tmpCell;
	}

	MPI_Barrier(MPI_COMM_WORLD);

	if (procID == 0) {
		for (int a = 0; a < dim_y; a++) {
			for (int b = 0; b < dim_x; b++) {
				printf("%d %d \n", cellDetails[a][b].parent_i, cellDetails[a][b].parent_j);
			}
		}
	}

	if (procID == 0) {
		tracePath((cell*)((void*)&cellDetails), dest, dim_x);
	}

	// Call TracePath here from root and use a gather to collect the celldetails from all cells
	// https://beowulf.beowulf.narkive.com/e6tVmOqX/mpi-programming-question-interleaved-mpi-gatherv

	// When the destination cell is not found and the open
	// list is empty, then we conclude that we failed to
	// reach the destination cell. This may happen when the
	// there is no way to destination cell (due to
	// blockages)
	if (foundDest == false)
		printf("Failed to find the Destination Cell\n");

	return;
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
	// Calculate time taken by a request
	// struct timespec requestStart, requestEnd;
	// clock_gettime(CLOCK_REALTIME, &requestStart);
    
    // StartTime after intialization
    double startTime = MPI_Wtime();
	for (int i = 0; i < 1; i++) {
		aStarSearch(map, src, dest, dim_x, dim_y, procID, nproc);
	}
	// clock_gettime(CLOCK_REALTIME, &requestEnd);
	// Calculate time it took
	// double accum = ( requestEnd.tv_sec - requestStart.tv_sec ) + ( requestEnd.tv_nsec - requestStart.tv_nsec )/ BILLION;

    // EndTime before I/O
    double endTime = MPI_Wtime();
	double computeTime = endTime - startTime;

    // Cleanup
    MPI_Finalize();

	printf("Elapsed time for proc %d: %f\n", procID, computeTime);
 
    return (0);
}


