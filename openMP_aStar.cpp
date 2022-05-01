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
#include <omp.h>

#define BILLION  1E9;

// #include "CycleTimer.h"
using namespace std::chrono;

using namespace std;

static int _argc;
static const char **_argv;

const char *get_option_string(const char *option_name, const char *default_value) {
    for (int i = _argc - 2; i >= 0; i -= 2)
        if (strcmp(_argv[i], option_name) == 0)
            return _argv[i + 1];
    return default_value;
}

int get_option_int(const char *option_name, int default_value) {
    for (int i = _argc - 2; i >= 0; i -= 2)
        if (strcmp(_argv[i], option_name) == 0)
            return atoi(_argv[i + 1]);
    return default_value;
}

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

	return;
}

// A Function to find the shortest path between
// a given source cell to a destination cell according
// to A* Search Algorithm
void aStarSearch(int *map, Pair src, Pair dest, int dim_x, int dim_y)
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
	int l, m;

	for (l = 0; l < dim_y; l++) {
		for (m = 0; m < dim_x; m++) {
			cellDetails[l][m].f = INT_MAX;
			cellDetails[l][m].g = INT_MAX;
			cellDetails[l][m].h = INT_MAX;
			cellDetails[l][m].parent_i = -1;
			cellDetails[l][m].parent_j = -1;
		}
	}

	// Initialising the parameters of the starting node
	l = src.first, m = src.second;
	cellDetails[l][m].f = 0.0;
	cellDetails[l][m].g = 0.0;
	cellDetails[l][m].h = 0.0;
	cellDetails[l][m].parent_i = l;
	cellDetails[l][m].parent_j = m;

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
	openList.push(make_pair(0.0, make_pair(l, m)));

	// We set this boolean value as false as initially
	// the destination is not reached.
	bool foundDest = false;

	omp_lock_t openListLock, closeListLock,cellLock;
	omp_init_lock(&closeListLock); //cas
	omp_init_lock(&openListLock); //concurrent PQ
	omp_init_lock(&cellLock); //cas 

	#pragma omp parallel
	{
		int e=0;
		// printf("Num threads:%d \n",omp_get_num_threads());
		while (true) {
			if (foundDest && openList.size() == 0) {
				break;
			}
			pPair p;
			int i, j, x, y, z;
			// To store the 'g', 'h' and 'f' of the 4 successors
			double gNew, hNew, fNew;
			// printf("here\n");
			omp_set_lock(&openListLock);
			// printf("here\n");
			if(openList.size() == 0){
				omp_unset_lock(&openListLock);
				continue;
			}
			p = openList.top();
			openList.pop();
			omp_unset_lock(&openListLock);
			// Add this vertex to the closed list
			i = p.second.first;
			j = p.second.second;

			omp_set_lock(&closeListLock);
			if(closedList[i][j]){
				omp_unset_lock(&closeListLock);
				continue;
			}
			closedList[i][j] = true;
			e++;
			omp_unset_lock(&closeListLock);
			// printf("Thread worked on %d, Elements worked on:%d \n",omp_get_thread_num(),e);
			/*	
			Generating all the 4 successor of this cell
			Cell-->Popped Cell (i, j)
			N --> North	 (i-1, j)
			S --> South	 (i+1, j)
			E --> East	 (i, j+1)
			W --> West   (i, j-1)*/
			bool notDone=true;
			for (z = 0; z < 4; z++) {
				// Only process this cell if this is a valid one
				// int x, y;
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
						foundDest = true;
						notDone=false;
					}
					omp_set_lock(&openListLock);
					omp_set_lock(&closeListLock);
					omp_set_lock(&cellLock);
					if (closedList[x][y] == false
							&& isUnBlocked(map, i, j, dim_y)
									== true && notDone) {
						gNew = cellDetails[i][j].g + 1.0;
						hNew = calculateHValue(x, y, dest);
						fNew = gNew + hNew;

						if (cellDetails[x][y].f == INT_MAX
							|| cellDetails[x][y].f > fNew) {
							
							openList.push(make_pair(fNew, make_pair(x, y)));

							// Update the details of this cell
							cellDetails[x][y].f = fNew;
							cellDetails[x][y].g = gNew;
							cellDetails[x][y].h = hNew;
							cellDetails[x][y].parent_i = i;
							cellDetails[x][y].parent_j = j;
						}
					}
					omp_unset_lock(&cellLock);
					omp_unset_lock(&closeListLock);
					omp_unset_lock(&openListLock);
				}
			}
		}
		printf("Thread num:%d, elements worked on:%d\n",omp_get_thread_num(),e);
	}

	printf("The destination cell is found\n");	
	tracePath((cell*)((void*)&cellDetails), dest, dim_y);

	// When the destination cell is not found and the open
	// list is empty, then we conclude that we failed to
	// reach the destination cell. This may happen when the
	// there is no way to destination cell (due to
	// blockages)
	if (foundDest == false)
		printf("Failed to find the Destination Cell\n");

	return;
}

int main(int argc, const char *argv[]) {
    int numProc;
    // int opt = 0;

	_argc = argc - 1;
    _argv = argv + 1;

	const char *inputFilename = get_option_string("-f", NULL);
    numProc = get_option_int("-p", 1);
    
    // // Read command line arguments
    // do {
    //     opt = getopt(argc, argv, "f:p");
    //     switch (opt) {
    //     case 'f':
    //         inputFilename = optarg;
    //         break;

    //     case 'p':
    //         numProc = atoi(optarg);
    //         break;

    //     case -1:
    //         break;

    //     default:
    //         break;
    //     }
    // } while (opt != -1);

    if (inputFilename == NULL) {
        printf("Usage: %s -f <filename> [-p <P>] [-i <N_iters>]\n", argv[0]);
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

    // A* Sequential Algorithm
	// Calculate time taken by a request
	struct timespec requestStart, requestEnd;
	clock_gettime(CLOCK_REALTIME, &requestStart);
	omp_set_num_threads(numProc);
	for (int i = 0; i < 1; i++) {
		aStarSearch(map, src, dest, dim_x, dim_y);
	}
	clock_gettime(CLOCK_REALTIME, &requestEnd);
	// Calculate time it took
	double accum = ( requestEnd.tv_sec - requestStart.tv_sec ) + ( requestEnd.tv_nsec - requestStart.tv_nsec )/ BILLION;

	printf("\n");
	printf("Total Execution Time: %lf\n", accum);
 
    return (0);
}


