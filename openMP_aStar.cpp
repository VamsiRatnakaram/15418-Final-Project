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

#include "tmp/oneTBB/include/oneapi/tbb/concurrent_priority_queue.h"

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

void outputFile(const char *input_filename, bool *closedList, cell *cellDetails, int dim_x, int dim_y, int num_of_threads, Pair dest) {
	// File outputs
	char resolved_path[PATH_MAX];
    realpath(input_filename, resolved_path);
    char *base = basename(resolved_path);
    std::string baseS = std::string(base);
    size_t lastindex = baseS.find_last_of("."); 
    string rawname = baseS.substr(0, lastindex); 

    std::stringstream Output;
    Output << "outputs//openMP	_" << rawname.c_str() << "_" << num_of_threads << ".txt";
    std::string OutputFile = Output.str();

    const char *ocf = OutputFile.c_str();
    FILE *outFile = fopen(ocf, "w+");
    if (!outFile) {
        printf("Unable to open file: %s.\n", ocf);
        return;
    }

	fprintf(outFile, "%d %d \n", dim_x, dim_y);

	for (int i = 0; i<dim_x; i++) {
		for (int j = 0; j<dim_y; j++) {
			fprintf(outFile, "%d ", closedList[i*dim_y +j]);
		}
		fprintf(outFile, "\n");
	}
	
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
		fprintf(outFile, "%d %d \n", p.first, p.second);
	}

	fclose(outFile);

	return;
}

// A Function to find the shortest path between
// a given source cell to a destination cell according
// to A* Search Algorithm
double aStarSearch(int *map, Pair src, Pair dest, int dim_x, int dim_y, int nproc, const char* input_filename)
{
	// Calculate time taken by a request
	struct timespec requestStart, requestEnd;
	clock_gettime(CLOCK_REALTIME, &requestStart);

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

	// Create a closed list and initialise it to false which
	// means that no cell has been included yet This closed
	// list is implemented as a boolean 2D array
	bool closedList[dim_x][dim_y];
	memset(closedList, false, sizeof(closedList));

	// Declare a 2D array of structure to hold the details
	// of that cell
	// cell cellDetails[dim_x][dim_y];
	cell *cellDetails = (cell*)calloc(dim_x*dim_y, sizeof(cell));
	int l, m;

	for (l = 0; l < dim_y; l++) {
		for (m = 0; m < dim_x; m++) {
			cellDetails[l*dim_y+m].f = INT_MAX;
			cellDetails[l*dim_y+m].g = INT_MAX;
			cellDetails[l*dim_y+m].h = INT_MAX;
			cellDetails[l*dim_y+m].parent_i = -1;
			cellDetails[l*dim_y+m].parent_j = -1;
		}
	}

	// Initialising the parameters of the starting node
	l = src.first, m = src.second;
	cellDetails[l*dim_y+m].f = 0.0;
	cellDetails[l*dim_y+m].g = 0.0;
	cellDetails[l*dim_y+m].h = 0.0;
	cellDetails[l*dim_y+m].parent_i = l;
	cellDetails[l*dim_y+m].parent_j = m;

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
		// printf("Num threads:%d \n",omp_get_num_threads());
		while (!foundDest) {
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
			else  {
				// printf("size %d\n", openList.size());
				p = openList.top();
				openList.pop();
				i = p.second.first;
				j = p.second.second;
				omp_set_lock(&closeListLock);
				if(closedList[i][j]){
					omp_unset_lock(&closeListLock);
					omp_unset_lock(&openListLock);
					continue;
				}
				closedList[i][j]=true;
				omp_unset_lock(&closeListLock);
				omp_unset_lock(&openListLock);

				/*	
				Generating all the 4 successor of this cell
				Cell-->Popped Cell (i, j)
				N --> North	 (i-1, j)
				S --> South	 (i+1, j)
				E --> East	 (i, j+1)
				W --> West   (i, j-1)*/
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
							omp_set_lock(&cellLock);
							cellDetails[x*dim_y+y].parent_i = i;
							cellDetails[x*dim_y+y].parent_j = j;
							omp_unset_lock(&cellLock);
							foundDest = true;
						}
						if (isUnBlocked(map, x, y, dim_y) == true) {
							omp_set_lock(&openListLock);
							omp_set_lock(&cellLock);
							omp_set_lock(&closeListLock);
							gNew = cellDetails[i*dim_y+j].g + 1.0;
							hNew = calculateHValue(x, y, dest);
							fNew = gNew + hNew;

							if (cellDetails[x*dim_y+y].f == INT_MAX || cellDetails[x*dim_y+y].f > fNew) {
								
								if (!closedList[x][y]){
									openList.push(make_pair(fNew, make_pair(x, y)));
								} 

								// Update the details of this cell
								cellDetails[x*dim_y+y].f = fNew;
								cellDetails[x*dim_y+y].g = gNew;
								cellDetails[x*dim_y+y].h = hNew;
								cellDetails[x*dim_y+y].parent_i = i;
								cellDetails[x*dim_y+y].parent_j = j;
							}
							omp_unset_lock(&closeListLock);
							omp_unset_lock(&cellLock);
							omp_unset_lock(&openListLock);
						}
					}
				}
			}
		}
		// printf("Thread num:%d, elements worked on:%d\n",omp_get_thread_num(),e);
	}

	printf("The destination cell is found\n");	
	tracePath(cellDetails, dest, dim_y);

	clock_gettime(CLOCK_REALTIME, &requestEnd);
	// Calculate time it took
	double accum = ( requestEnd.tv_sec - requestStart.tv_sec ) + ( requestEnd.tv_nsec - requestStart.tv_nsec )/ BILLION;

	outputFile(input_filename, (bool*)closedList, cellDetails, dim_x, dim_y, nproc, dest);

	return accum;
}

int main(int argc, const char *argv[]) {
    int numProc;

	_argc = argc - 1;
    _argv = argv + 1;

	const char *inputFilename = get_option_string("-f", NULL);
    numProc = get_option_int("-p", 1);

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
	double accum;
	omp_set_num_threads(numProc);
	for (int i = 0; i < 1; i++) {
		accum = aStarSearch(map, src, dest, dim_x, dim_y, numProc, inputFilename);
	}

	printf("\n");
	printf("Total Execution Time: %lf\n", accum);
 
    return (0);
}