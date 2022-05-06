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

#define BILLION  1E9;

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
	// return abs(dest.second - col) + abs(dest.first - row);
	return 0;
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
	int i=0;
	while (!Path.empty()) {
		pair<int, int> p = Path.top();
		Path.pop();
		printf("-> (%d,%d) ", p.first, p.second);
		i++;
	}
	printf("The path cost:%d\n",i);

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
    Output << "outputs//seq_" << rawname.c_str() << "_" << num_of_threads << ".txt";
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

	int pathLength=0;
	while (!Path.empty()) {
		pathLength++;
		pair<int, int> p = Path.top();
		Path.pop();
		fprintf(outFile, "%d %d \n", p.first, p.second);
	}

	fprintf(outFile, "LENGTH %d", pathLength);

	fclose(outFile);

	return;
}

// A Function to find the shortest path between
// a given source cell to a destination cell according
// to A* Search Algorithm
double aStarSearch(int *map, Pair src, Pair dest, int dim_x, int dim_y, const char *input_filename)
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
	openList.push(make_pair(0.0, make_pair(i, j)));

	// We set this boolean value as false as initially
	// the destination is not reached.
	bool foundDest = false;
	while (openList.size() != 0 && !foundDest) {

		pPair p = openList.top();
        openList.pop();

		// Add this vertex to the closed list
		i = p.second.first;
		j = p.second.second;
		if (closedList[i][j]) {
			continue;
		}

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
                    cellDetails[x*dim_y+y].parent_i = i;
                    cellDetails[x*dim_y+y].parent_j = j;
                    foundDest = true;
                    break;
                }
                else if (isUnBlocked(map, x, y, dim_y) == true) {
                    gNew = cellDetails[i*dim_y+j].g + 1.0;
                    hNew = calculateHValue(x, y, dest);
                    fNew = gNew + hNew;

                    if (cellDetails[x*dim_y+y].f == INT_MAX
					
                        || cellDetails[x*dim_y+y].f > fNew) {
                        openList.push(make_pair(
                            fNew, make_pair(x, y)));

                        // Update the details of this cell
                        cellDetails[x*dim_y+y].f = fNew;
                        cellDetails[x*dim_y+y].g = gNew;
                        cellDetails[x*dim_y+y].h = hNew;
                        cellDetails[x*dim_y+y].parent_i = i;
                        cellDetails[x*dim_y+y].parent_j = j;
                    }
                }
            }
        }

		closedList[i][j] = true;
	}
	clock_gettime(CLOCK_REALTIME, &requestEnd);
	printf("The destination cell is found\n");
	tracePath(cellDetails, dest, dim_y);
	// Calculate time it took
	double accum = ( requestEnd.tv_sec - requestStart.tv_sec ) + ( requestEnd.tv_nsec - requestStart.tv_nsec )/ BILLION;

	outputFile(input_filename, (bool*)closedList, cellDetails, dim_x, dim_y, 1, dest);

	// When the destination cell is not found and the open
	// list is empty, then we conclude that we failed to
	// reach the destination cell. This may happen when the
	// there is no way to destination cell (due to
	// blockages)
	if (foundDest == false)
		printf("Failed to find the Destination Cell\n");

	free(cellDetails);

	return accum;
}

int main(int argc, char *argv[]) {
    char *inputFilename = NULL;
    int numProc;
    int opt = 0;
    
    // Read command line arguments
    do {
        opt = getopt(argc, argv, "f:p");
        switch (opt) {
        case 'f':
            inputFilename = optarg;
            break;

        case 'p':
            numProc = atof(optarg);
            break;

        case -1:
            break;

        default:
            break;
        }
    } while (opt != -1);

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
	for (int i = 0; i < 1; i++) {
		accum = aStarSearch(map, src, dest, dim_x, dim_y, inputFilename);
	}

	printf("\n");
	printf("Total Execution Time: %lf\n", accum);
	printf("Size of cell %ld\n",sizeof(cell));
	printf("Size of pPair %ld\n",sizeof(pPair));
	printf("Size of Pair %ld\n",sizeof(Pair));

	fclose(input);
 
    return (0);
}


