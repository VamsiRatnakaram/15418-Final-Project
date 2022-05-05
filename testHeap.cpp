#include <iostream>
#include <deque>
#include <vector>
#include <mutex>
#include <shared_mutex>
#include <vector>
#include <algorithm>
#include <bits/stdc++.h>
#include "minheap.hpp"

// Creating a shortcut for int, int pair type
typedef pair<int, int> Pair;

// Creating a shortcut for pair<int, pair<int, int>> type
typedef pair<double, pair<int, int> > pPair;

int main(){
    MinHeap heap;

    pPair tmp = make_pair(4.0, make_pair(0, 0));
    heap.insert(tmp);
    pPair tmp2 = make_pair(3.0, make_pair(0, 0));
    heap.insert(tmp2);
    pPair tmp3 = make_pair(6.0, make_pair(0, 0));
    heap.insert(tmp3);

    pPair tmp4 = heap.remove();

    heap.print();
    printf("%lf \n", tmp4.first);
    return 0;
}