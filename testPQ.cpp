#include <tbb/concurrent_priority_queue.h>
#include <bits/stdc++.h>
#include <vector>
#include <algorithm>

using namespace std;

// Creating a shortcut for int, int pair type
typedef pair<int, int> Pair;

// Creating a shortcut for pair<int, pair<int, int>> type
typedef pair<double, pair<int, int> > pPair;

struct compare{
    bool operator() (const pPair &p1,const pPair &p2 ){
         return p1.first>p2.first;
    }
};

int main() {
    concurrent_priority_queue<pPair, compare> tbbQ;

    pPair tmp = make_pair(1.0, make_pair(0.0,0.0));
    tbbQ.push(tmp);
    pPair tmp = make_pair(3.0, make_pair(0.0,0.0));
    tbbQ.push(tmp);
    pPair tmp = make_pair(2.0, make_pair(0.0,0.0));
    tbbQ.push(tmp);

    pPair tmp4;
    tbbQ.try_pop(&tmp4);
    pPair tmp5;
    tbbQ.try_pop(&tmp5);
    pPair tmp6;
    tbbQ.try_pop(&tmp6);

    printf("%f %f %f", tmp4.first, tmp5.first, tmp6.first);
    return 0;
}