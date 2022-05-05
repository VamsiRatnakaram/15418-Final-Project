#include "minheap.hpp"

int main(){
    MinHeap heap;

    pPair tmp = make_pair(1.0, make_pair(0, 0));
    heap.insert(&tmp);
    pPair tmp2 = make_pair(4.0, make_pair(0, 0));
    heap.insert(&tmp2);
    pPair tmp3 = make_pair(6.0, make_pair(0, 0));
    heap.insert(&tmp3);

    heap.print();
    return;
}