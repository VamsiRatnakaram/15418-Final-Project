class MDList{
    const int D;
    const int N;
    Node* head;
}

struct Node{
    int key, k[D];
    void* val;
    Node* child[D];
}