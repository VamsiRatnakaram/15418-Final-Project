#define SetMark(p,m) ((p)|(m))
#define ClearMark(p,m) ((p)&~(uintptr_t)(m))
#define IsMarked(p,m) ((p)&(uintptr_t)(m))
#define F_ADP 0x1U
#define F_PRG 0x2U
#define F_DEL 0x1U
#define F_ALL 0x3U
#define Clear(p) ClearMark(p,F_ALL)


class MDList{
    const int D;
    const int N;
    Node* head;
};

struct AdoptDesc{
    Node* curr;
    int dp,dc;
};

struct Node{
    atomic<uintptr_t> val;
    int key;
    int k[D];
    atomic<uintptr_t> child[D];
    atomic<AdoptDesc*> adesc;
};

struct HeadNode:Node{
    int ver;
};

struct Stack{
    Node* node[D];
    HeadNode* head;
}

class PriorityQueue {
    int N;
    int R;
    atomic_bool notPurging{true};
    atomic<int> nMarkedNodes{0};
    atomic<uintptr_t> head;
    atomic<Stack*> stack;
    HeadNode firstHeadNode;
    Stack firstStack;

    PriorityQueue(int N,int R): N(N),R(R){
        firstHeadNode.val=F_DEL;
        firstHeadNode.adesc=NULL;
        firstHeadNode.key=0;
        setCoords(&firstHeadNode,0);
        firstHeadNode.ver=1;
        for(int i=0;i<D;i++){
            firstStack.node[i]=&firstHeadNode;
        }
        stack.store(&firstStack);
    }
    
    void setCoords(Node* n,int key){
        int basis=ceil(pow(N,1.0/D));
        int quotient=key;
        int* k=n->k;
        for(int i=D-1;i>=0;i--){
            k[i]=quotient%basis;
            quotient=quotient/basis;
        }
    }
};



vector<int> keyToCoord(int key){
    int basis=ceil(pow(N,1.0/D));
    int quotient=key;
    vector<int> k;
    k.resize(D);
    for(int i=D-1;i>=0;i--){
        k[i]=quotient%basis;
        quotient=quotient/basis;
    }
    return k;
}

Node* searchNode(vector<int> k){
    Node* cur=head;
    int d=0;
    while(d<D){
        while(cur != NULL && k[D] > cur->k[d]){
            cur=cur->child[d];
        }
        if(cur==NULL || k[d] < cur->k[D])
            return NULL;
        d++;
    }
    return cur;
}