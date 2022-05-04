#include <bits/stdc++.h>
#include <atomic>
#include <stdint.h>

#define SetMark(p,m) ((p)|(m))
#define ClearMark(p,m) ((p)&~(uintptr_t)(m))
#define IsMarked(p,m) ((p)&(uintptr_t)(m))
#define F_ADP 0x1U
#define F_PRG 0x2U
#define F_DEL 0x1U
#define F_ALL 0x3U
#define Clear(p) ClearMark(p,F_ALL)
const int D = 4;
using namespace std;

struct AdoptDesc;

struct Node{
    std::atomic <uintptr_t> val;
    int key;
    int k[D];
    std::atomic <uintptr_t> child[D];
    std::atomic <AdoptDesc*> adesc;
};

struct AdoptDesc{
    Node* curr;
    int dp,dc;
};

struct HeadNode:Node{
    int ver;
};

struct Stack{
    Node* node[D];
    HeadNode* head;
};

Node *newNode(){
    Node* temp=(Node*)calloc(1,sizeof(Node));
    return temp;
}

AdoptDesc *newDesc(){
    return (AdoptDesc*)calloc(1,sizeof(AdoptDesc));
}

Stack *newStack(){
    return (Stack*)calloc(1,sizeof(Stack));
}

HeadNode *newHeadNode(){
    return (HeadNode*)calloc(1,sizeof(HeadNode));
}

void setDupNode(Node *n, Node *o) {
    memcpy((void*)n,(void*)o,sizeof(Node));
}

class PriorityQueue {
    int N;
    int R;
    std::atomic_bool notPurging{true};
    std::atomic <int> nMarkedNodes{0};
    std::atomic <uintptr_t> head;
    std::atomic <Stack*> stack;
    HeadNode firstHeadNode;
    Stack firstStack;
    public:
    PriorityQueue(int N,int R): N(N),R(R){
        firstHeadNode.val=F_DEL;
        firstHeadNode.adesc=NULL;
        firstHeadNode.key=0;
        setCoords(&firstHeadNode,0);
        firstHeadNode.ver=1;
        for(int i=0;i<D;i++){
            firstHeadNode.child[i].store((uintptr_t)NULL);
        }
        head.store((uintptr_t)(&firstHeadNode));
        firstStack.head=&firstHeadNode;
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
        uintptr_t tempCur=head;
        Node* cur = (Node*)(void*)tempCur;
        int d=0;
        while(d<D){
            while(cur != NULL && k[D] > cur->k[d]){
                tempCur=cur->child[d];
                cur=(Node*)(void*)tempCur;
            }
            if(cur==NULL || k[d] < cur->k[D])
                return NULL;
            d++;
        }
        return cur;
    }

    bool insert(int key, uintptr_t val){
        //printf("Reached here\n");
        Stack *s = newStack();
        Node *n = newNode();
        n->key = key;
        n->val = val;
        setCoords(n,key);
        //printf("Reached here 2\n");
        for (int i = 0; i<D; i++) n->child[i].store((uintptr_t)NULL);
        //printf("Reached here 3\n");
        while (true) {
            Node* pred = NULL;
            int dp = 0, dc = 0;
            s->head = (HeadNode*) (head.load());
            Node* curr = s->head;
            //printf("Reached here 4\n");
            LocatePlace(dp, dc, pred, curr, n, s);
            //printf("Reached here 5\n");
            if (dc == D) {
                // this key is already in the queue
                return false;
            }
            finishInserting(curr, dp, dc);
            //printf("Reached here 6\n");
            FillNewNode(dp, dc,n, curr);
            //printf("Reached here 7\n");
            uintptr_t temp = (uintptr_t)curr;
            if (pred->child[dp].compare_exchange_strong(temp, (uintptr_t)n)) {
                finishInserting(n, dp, dc);
                RewindStack(s, n, pred, dp);
                return true;
            }
        }
    }

    inline void LocatePlace(int &dp, int &dc, Node *&pred, Node *&curr, Node *n, Stack *s) {
        while (dc < D) {
            while (curr != NULL && n->k[dc] > curr->k[dc]) {
                pred = curr;
                dp = dc;
                finishInserting(curr, dc, dc);
                curr = (Node*)(Clear(curr->child[dc].load()));
            }
            if (curr == NULL || n->k[dc] < curr->k[dc]) {
                break;
            }
            s->node[dc] = curr;
            dc++;
        }
    }

    inline void FillNewNode(int dp, int dc, Node *n, Node *curr) {
        if (dp < dc) {
            AdoptDesc *desc = newDesc();
            desc->curr = curr;
            desc->dc = dc;
            desc->dp = dp;
            n->adesc.store(desc);
        }
        else {
            n->adesc.store(NULL);
        }
        for (int i = 0; i < dp; i++) n->child[i] = F_ADP;
        for (int i = dp; i < D; i++) n->child[i] = (uintptr_t)NULL;
        n->child[dc] = (uintptr_t) curr;
    }

    void finishInserting(Node *n, int dp, int dc) {
        if (n == NULL) return;
        AdoptDesc *ad = n->adesc;
        if (ad == NULL || dc < ad->dp || dp > ad->dc) return;
        uintptr_t child;
        Node* curr = ad->curr;
        for (int i = ad->dp; i < ad->dc; i++) {
            child = Clear(curr->child[i].fetch_or(F_ADP));
            uintptr_t temp = (uintptr_t)NULL;
            n->child[i].compare_exchange_strong(temp,child);
        }
        n->adesc = NULL;
    }

    uintptr_t deleteMin() {
        Stack *sOld = stack.load();
        Stack *s = newStack();
        *s = *sOld;
        int d = D-1;
        //printf("Reached here 1\n");
        while (true) {
            Node *last = s->node[d];
            finishInserting(last, d, d);
            printf("Reached here 1\n");
            Node *child = (Node*) (Clear(last->child[d].load()));
            if (child==NULL) {
                if (d == 0) return (uintptr_t)NULL;
                d--;
                continue;
            }
            printf("Reached here 2\n");
            uintptr_t val = child->val;
            if (IsMarked(val, F_DEL)) {
                if (Clear(val) == (uintptr_t)NULL) {
                    for (int i = 0; i <D; i++) s->node[i] = child;
                }
                else {
                    s->head = (HeadNode*) (Clear(val));
                    for (int i = 0; i< D; i++) s->node[i] = s->head;
                }
                d = D-1;
            }
            else {
                printf("Reached here 3\n");
                if (child->val.compare_exchange_strong(val, F_DEL)) {
                    for (int i = d; i<D; i++) s->node[i] = child;
                    stack.compare_exchange_strong(sOld, s);
                    int marked = nMarkedNodes.fetch_add(1);
                    if (marked > R) purge(s->head, s->node[D-1]);
                    return val;
                }
            }
            printf("Reached here 4\n");
        }
    }

    void purge(HeadNode *hn, Node *prg) {
        if (!notPurging.load()) return;
        bool temp = true;
        if (!notPurging.compare_exchange_strong(temp, false)) return;
        if ((uintptr_t) (hn) != head.load()) {
            notPurging.store(true);
            return;
        }
        nMarkedNodes.store(0);
        HeadNode *hnNew = newHeadNode();
        Node *prgNew = newNode();
        setDupNode(prgNew,prg);
        
        hnNew->val = F_DEL;
        hnNew->ver = hn->ver + 1;
        hnNew->key = hn->key;
        setCoords(hnNew, 0);
        for (int i = 0; i<D; i++) hnNew->child[i].store((uintptr_t)NULL);
        int d = 0;
        Node *pvt = hn;
        uintptr_t child;
        while(d<D) {
            if (!LocatePivot(prg, pvt, d, child)) {
                pvt = hn;
                d = 0;
                continue;
            }
            if (hn == pvt) {
                hnNew->child[d].store(child);
                prgNew->child[d].store(child);
            }
            else {
                prgNew->child[d].store(child);
                if (d == 0 || prgNew->child[d-1].load() == F_ADP) {
                    hnNew->child[d].store((uintptr_t) prgNew);
                }
                else {
                    hnNew->child[d].store((uintptr_t)NULL);
                }
            }
            d++;
        }
        hn->val.store(SetMark((uintptr_t) prg, F_DEL));
        prg->val.store(SetMark((uintptr_t) hnNew, F_DEL));
        head.store((uintptr_t) hnNew);
        Stack *s = newStack();
        UpdateStackAfterPurge(s, hnNew);
        notPurging.store(true);
        return;
    }

    inline bool LocatePivot(Node *prg, Node *&pvt, int d, uintptr_t &child) {
        while (pvt->k[d] < prg->k[d]) {
            finishInserting(pvt, d, d);
            pvt = (Node*) (Clear(pvt->child[d])); 
        }
        do {
            child = pvt->child[d];
        } while (!IsMarked(child, F_ALL) && !pvt->child[d].compare_exchange_weak(child, SetMark(child, F_PRG)));
        if (IsMarked(child, F_ADP)) {
            return false;
        }
        else {
            child = ClearMark(child, F_PRG);
            return true;
        }
    }

    inline void RewindStack(Stack *s, Node *n, Node *pred, int dp) {
        // Note: no need to rewind stack if node is already deleted...
        for (bool first_iteration = true; !IsMarked(n->val, F_DEL); first_iteration = false) {
            Stack *sNow = stack.load();
            if (s->head->ver == sNow->head->ver) {
                if (n->key > sNow->node[D-1]->key) {
                    if (!first_iteration) break;
                    *s= *sNow;
                }
                else {
                    for (int i=dp; i<D; i++) s->node[i] = pred;
                }
            }
            else if (s->head->ver > sNow->head->ver) {
                Node *prg = (Node*) (ClearMark(sNow->head->val, F_DEL));
                if (prg->key <= sNow->node[D-1]->key) {
                    s->head = (HeadNode*) (ClearMark(prg->val, F_DEL));
                    for (size_t i=0; i <D; i++) s->node[i] = s->head;
                }
                else {
                    if (!first_iteration) break;
                    *s = *sNow;
                }
            }
            else { // s->head->ver < sNow->head->ver
                Node *prg = (Node*) (ClearMark(sNow->head->val, F_DEL));
                if (prg->key > n->key) {
                    for (int i =dp; i<D; i++) s->node[i] = pred;
                }
                else {
                    s->head = (HeadNode*) (ClearMark(prg->val, F_DEL));
                    for (int i=0; i<D; i++) s->node[i] = s->head;
                }
            }
            if (stack.compare_exchange_strong(sNow, s)) {
                break;
            }
        }
    }

    inline void UpdateStackAfterPurge(Stack *s, HeadNode *hnNew) {
        for (bool firstIteration = true; true; firstIteration = false) {
            Stack *sNow = stack.load();
            if (hnNew->ver <= sNow->head->ver) {
                // the stack has been updated already
                return;
            }

            Node *prg = (Node*) (ClearMark(sNow->head->val, F_DEL));
            if (prg->key <= sNow->node[D-1]->key) {
                s->head = (HeadNode*) (ClearMark(prg->val, F_DEL));
                for (size_t i=0; i<D; i++) s->node[i] = s->head;
            }
            else {
                if (!firstIteration) break;
                *s = *sNow;
            }
            if (stack.compare_exchange_strong(sNow, s)) {
                break;
            }
        }
    }
};

int main() {
    PriorityQueue pq(100,10);
    bool key=false;
    while(!key){
        key=pq.insert(10,10);
    }
    key=false;
    while(!key){
        key=pq.insert(5,5);
    }
    key=false;
    while(!key){
        key=pq.insert(20,20);
    }
    key=false;
    while(!key){
        key=pq.insert(30,30);
    }

    printf("Deleted Value:%ld\n",pq.deleteMin());
    // printf("Deleted Value:%ld\n",pq.deleteMin());
    // printf("Deleted Value:%ld\n",pq.deleteMin());
    // printf("Deleted Value:%ld\n",pq.deleteMin());
}