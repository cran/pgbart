#ifndef GUARD_Queue
#define GUARD_Queue

#include "List.h"

class Queue {
public:
    Cell *front;
    Cell *rear;
    Queue();
    ~Queue();

    bool empty();
    void append(void *p);
    void* pop();

    void CopyFrom(Queue* src);
};

#endif
