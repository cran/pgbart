#include <cstdio>
#include <iostream>
#include <exception>
//#include <math.h>
#include <cmath>
//#include <stdlib.h>
#include <cstdlib>
//#include <time.h>
#include <ctime>

// extern "C" {
#include <R.h>
#include <Rmath.h>

#include "Queue.h"

Queue::Queue(){
    this->front = this->rear = new Cell;
    this->rear->after = NULL;
}

Queue::~Queue(){
    while(!(this->empty()))
        this->pop();
    delete this->front;
}

bool Queue::empty(){
    return this->front->after == NULL;
}

void Queue::append(void* p){
    Cell* cell = new Cell;
    cell->contents = p;
    cell->before = this->rear;
    cell->after = NULL;

    this->rear->after = cell;
    this->rear = cell;
}

void* Queue::pop(){
    if(this->empty())
        return 0;
    else{
        Cell* cur = this->front->after;

        this->front->after = cur->after;

        if(cur == this->rear){
            this->front->after = NULL;
            this->rear = this->front;
        }else{
            cur->after->before = this->front;
        }
        void* content = cur->contents;
        delete cur;
        return content;
    }
}

void Queue::CopyFrom(Queue* src){
    while(!(this->empty()))
        this->pop();

    Cell* src_cell = src->front->after;
    while(src_cell){
        this->append(src_cell->contents);
        src_cell = src_cell->after;
    }
}

