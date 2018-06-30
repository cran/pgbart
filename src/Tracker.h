#ifndef H_TRACKER
#define H_TRACKER

#include "Node.h"

class State{
public:
	bool isgrow;
	int var;
	int split_idx;
        int LeftEx;
        int RightEx;
	State* next;
};

class Tracker{
public:
	State* front;
	State* rear;
	State* cur;

	Tracker();
	~Tracker();
	void append(bool misgrow, int mvar, int msplit_idx, int mLeftEx, int mRightEx);
	void reset();
	bool gonext();
	void CopyFrom(Tracker* src);
};

#endif
