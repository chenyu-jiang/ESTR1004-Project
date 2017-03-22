#pragma once

typedef struct Descripter {
	int x;
	int y;
	int direc;
	double vector[32];
} descripter;

typedef struct Linkedlist {
	descripter vector;
	struct Linkedlist* next;
} linkedlist;

linkedlist NewLinkedlist(descripter v);
