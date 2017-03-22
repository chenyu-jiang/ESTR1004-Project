#pragma once

typedef struct Descripter {
	int x;
	int y;
	int direc;
	int vector[128];
} descripter;

typedef struct Linkedlist {
	descripter vector;
	struct Linkedlist* next;
} linkedlist;

linkedlist NewLinkedlist(descripter v);
