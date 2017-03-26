#pragma once

typedef struct Descripter {
	int x;
	int y;
	int direc;
	double vector[128];
} descripter;

typedef struct Point {
	int x;
	int y;
} point;


typedef struct PointPair {
	struct Point p1;
	struct Point p2;
} pointpair;