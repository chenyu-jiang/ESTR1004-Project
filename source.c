#include <stdio.h>
#include "qdbmp.h"
#include <math.h>
#include <stdlib.h>
#include <string.h>
#include "Descriptor.h"
#include <time.h>

//#define DEBUG
// REGARDS TO http://qdbmp.sourceforge.net/

#define PI 3.14159
#define LOGE 2.71828
#define MINSECTION 10
#define POINTPAIR_THRHLD 0.65
#define RANSAC_RADIUS 25
#define RANSAC_TIMES 1000000

void RGBtoBW(BMP* bmp);
void invert(BMP* bmp);
BMP* zoom(BMP* bmp, int k);
BMP* downSample(BMP* bmp, int k);
void reduceLevel(BMP* bmp, int level);
void BritandCntr(BMP* bmp, int brightness, double contrast);
BMP* NaiveGaussianBlur(BMP* bmp, double theta, int radius);
BMP* FastGaussianBlur(BMP* bmp, double theta, int radius);
BMP* naiveRotate(double degrees, BMP* bmp);
BMP* shearRotate(double degrees, BMP* bmp);
BMP* shearRotateShell(double degrees, BMP* bmp);
BMP* SobelEdgeDetection(BMP* bmp, int thrhld, int coefficient);
BMP* AMRotation(BMP* bmp, double theta);
BMP* CannyEdgeDetection(BMP* bmp, int thrhld, int coefficient, double highShr, double lowShr, int radius);
BMP* DoubleThreshold(BMP* bmp, double high, double low, int radius);
int* EgdeConnection(int* img, int width, int height, int radius);
int* IsConnected(int* img, int width, int height, int x, int y, int radius);
BMP* DINTRotation(BMP* bmp, double theta);
int* NonMaximumSpr(double* img, int width, int height, int windowSize);
descripter* HarrisCorner(BMP* bmp, double thrhld, int* poic);
BMP* HarrisCornerDetector(BMP* bmp, double thrhld);
descripter* GenerateDescripter(int* Harris, int width, int height, int* GradY, int* GradX, BMP* bmp, int* GrA, int* point);
BMP* AMRotationParts(BMP* bmp, int x, int y, int radius, double theta, BMP** bmpt);
double distof128(descripter desc1, descripter desc2);
void DrawLine(BMP* bmp, int inix, int iniy, int desx, int desy, int r, int g, int b);
BMP* ImageStitching(BMP* img1, BMP* img2, double thrhld);
pointpair* GeneratePointpair(BMP* img1, BMP* img2, double thrhld, int *paircount);
double determinant(double a11, double a12, double a13, double a21, double a22, double a23, double a31, double a32, double a33);
affinematrix RANSACAffm(BMP* img1, BMP* img2, double thrhld, int times, double radius);
pointpair* RANSACMatch(BMP* img1, BMP* img2, double thrhld, int* RANSACcount);
int matchjudgement(affinematrix matrix, pointpair *naivepair, int paircount, int* match);//Judge how many point pairs has been matched under certain affine matrix;


int main()
{
	BMP *bmp, *bmp1;
	char filename[] = "C:\\Users\\HP\\Documents\\Visual Studio 2015\\Projects\\ImageProcessing\\Debug\\rm3.bmp";
	char filename1[] = "C:\\Users\\HP\\Documents\\Visual Studio 2015\\Projects\\ImageProcessing\\Debug\\rm4.bmp";
	bmp = BMP_ReadFile(filename);
	BMP_CHECK_ERROR(stderr, -1);
	bmp1 = BMP_ReadFile(filename1);
	BMP_CHECK_ERROR(stderr, -1);
	/////////////////////////////////////////////////////////////////////////
	//Your code in between
	srand(time(NULL));
	BMP* newbmp = ImageStitching(bmp, bmp1, 200000);
	char filepath[100] = "C:\\Users\\HP\\Documents\\Visual Studio 2015\\Projects\\ImageProcessing\\Debug\\Tst\\DOWS\\PointPair.bmp";
	BMP_WriteFile(newbmp, filepath);
	newbmp = HarrisCornerDetector(bmp, 200000);
	BMP_WriteFile(newbmp, "C:\\Users\\HP\\Documents\\Visual Studio 2015\\Projects\\ImageProcessing\\Debug\\Tst\\DOWS\\HarrisL.bmp");
	newbmp = HarrisCornerDetector(bmp1, 200000);
	BMP_WriteFile(newbmp, "C:\\Users\\HP\\Documents\\Visual Studio 2015\\Projects\\ImageProcessing\\Debug\\Tst\\DOWS\\HarrisR.bmp");
	BMP_Free(newbmp);
	/////////////////////////////////////////////////////////////////////////
	BMP_Free(bmp);
	BMP_Free(bmp1);
	//BMP_Free(output);
	return 0;
}

void RGBtoBW(BMP* bmp)
{
	UCHAR   r, g, b;
	UINT    width, height;
	UINT    x, y;
	width = BMP_GetWidth(bmp);
	height = BMP_GetHeight(bmp);
	for (x = 0; x < width; x++)
	{
		for (y = 0; y < height; y++)
		{
			BMP_GetPixelRGB(bmp, x, y, &r, &g, &b);
			BMP_SetPixelRGB(bmp, x, y, 0.3*r + 0.59*g + 0.11*b, 0.3*r + 0.59*g + 0.11*b, 0.3*r + 0.59*g + 0.11*b);
		}
	}
}

void invert(BMP* bmp)
{
	UCHAR   r, g, b;
	UINT    width, height;
	UINT    x, y;
	width = BMP_GetWidth(bmp);
	height = BMP_GetHeight(bmp);
	for (x = 0; x < width; x++)
	{
		for (y = 0; y < height; y++)
		{
			BMP_GetPixelRGB(bmp, x, y, &r, &g, &b);
			BMP_SetPixelRGB(bmp, x, y, 255 - r, 255 - g, 255 - b);
		}
	}
}

BMP* zoom(BMP* bmp, int k)
{
	UINT    width, height;
	width = BMP_GetWidth(bmp);
	height = BMP_GetHeight(bmp);
	BMP* newim = BMP_Create(width*k, height*k, BMP_GetDepth(bmp));
	for (UINT iX = 0; iX < width; iX++)
	{
		for (UINT iY = 0; iY < height; iY++)
		{
			UCHAR difr, difg, difb;
			UCHAR tmpr, tmpg, tmpb;
			BMP_GetPixelRGB(bmp, iX, iY, &tmpr, &tmpg, &tmpb);
			BMP_GetPixelRGB(bmp, iX, iY + 1, &difr, &difg, &difb);
			difr = -(tmpr - difr) / k;
			difg = -(tmpg - difg) / k;
			difb = -(tmpb - difb) / k;
			for (int i = 0; i < k; i++)
			{
				BMP_SetPixelRGB(newim, iX*k, iY*k + i, tmpr + difr*i, tmpg + difg*i, tmpb + difb*i);
			}
		}
	}
	for (UINT iY = 0; iY < height*k; iY++)
	{
		for (UINT iX = 0; iX < width; iX++)
		{
			UCHAR difr, difg, difb;
			UCHAR tmpr, tmpg, tmpb;
			BMP_GetPixelRGB(newim, iX*k, iY, &tmpr, &tmpg, &tmpb);
			BMP_GetPixelRGB(newim, iX*k + k, iY, &difr, &difg, &difb);
			difr = -(tmpr - difr) / k;
			difg = -(tmpg - difg) / k;
			difb = -(tmpb - difb) / k;
			for (int i = 0; i < k; i++)
			{
				BMP_SetPixelRGB(newim, iX*k + i, iY, tmpr + difr*i, tmpg + difg*i, tmpb + difb*i);
			}
		}
	}
	return newim;
}

BMP* downSample(BMP* bmp, int k)
{
	int width, height;
	int depth;
	width = BMP_GetWidth(bmp);
	height = BMP_GetHeight(bmp);
	depth = BMP_GetDepth(bmp);
	BMP* newim = BMP_Create(width / k, height / k, depth);
	for (int i= 0; i < width/k; i++)
	{
		for (int j = 0; j < height/k; j++)
		{
			int r = 0, g = 0, b = 0;
			BMP_GetPixelRGB(bmp, i*k, j*k, &r, &g, &b);
			BMP_SetPixelRGB(newim, i, j, r, g, b);
		}
	}
	return newim;
}

void reduceLevel(BMP* bmp, int level)
{
	UINT    width, height;
	USHORT depth;
	width = BMP_GetWidth(bmp);
	height = BMP_GetHeight(bmp);
	depth = BMP_GetDepth(bmp);
	for (UINT i = 0; i < width; i++)
	{
		for (UINT j = 0; j < height; j++)
		{
			UCHAR r, g, b;
			BMP_GetPixelRGB(bmp, i, j, &r, &g, &b);
			BMP_SetPixelRGB(bmp, i, j, (UCHAR)(r / (255 / level) / (1.0*level - 1) * 255), (UCHAR)(g / (255 / level) / (1.0*level - 1) * 255), (UCHAR)(b / (255 / level) / (1.0*level - 1) * 255));
			BMP_CHECK_ERROR(stdout, -1);
		}
	}
}

void BritandCntr(BMP* bmp, int brightness, double contrast) //Brightness: amount of RGB value you want to add
//Contrast: | 1 means original | <1 less contrast | >1 more contrast
{
	UINT    width, height;
	width = BMP_GetWidth(bmp);
	height = BMP_GetHeight(bmp);
	for (UINT i = 0; i < width; i++)
	{
		for (UINT j = 0; j < height; j++)
		{
			UCHAR Ur, Ug, Ub;
			BMP_GetPixelRGB(bmp, i, j, &Ur, &Ug, &Ub);
			int r, g, b;
			r = (int)(contrast*(Ur - 128)) + 128 + brightness;
			g = (int)(contrast*(Ug - 128)) + 128 + brightness;
			b = (int)(contrast*(Ub - 128)) + 128 + brightness;
			if (r < 0) r = 0;
			if (g < 0) g = 0;
			if (b < 0) b = 0;
			if (r > 255) r = 255;
			if (g > 255) g = 255;
			if (b > 255) b = 255;
			BMP_SetPixelRGB(bmp, i, j, (UCHAR)r, (UCHAR)g, (UCHAR)b);
			BMP_CHECK_ERROR(stdout, -1);
		}
	}
}

double GaussFunction(double theta, int dimension, double* gaussmatrix)
{
	double sum = 0;
	for (int i = 0; i < dimension; i++)
	{
		for (int j = 0; j < dimension; j++)
		{
			int x = i - (dimension - 1) / 2;
			int y = j - (dimension - 1) / 2;
			double gauss = 1.0 / (2 * PI*theta*theta)*pow(LOGE, 1.0*(-x*x - y*y) / (2 * theta*theta));
			gaussmatrix[i*dimension + j] = gauss;
			sum += gauss;
		}
	}
	return sum;
}

BMP* NaiveGaussianBlur(BMP* bmp, double theta, int radius)
{
	double* gaussmatrix = malloc(sizeof(double)*radius*radius);
	double sum = GaussFunction(theta, radius, gaussmatrix);
#ifdef DEBUG
	for (int i = 0; i < radius; i++)
	{
		for (int j = 0; j < radius; j++)
		{
			printf("%10.4lf", gaussmatrix[i*radius + j]);
		}
		printf("\n");
	}
#endif // DEBUG

	UINT    width, height;
	USHORT depth;
	width = BMP_GetWidth(bmp);
	height = BMP_GetHeight(bmp);
	depth = BMP_GetDepth(bmp);
	BMP* newbmp = BMP_Create(width, height, depth);
	//Go through all pixels
	for (UINT x = 0; x < width; x++)
	{
		for (UINT y = 0; y < height; y++)
		{
			//Fill in colors
			double r = 0, g = 0, b = 0;
			for (int i = 0; i < radius; i++)
			{
				for (int j = 0; j < radius; j++)
				{
					//Reflect boundary
					int xi = x, yi = y;
					if (xi - (radius - 1) / 2 + i < 0)
					{
						xi = -(xi - (radius - 1) / 2 + i);
					}
					else
					{
						xi = (xi - (radius - 1) / 2 + i);
					}
					if (yi - (radius - 1) / 2 + j < 0)
					{
						yi = -(yi - (radius - 1) / 2 + j);
					}
					else
					{
						yi = (yi - (radius - 1) / 2 + j);
					}
					if (xi >= width)
					{
						xi = width - (xi - width) - 1;
					}
					if (yi >= height)
					{
						yi = height - (yi - height) - 1;
					}
					UCHAR R, G, B;
					BMP_GetPixelRGB(bmp, xi, yi, &R, &G, &B);
					r += R*gaussmatrix[i*radius + j];
					g += G*gaussmatrix[i*radius + j];
					b += B*gaussmatrix[i*radius + j];
				}
			}

			//Calculate total
			r /= sum;
			g /= sum;
			b /= sum;
			BMP_SetPixelRGB(newbmp, x, y, (UCHAR)r, (UCHAR)g, (UCHAR)b);
		}
	}
	free(gaussmatrix);
	return newbmp;
}

BMP* FastGaussianBlur(BMP* bmp, double theta, int radius)
{
	double* gaussmatrix = malloc(sizeof(double)*radius*radius);
	double sum = GaussFunction(theta, radius, gaussmatrix);
	sum = 0;
	for (int i = 0; i < radius; i++)
	{
		sum += gaussmatrix[i + radius*(radius / 2)];
	}
#ifdef DEBUG
	for (int i = 0; i < radius; i++)
	{
		for (int j = 0; j < radius; j++)
		{
			printf("%10.4lf", gaussmatrix[i*radius + j]);
		}
		printf("\n");
	}
#endif // DEBUG

	UINT    width, height;
	USHORT depth;
	width = BMP_GetWidth(bmp);
	height = BMP_GetHeight(bmp);
	depth = BMP_GetDepth(bmp);
	BMP* newbmp = BMP_Create(width, height, depth);
	BMP* newbmp1 = BMP_Create(width, height, depth);
	for (UINT x = 0; x < width; x++)
	{
		for (UINT y = 0; y < height; y++)
		{
			//Fill in colors
			double r = 0, g = 0, b = 0;
			for (int i = 0; i < radius; i++)
			{
				//Reflect boundary
				int xi = x;
				if (xi - (radius - 1) / 2 + i < 0)
				{
					xi = -(xi - (radius - 1) / 2 + i);
				}
				else
				{
					xi = (xi - (radius - 1) / 2 + i);
				}
				if (xi >= width)
				{
					xi = width - (xi - width) - 1;
				}
				UCHAR R, G, B;
				BMP_GetPixelRGB(bmp, xi, y, &R, &G, &B);
				r += R*gaussmatrix[i + radius*(radius / 2)];
				g += G*gaussmatrix[i + radius*(radius / 2)];
				b += B*gaussmatrix[i + radius*(radius / 2)];
			}
			//Calculate total
			r /= sum;
			g /= sum;
			b /= sum;
			BMP_SetPixelRGB(newbmp, x, y, (UCHAR)r, (UCHAR)g, (UCHAR)b);
		}
	}

	for (UINT x = 0; x < width; x++)
	{
		for (UINT y = 0; y < height; y++)
		{
			//Fill in colors
			double r = 0, g = 0, b = 0;
			for (int i = 0; i < radius; i++)
			{
				//Reflect boundary
				int yi = y;
				if (yi - (radius - 1) / 2 + i < 0)
				{
					yi = -(yi - (radius - 1) / 2 + i);
				}
				else
				{
					yi = (yi - (radius - 1) / 2 + i);
				}
				if (yi >= height)
				{
					yi = height - (yi - height) - 1;
				}
				UCHAR R, G, B;
				BMP_GetPixelRGB(newbmp, x, yi, &R, &G, &B);
				r += R*gaussmatrix[i + radius*(radius / 2)];
				g += G*gaussmatrix[i + radius*(radius / 2)];
				b += B*gaussmatrix[i + radius*(radius / 2)];
			}
			//Calculate total
			r /= sum;
			g /= sum;
			b /= sum;
			BMP_SetPixelRGB(newbmp1, x, y, (UCHAR)r, (UCHAR)g, (UCHAR)b);
		}
	}
	free(gaussmatrix);
#ifdef DEBUG
	BMP_WriteFile(newbmp, "C:\\Users\\HP\\Documents\\Visual Studio 2015\\Projects\\ImageProcessing\\Debug\\Output\\Test.bmp");
#endif // DEBUG

	BMP_Free(newbmp);
	return newbmp1;
}

double* GaussianKernelDouble(double* img, double theta, int radius,int width,int height)
{
	double* gaussmatrix = malloc(sizeof(double)*radius*radius);
	double sum = GaussFunction(theta, radius, gaussmatrix);
	sum = 0;
	for (int i = 0; i < radius; i++)
	{
		sum += gaussmatrix[i + radius*(radius / 2)];
	}
	double *ans = malloc(sizeof(double)*width*height);
	double *img1 = malloc(sizeof(double)*width*height);
	memset(img1, 0, sizeof(double)*width*height);
	for (int x = 0; x < width; x++)
	{
		for (int y = 0; y < height; y++)
		{
			//Fill in colors
			double val = 0;
			for (int i = 0; i < radius; i++)
			{
				//Reflect boundary
				int xi = x;
				if (xi - (radius - 1) / 2 + i < 0)
				{
					xi = -(xi - (radius - 1) / 2 + i);
				}
				else
				{
					xi = (xi - (radius - 1) / 2 + i);
				}
				if (xi >= width)
				{
					xi = width - (xi - width) - 1;
				}
				double value = img[xi + y*width];
				val += value*gaussmatrix[i + radius*(radius / 2)];
			}
			//Calculate total
			val /= sum;
			img1[x + y*width] = val;
		}
	}

	for (int x = 0; x < width; x++)
	{
		for (int y = 0; y < height; y++)
		{
			//Fill in colors
			double val = 0;
			for (int i = 0; i < radius; i++)
			{
				//Reflect boundary
				int yi = y;
				if (yi - (radius - 1) / 2 + i < 0)
				{
					yi = -(yi - (radius - 1) / 2 + i);
				}
				else
				{
					yi = (yi - (radius - 1) / 2 + i);
				}
				if (yi >= height)
				{
					yi = height - (yi - height) - 1;
				}
				double value = img1[x + yi*width];
				val += value*gaussmatrix[i + radius*(radius / 2)];
			}
			//Calculate total
			val /= sum;
			ans[x + y*width] = val;
		}
	}
	free(gaussmatrix);
	free(img);
	free(img1);
	return ans;
}

double* GaussianKernel(int* img, double theta, int radius, int width, int height)
{
	double* gaussmatrix = malloc(sizeof(double)*radius*radius);
	double sum = GaussFunction(theta, radius, gaussmatrix);
	sum = 0;
	for (int i = 0; i < radius; i++)
	{
		sum += gaussmatrix[i + radius*(radius / 2)];
	}
	double *ans = malloc(sizeof(double)*width*height);
	double *img1 = malloc(sizeof(double)*width*height);
	memset(img1, 0, sizeof(double)*width*height);
	for (int x = 0; x < width; x++)
	{
		for (int y = 0; y < height; y++)
		{
			//Fill in colors
			double val = 0;
			for (int i = 0; i < radius; i++)
			{
				//Reflect boundary
				int xi = x;
				if (xi - (radius - 1) / 2 + i < 0)
				{
					xi = -(xi - (radius - 1) / 2 + i);
				}
				else
				{
					xi = (xi - (radius - 1) / 2 + i);
				}
				if (xi >= width)
				{
					xi = width - (xi - width) - 1;
				}
				int value = img[xi + y*width];
				val += value*gaussmatrix[i + radius*(radius / 2)];
			}
			//Calculate total
			val /= sum;
			img1[x + y*width] = val;
		}
	}

	for (int x = 0; x < width; x++)
	{
		for (int y = 0; y < height; y++)
		{
			//Fill in colors
			double val = 0;
			for (int i = 0; i < radius; i++)
			{
				//Reflect boundary
				int yi = y;
				if (yi - (radius - 1) / 2 + i < 0)
				{
					yi = -(yi - (radius - 1) / 2 + i);
				}
				else
				{
					yi = (yi - (radius - 1) / 2 + i);
				}
				if (yi >= height)
				{
					yi = height - (yi - height) - 1;
				}
				int value = img1[x + yi*width];
				val += value*gaussmatrix[i + radius*(radius / 2)];
			}
			//Calculate total
			val /= sum;
			ans[x + y*width] = val;
		}
	}
	free(gaussmatrix);
	free(img);
	free(img1);
	return ans;
}

BMP* naiveRotate(double degrees, BMP* bmp)
{
	UINT width, height;
	USHORT depth;
	width = BMP_GetWidth(bmp);
	height = BMP_GetHeight(bmp);
	depth = BMP_GetDepth(bmp);
	UINT lenth = (UINT)(sqrt(width*width + height*height) + 0.5);
	BMP* newim = BMP_Create(lenth, lenth, depth);
	UINT origin = lenth / 2;
	UINT picX = origin - width / 2;
	UINT picY = origin - height / 2;
	UINT cntX = (UINT)(1.0*width / 2 * cos(degrees) + 1.0*height / 2 * sin(degrees));
	UINT cntY = (UINT)(-1.0*width / 2 * sin(degrees) + 1.0*height / 2 * cos(degrees));
	long int disX = origin - cntX;
	long int disY = origin - cntY;
	for (UINT i = 0; i < width; i++)
	{
		for (UINT j = 0; j < height; j++)
		{
			UCHAR r, g, b;
			BMP_GetPixelRGB(bmp, i, j, &r, &g, &b);
			BMP_SetPixelRGB(newim, disX + i*cos(degrees) + j*sin(degrees), disY - i*sin(degrees) + j*cos(degrees), r, g, b);
		}
	}
	return newim;
}

BMP* transpose(BMP* bmp)
{
	UINT width, height;
	USHORT depth;
	width = BMP_GetWidth(bmp);
	height = BMP_GetHeight(bmp);
	depth = BMP_GetDepth(bmp);
	BMP* ans = BMP_Create(height, width, depth);
	for (int i = 0; i < width; i++)
	{
		for (int j = 0; j < height; j++)
		{
			UCHAR r, g, b;
			BMP_GetPixelRGB(bmp, i, j, &r, &g, &b);
			BMP_SetPixelRGB(ans, j, width - i, r, g, b);
		}
	}
	return ans;
}

BMP* r180(BMP* bmp)
{
	UINT width, height;
	USHORT depth;
	width = BMP_GetWidth(bmp);
	height = BMP_GetHeight(bmp);
	depth = BMP_GetDepth(bmp);
	BMP* ans = BMP_Create(width, height, depth);
	for (int i = 0; i < width; i++)
	{
		for (int j = 0; j < height; j++)
		{
			UCHAR r, g, b;
			BMP_GetPixelRGB(bmp, i, j, &r, &g, &b);
			BMP_SetPixelRGB(ans, width - i - 1, height - j - 1, r, g, b);
		}
	}
	return ans;
}

BMP* shearRotate(double degrees, BMP* bmp)
{
	UINT width, height;
	USHORT depth;
	width = BMP_GetWidth(bmp);
	height = BMP_GetHeight(bmp);
	depth = BMP_GetDepth(bmp);
	int outWidth = (int)(tan(degrees / 2)*(height)) + width;
	int outHeight = (int)(sin(degrees)*(width)) + height;
	BMP* newim = BMP_Create(outWidth, outHeight, depth);
	UINT originX = outWidth / 2;
	UINT originY = outHeight / 2;
	UINT picX = originX - width / 2;
	UINT picY = originY - height / 2;
	//perform first X shift
	int disX1 = (int)(-tan(degrees / 2)*(height)) / 2;
	//initialize flag to 0
	int* flag = malloc(sizeof(int)*outHeight*outWidth);
	memset(flag, 0, sizeof(int)*outHeight*outWidth);
	for (UINT i = 0; i < height; i++)
	{
		double prshX = (-tan(degrees / 2));
		int shX = (int)(prshX*i);
		for (UINT j = 0; j < width; j++)
		{
			UCHAR r, g, b;
			BMP_GetPixelRGB(bmp, j, i, &r, &g, &b);
			BMP_SetPixelRGB(newim, (UINT)(picX + j - shX + disX1), (UINT)(picY + i), r, g, b);
			flag[(picX + j - shX + disX1) + (picY + i)*outWidth] = 1;
		}
	}
#ifdef DEBUG
	BMP* debug0 = BMP_Create(outWidth, outHeight, depth);
	for (int i = 0; i < outWidth; i++)
	{
		for (int j = 0; j < outHeight; j++)
		{
			if (flag[i + outWidth*j] == 1) BMP_SetPixelRGB(debug0, (UINT)(i), (UINT)(j), 100, 100, 100);
		}
	}
	BMP_WriteFile(debug0, "C:\\Users\\HP\\Documents\\Visual Studio 2015\\Projects\\ImageProcessing\\Debug\\Output\\debug0.bmp");
	char filepath1[100] = "C:\\Users\\HP\\Documents\\Visual Studio 2015\\Projects\\ImageProcessing\\Debug\\Output\\rotate1.bmp";
	BMP_WriteFile(newim, filepath1);
	BMP_Free(debug0);
#endif // DEBUG
	//perform Y shift
	int disY = (int)(sin(degrees)*(outWidth)) / 2;
	for (int i = 0; i < outWidth / 2; i++)
	{
		double preshY = sin(degrees);
		int shY = (int)(preshY*i);
		for (int j = outHeight - 1; j >= 0; j--)
		{
			UCHAR r, g, b;
			BMP_GetPixelRGB(newim, i, j, &r, &g, &b);
			if (flag[j*outWidth + i] == 1)
			{
				BMP_SetPixelRGB(newim, (UINT)(i), (UINT)(j - shY + disY), r, g, b);
				if ((outWidth*(j - shY + disY) + i) >= 0 && (outWidth*(j - shY + disY) + i) < outWidth*outHeight) flag[outWidth*(j - shY + disY) + i] = 2;
			}
		}
	}
	for (int i = outWidth / 2; i < outWidth; i++)
	{
		double preshY = sin(degrees);
		int shY = (int)(preshY*i);
		for (int j = 0; j < outHeight; j++)
		{
			UCHAR r, g, b;
			BMP_GetPixelRGB(newim, i, j, &r, &g, &b);
			if (flag[j*outWidth + i] == 1)
			{
				BMP_SetPixelRGB(newim, (UINT)(i), (UINT)(j - shY + disY), r, g, b);
				if ((outWidth*(j - shY + disY) + i) >= 0 && (outWidth*(j - shY + disY) + i) < outWidth*outHeight) flag[outWidth*(j - shY + disY) + i] = 2;
			}
		}
	}
#ifdef DEBUG
	BMP* debug1 = BMP_Create(outWidth, outHeight, depth);
#endif // DEBUG
	for (int i = 0; i < outWidth; i++)
	{
		for (int j = 0; j < outHeight; j++)
		{
			if (flag[outWidth*j + i] != 2)
			{
				BMP_SetPixelRGB(newim, (UINT)(i), (UINT)(j), 0, 0, 0);
			}
#ifdef DEBUG
			if (flag[outWidth*j + i] == 2)
			{
				BMP_SetPixelRGB(debug1, (UINT)i, (UINT)j, (UCHAR)100, (UCHAR)100, (UCHAR)100);
			}
#endif // DEBUG
		}
	}
#ifdef DEBUG
	BMP_WriteFile(debug1, "C:\\Users\\HP\\Documents\\Visual Studio 2015\\Projects\\ImageProcessing\\Debug\\Output\\debug1.bmp");
	BMP_Free(debug1);
	char filepath2[100] = "C:\\Users\\HP\\Documents\\Visual Studio 2015\\Projects\\ImageProcessing\\Debug\\Output\\rotate2.bmp";
	BMP_WriteFile(newim, filepath2);
#endif // DEBUG
	//perform second X shift
	int disX2 = (int)(-tan(degrees / 2)*(outHeight / 2));
	for (int i = 0; i < outHeight / 2; i++)
	{
		double prshX = (-tan(degrees / 2));
		int shX = (int)(prshX*i);
		for (int j = 0; j < outWidth; j++)
		{
			UCHAR r, g, b;
			BMP_GetPixelRGB(newim, j, i, &r, &g, &b);
			if (flag[i*outWidth + j] == 2)
			{
				BMP_SetPixelRGB(newim, (UINT)(j - shX + disX2), (UINT)(i), r, g, b);
				if ((outWidth*i + (j - shX + disX2)) >= 0 && (outWidth*i + (j - shX + disX2)) < outWidth*outHeight) flag[(outWidth*i + (j - shX + disX2))] = 3;
			}
		}
	}
	for (int i = outHeight / 2; i < outHeight; i++)
	{
		double prshX = (-tan(degrees / 2));
		int shX = (int)(prshX*i);
		for (int j = outWidth - 1; j >= 0; j--)
		{
			UCHAR r, g, b;
			BMP_GetPixelRGB(newim, j, i, &r, &g, &b);
			if (flag[i*outWidth + j] == 2)
			{
				BMP_SetPixelRGB(newim, (UINT)(j - shX + disX2), (UINT)(i), r, g, b);
				if ((outWidth*i + (j - shX + disX2)) >= 0 && (outWidth*i + (j - shX + disX2)) < outWidth*outHeight) flag[(outWidth*i + (j - shX + disX2))] = 3;
			}
		}
	}
	for (int i = 0; i < outWidth; i++)
	{
		for (int j = 0; j < outHeight; j++)
		{
			if (flag[outWidth*j + i] != 3) BMP_SetPixelRGB(newim, (UINT)(i), (UINT)(j), 0, 0, 0);
		}
	}
#ifdef DEBUG
	BMP* debug2 = BMP_Create(outWidth, outHeight, depth);
	for (int i = 0; i < outWidth; i++)
	{
		for (int j = 0; j < outHeight; j++)
		{
			if (flag[i + outWidth*j] == 3) BMP_SetPixelRGB(debug2, (UINT)(i), (UINT)(j), 100, 100, 100);
		}
	}
	BMP_WriteFile(debug2, "C:\\Users\\HP\\Documents\\Visual Studio 2015\\Projects\\ImageProcessing\\Debug\\Output\\debug2.bmp");
	BMP_Free(debug2);
	char filepath3[100] = "C:\\Users\\HP\\Documents\\Visual Studio 2015\\Projects\\ImageProcessing\\Debug\\Output\\rotate3.bmp";
	BMP_WriteFile(newim, filepath3);
#endif // DEBUG
	free(flag);
	return newim;
}

BMP* shearRotateShell(double degrees, BMP* bmp)
{
	BMP* tmp = bmp;
	//if rotates a large degree, decomposite it into several small rotations
	if (degrees > 30 * 0.0174532925)
	{
		if (degrees <= 90 * 0.0174532925)
		{
			int times = (int)(degrees / (30 * 0.0174532925));
			double left = degrees - (30 * 0.0174532925) * times;
			for (int i = 0; i < times; i++)
			{
				tmp = shearRotate(0.0174532925 * 30, tmp);
				UINT width, height;
				USHORT depth;
				width = BMP_GetWidth(tmp);
				height = BMP_GetHeight(tmp);
				depth = BMP_GetDepth(tmp);
				int minX = width, maxX = -1, minY = height, maxY = -1;
				for (int i = 0; i < width; i++)
				{
					for (int j = 0; j < height; j++)
					{
						UCHAR r, g, b;
						BMP_GetPixelRGB(tmp, i, j, &r, &g, &b);
						if (r != 0 || g != 0 || b != 0)
						{
							if (i < minX) minX = i;
							if (i > maxX) maxX = i;
							if (j < minY) minY = j;
							if (j > maxY) maxY = j;
						}
					}
				}
				BMP* ntmp = BMP_Create(maxX - minX + 1, maxY - minY + 1, depth);
				for (int i = minX; i < maxX + 1; i++)
				{
					for (int j = minY; j < maxY + 1; j++)
					{
						UCHAR r, g, b;
						BMP_GetPixelRGB(tmp, i, j, &r, &g, &b);
						BMP_SetPixelRGB(ntmp, i - minX, j - minY, r, g, b);
					}
				}
#ifdef DEBUG
				BMP_WriteFile(ntmp, "C:\\Users\\HP\\Documents\\Visual Studio 2015\\Projects\\ImageProcessing\\Debug\\Output\\ntmp.bmp");
#endif // DEBUG
				BMP_Free(tmp);
				tmp = ntmp;
			}
			tmp = shearRotate(left, tmp);
			UINT width, height;
			USHORT depth;
			width = BMP_GetWidth(tmp);
			height = BMP_GetHeight(tmp);
			depth = BMP_GetDepth(tmp);
			int minX = width, maxX = -1, minY = height, maxY = -1;
			for (int i = 0; i < width; i++)
			{
				for (int j = 0; j < height; j++)
				{
					UCHAR r, g, b;
					BMP_GetPixelRGB(tmp, i, j, &r, &g, &b);
					if (r != 0 || g != 0 || b != 0)
					{
						if (i < minX) minX = i;
						if (i > maxX) maxX = i;
						if (j < minY) minY = j;
						if (j > maxY) maxY = j;
					}
				}
			}
			BMP* ntmp = BMP_Create(maxX - minX + 1, maxY - minY + 1, depth);
			for (int i = minX; i < maxX + 1; i++)
			{
				for (int j = minX; j < maxY + 1; j++)
				{
					UCHAR r, g, b;
					BMP_GetPixelRGB(tmp, i, j, &r, &g, &b);
					BMP_SetPixelRGB(ntmp, i - minX, j - minY, r, g, b);
				}
			}
			BMP_Free(tmp);
			tmp = ntmp;
		}
		else
		{
			if (degrees == 90 * 0.0174532925)
			{
				return transpose(tmp);
			}
			else
			{
				if (degrees < 180 * 0.0174532925)
				{
					tmp = transpose(tmp);
					tmp = shearRotateShell(degrees - 90 * 0.0174532925, tmp);
				}
				else
				{
					if (degrees == 180 * 0.0174532925)
					{
						return r180(tmp);
					}
					else
					{
						if (degrees < 270 * 0.0174532925)
						{
							tmp = r180(tmp);
							tmp = shearRotateShell(degrees - 180 * 0.0174532925, tmp);
						}
						else
						{
							if (degrees == 270 * 0.0174532925)
							{
								tmp = r180(tmp);
								tmp = transpose(tmp);
							}
							else
							{
								tmp = r180(tmp);
								tmp = transpose(tmp);
								tmp = shearRotateShell(degrees - 270 * 0.0174532925, tmp);
							}
						}
					}
				}
			}
		}

	}
	else
	{
		tmp = shearRotate(degrees, bmp);
	}
	return tmp;
}

BMP* SobelEdgeDetection(BMP* bmp, int thrhld, int coefficient)
{
	int SobelY[3][3] = { {-1,0,1},{-2,0,2},{-1,0,1} };
	int SobelX[3][3] = { {-1,-2,-1},{0,0,0},{1,2,1} };
	//Convert to Greyscale
	RGBtoBW(bmp);
	BMP_WriteFile(bmp, "C:\\Users\\HP\\Documents\\Visual Studio 2015\\Projects\\ImageProcessing\\Debug\\Tst\\BW.bmp");
	//Gaussian blur
	BMP* tmpbmp = FastGaussianBlur(bmp, 1, 3);
	BMP* swap = bmp;
	bmp = tmpbmp;
	tmpbmp = swap;
	BMP_WriteFile(bmp, "C:\\Users\\HP\\Documents\\Visual Studio 2015\\Projects\\ImageProcessing\\Debug\\Tst\\Blur.bmp");
	//BMP_Free(tmpbmp);
	//Get basic information
	UINT width, height;
	USHORT depth;
	width = BMP_GetWidth(bmp);
	height = BMP_GetHeight(bmp);
	depth = BMP_GetDepth(bmp);
	int* EdgeX = malloc(sizeof(int)*height*width);
	int* EdgeY = malloc(sizeof(int)*height*width);
	memset(EdgeX, 0, sizeof(int)*height*width);
	memset(EdgeY, 0, sizeof(int)*height*width);
	//Convolve with Sobel matrix
	for (int i = 1; i < width - 1; i++)
	{
		for (int j = 1; j < height - 1; j++)
		{
			int sumX = 0, sumY = 0;
			//Calculate gradient
			for (int iX = 0; iX < 3; iX++)
			{
				for (int iY = 0; iY < 3; iY++)
				{
					UCHAR r = 0, g = 0, b = 0;
					BMP_GetPixelRGB(bmp, i + iX - 1, j + iY - 1, &r, &g, &b);
					sumX += r*SobelX[iX][iY];
					sumY += r*SobelY[iX][iY];
				}
			}
			if (abs(sumX) > thrhld && abs(sumX) < 1020) EdgeX[i + j*width] = sumX;
			if (abs(sumY) > thrhld && abs(sumY) < 1020) EdgeY[i + j*width] = sumY;
		}
	}
	BMP* gradX = BMP_Create(width, height, depth);
	BMP* gradY = BMP_Create(width, height, depth);
	BMP* gradA = BMP_Create(width, height, depth);
	for (int i = 0; i < width; i++)
	{
		for (int j = 0; j < height; j++)
		{
			int b = (int)sqrt(EdgeX[i + j*width] * EdgeX[i + j*width] + EdgeY[i + j*width] * EdgeY[i + j*width]);
			BMP_SetPixelRGB(gradX, i, j, 128 + EdgeX[i + j*width] / coefficient, 128 + EdgeX[i + j*width] / coefficient, 128 + EdgeX[i + j*width] / coefficient);
			BMP_SetPixelRGB(gradY, i, j, 128 + EdgeY[i + j*width] / coefficient, 128 + EdgeY[i + j*width] / coefficient, 128 + EdgeY[i + j*width] / coefficient);
			BMP_SetPixelRGB(gradA, i, j, b / coefficient, b / coefficient, b / coefficient);
		}
	}
	BMP_WriteFile(gradX, "C:\\Users\\HP\\Documents\\Visual Studio 2015\\Projects\\ImageProcessing\\Debug\\Tst\\gradX.bmp");
	BMP_WriteFile(gradY, "C:\\Users\\HP\\Documents\\Visual Studio 2015\\Projects\\ImageProcessing\\Debug\\Tst\\gradY.bmp");
	BMP_WriteFile(gradA, "C:\\Users\\HP\\Documents\\Visual Studio 2015\\Projects\\ImageProcessing\\Debug\\Tst\\gradA.bmp");
	BMP_Free(gradX);
	BMP_Free(gradY);
	return gradA;
}

BMP* CannyEdgeDetection(BMP* bmp, int thrhld, int coefficient, double highShr, double lowShr, int radius)
{
	int SobelX[3][3] = { { -1,0,1 },{ -2,0,2 },{ -1,0,1 } };
	int SobelY[3][3] = { { -1,-2,-1 },{ 0,0,0 },{ 1,2,1 } };
	//Convert to Greyscale
	RGBtoBW(bmp);
	//BMP_WriteFile(bmp, "C:\\Users\\HP\\Documents\\Visual Studio 2015\\Projects\\ImageProcessing\\Debug\\Tst\\BW.bmp");
	//Gaussian blur
	BMP* tmpbmp = FastGaussianBlur(bmp, 1, 3);
	BMP* swap = bmp;
	bmp = tmpbmp;
	tmpbmp = swap;
	//BMP_WriteFile(bmp, "C:\\Users\\HP\\Documents\\Visual Studio 2015\\Projects\\ImageProcessing\\Debug\\Tst\\Blur.bmp");
	//BMP_Free(tmpbmp);
	//Get basic information
	UINT width, height;
	USHORT depth;
	width = BMP_GetWidth(bmp);
	height = BMP_GetHeight(bmp);
	depth = BMP_GetDepth(bmp);
	int* EdgeX = malloc(sizeof(int)*height*width);
	int* EdgeY = malloc(sizeof(int)*height*width);
	int* Direc = malloc(sizeof(int)*height*width);
	int* GradA = malloc(sizeof(int)*height*width);
	memset(EdgeX, 0, sizeof(int)*height*width);
	memset(EdgeY, 0, sizeof(int)*height*width);
	memset(Direc, -1, sizeof(int)*height*width);
	memset(GradA, 0, sizeof(int)*height*width);
	BMP* Canny = BMP_Create(width, height, depth);
	//Convolve with Sobel matrix
	for (int i = 1; i < width - 1; i++)
	{
		for (int j = 1; j < height - 1; j++)
		{
			int sumX = 0, sumY = 0;
			//Calculate gradient
			for (int iX = 0; iX < 3; iX++)
			{
				for (int iY = 0; iY < 3; iY++)
				{
					UCHAR r = 0, g = 0, b = 0;
					BMP_GetPixelRGB(bmp, i + iY - 1, j + iX - 1, &r, &g, &b);
					sumX += r*SobelX[iX][iY];
					sumY += r*SobelY[iX][iY];
				}
			}
			if (abs(sumX) > thrhld && abs(sumX) < 1020) EdgeX[i + j*width] = sumX;
			if (abs(sumY) > thrhld && abs(sumY) < 1020) EdgeY[i + j*width] = sumY;
			GradA[i + j*width] = sqrt(EdgeX[i + j*width] * EdgeX[i + j*width] + EdgeY[i + j*width] * EdgeY[i + j*width]);
			double theta = 0;
			if (GradA[i + j*width] != 0)
			{
				theta = atan2(EdgeY[i + j*width], EdgeX[i + j*width]);
				//Angle rounding ( long code here )
				{
					if (theta < PI / 8 && theta >= -PI / 8)
					{
						Direc[i + j*width] = 0; // 0 for East
					}
					else
					{
						if (theta >= PI / 8 && theta < PI * 3 / 8)
						{
							Direc[i + j*width] = 1;// 1 for North-East
						}
						else
						{
							if (theta >= PI * 3 / 8 && theta < PI * 5 / 8)
							{
								Direc[i + j*width] = 2;// 2 for North
							}
							else
							{
								if (theta >= PI * 5 / 8 && theta < PI * 7 / 8)
								{
									Direc[i + j*width] = 3;// 3 for North-West
								}
								else
								{
									if (theta >= PI * 7 / 8 || theta < -PI * 7 / 8)
									{
										Direc[i + j*width] = 4;// 4 for West
									}
									else
									{
										if (theta >= -PI * 7 / 8 && theta < -PI * 5 / 8)
										{
											Direc[i + j*width] = 5;// 5 for South-West
										}
										else
										{
											if (theta >= -PI * 5 / 8 && theta < -PI * 3 / 8)
											{
												Direc[i + j*width] = 6;// 6 for South
											}
											else
											{
												if (theta >= -PI * 3 / 8 && theta < -PI * 1 / 8)
												{
													Direc[i + j*width] = 7;// 7 for South-East
												}
											}
										}
									}
								}
							}
						}
					}
				}
			}
			//
			//printf("%3d", Direc[i + j*width]);
		}
		//	printf("\n");
	}
	// Second Iteration
	for (int i = 1; i < width - 1; i++)
	{
		for (int j = 0; j < height - 1; j++)
		{
			if (Direc[i + j*width] != -1)
			{
				//check neibours
				if (Direc[i + j*width] == 0 || Direc[i + j*width] == 4)
				{
					if (GradA[i + j*width] > GradA[i + 1 + j*width] && GradA[i + j*width] > GradA[i - 1 + j*width] && GradA[i + j*width] != 0)
					{
						BMP_SetPixelRGB(Canny, i, j, (UCHAR)(GradA[i + j*width] / coefficient), (UCHAR)(GradA[i + j*width] / coefficient), (UCHAR)(GradA[i + j*width] / coefficient));
					}
				}
				else
				{
					if (Direc[Direc[i + j*width] == 1 || Direc[i + j*width] == 5])
					{
						if (GradA[i + j*width] > GradA[i + 1 + (j - 1)*width] && GradA[i + j*width] > GradA[i - 1 + (j + 1)*width] && GradA[i + j*width] != 0)
						{
							BMP_SetPixelRGB(Canny, i, j, (UCHAR)(GradA[i + j*width] / coefficient), (UCHAR)(GradA[i + j*width] / coefficient), (UCHAR)(GradA[i + j*width] / coefficient));
						}
					}
					else
					{
						if (Direc[Direc[i + j*width] == 2 || Direc[i + j*width] == 6])
						{
							if (GradA[i + j*width] > GradA[i + (j + 1)*width] && GradA[i + j*width] > GradA[i + (j - 1)*width] && GradA[i + j*width] != 0)
							{
								BMP_SetPixelRGB(Canny, i, j, (UCHAR)(GradA[i + j*width] / coefficient), (UCHAR)(GradA[i + j*width] / coefficient), (UCHAR)(GradA[i + j*width] / coefficient));
							}
						}
						else
						{
							if (Direc[Direc[i + j*width] == 3 || Direc[i + j*width] == 7])
							{
								if (GradA[i + j*width] > GradA[i - 1 + (j - 1)*width] && GradA[i + j*width] > GradA[i + 1 + (j + 1)*width] && GradA[i + j*width] != 0)
								{
									BMP_SetPixelRGB(Canny, i, j, (UCHAR)(GradA[i + j*width] / coefficient), (UCHAR)(GradA[i + j*width] / coefficient), (UCHAR)(GradA[i + j*width] / coefficient));
								}
							}
						}
					}
				}
			}
		}
	}
	free(GradA);
	free(Direc);
	free(EdgeX);
	free(EdgeY);
	Canny = DoubleThreshold(Canny, highShr, lowShr, radius);
	//BMP_WriteFile(Canny, "C:\\Users\\HP\\Documents\\Visual Studio 2015\\Projects\\ImageProcessing\\Debug\\Tst\\Canny.bmp");
	return Canny;
}

BMP* DoubleThreshold(BMP* bmp, double high, double low, int radius)
{
	//CONSTANTS
	double HighShrRate = high;
	double LowShrRate = low;
	//  Generate histogram
	UINT width, height;
	USHORT depth;
	width = BMP_GetWidth(bmp);
	height = BMP_GetHeight(bmp);
	depth = BMP_GetDepth(bmp);
	int* image = malloc(sizeof(int)*height*width);
	memset(image, -1, sizeof(int)*height*width);
	int histogram[256] = { 0 };
	int maxHeight = 0;
	int totalPixel = 0;
	for (int i = 0; i < width; i++)
	{
		for (int j = 0; j < height; j++)
		{
			UCHAR r = 0, g = 0, b = 0;
			BMP_GetPixelRGB(bmp, i, j, &r, &g, &b);
			if (r != 0)
			{
				histogram[r]++;
				if (histogram[r] > maxHeight) maxHeight = histogram[r];
				totalPixel++;
			}
		}
	}
	//Write the histogram to file
	/*BMP* histoBMP = BMP_Create(256, (UINT)(maxHeight*1.2), depth);
	for (int i = 0; i < 256; i++)
	{
		for (int j = 0; j < (UINT)(maxHeight*1.2); j++)
		{
			if (j <= histogram[i]) BMP_SetPixelRGB(histoBMP, i, (UINT)(maxHeight*1.2) - 1 - j, i, i, i);
		}
	}
	BMP_WriteFile(histoBMP, "C:\\Users\\HP\\Documents\\Visual Studio 2015\\Projects\\ImageProcessing\\Debug\\Tst\\Histo.bmp");
	BMP_Free(histoBMP);
	*/
	//Calculate shreshold according to the histogram
	int countPixel = 0;
	int maxPix = (int)(totalPixel*(1 - HighShrRate));
	int minPix = (int)(totalPixel*(1 - LowShrRate));
	int maxVal = 0;
	int minVal = 0;
	int maxflag = 0, minflag = 0;
	for (int i = 255; i > 0; i--)
	{
		countPixel += histogram[i];
		if (countPixel > maxPix && maxflag == 0)
		{
			maxVal = i;
			maxflag = 1;
		}
		if (countPixel > minPix && minflag == 0)
		{
			minVal = i;
			minflag = 1;
		}
		if (minflag == 1 && maxflag == 1) break;
	}
	printf("totoal pixel= %d\nmax= %d\nmin= %d\n", totalPixel, maxVal, minVal); // DEBUG
	//Double Shreshold
	for (int i = 0; i < width; i++)
	{
		for (int j = 0; j < height; j++)
		{
			UCHAR r, g, b;
			BMP_GetPixelRGB(bmp, i, j, &r, &g, &b);
			if (r == 0) continue;
			if (r >= maxVal)
			{
				image[i + j*width] = 2; // 2 for marked edge
			}
			else
			{
				if (r <= minVal)
				{
					image[i + j*width] = 1; // 1 for marked noise
				}
				else
				{
					image[i + j*width] = 0; // 0 for marked as between
				}
			}
		}
	}
	///////////////////////////D//E//B//U//G//////////////////////////////
	/*for (int i = 0; i < height; i++)
	{
		for (int j = 0; j < width; j++)
		{
			printf("%3d", image[i*width + j]);
		}
		printf("\n");
	}*/
	//////////////////////////////////////////////////////////////////////
	//Judge Connectivity
	EgdeConnection(image, width, height, radius);
	//Write BMP
	for (int i = 0; i < width; i++)
	{
		for (int j = 0; j < height; j++)
		{
			if (image[i + j*width] != 2)
			{
				BMP_SetPixelRGB(bmp, i, j, 0, 0, 0);
			}
			else
			{
				if (image[i + j*width] == 2) BMP_SetPixelRGB(bmp, i, j, 255, 255, 255);
			}
		}
	}
	return bmp;
}

int* EgdeConnection(int* img, int width, int height, int radius)
{
	for (int i = 0; i < width; i++)
	{
		for (int j = 0; j < height; j++)
		{
			if (img[i + j*width] == 0)
			{
				IsConnected(img, width, height, i, j, radius);
			}
		}
	}
	/////////////////////////////D//E//B//U//G////////////////////////////////
	/*printf("\n\n");
	for (int i = 0; i < height; i++)
	{
		for (int j = 0; j < width; j++)
		{
			printf("%3d", img[i*width + j]);
		}
		printf("\n");
	}*/
	/////////////////////////////////////////////////////////////////////////
	return img;
}

int* IsConnected(int* img, int width, int height, int x, int y, int radius)
{
	img[x + y*width] = 3; // 3 for marked temporaryly
	for (int i = 0; i < radius; i++)
	{
		for (int j = 0; j < radius; j++)
		{
			if (!(x - (radius / 2) + i + (y - (radius / 2) + j)*width >= 0 && x - (radius / 2) + i + (y - (radius / 2) + j)*width < height*width)) continue;
			if (img[x - (radius / 2) + i + (y - (radius / 2) + j)*width] == 2)
			{
				img[x + y*width] = 2;
				return img;
			}
		}
	}
	for (int i = 0; i < radius; i++)
	{
		for (int j = 0; j < radius; j++)
		{
			if (!(x - (radius / 2) + i + (y - (radius / 2) + j)*width >= 0 && x - (radius / 2) + i + (y - (radius / 2) + j)*width < height*width)) continue;
			if (img[x - (radius / 2) + i + (y - (radius / 2) + j)*width] == 0)
			{
				img = IsConnected(img, width, height, x - (radius / 2) + i, y - (radius / 2) + j, radius);
			}
			if (img[x - (radius / 2) + i + (y - (radius / 2) + j)*width] == 2)
			{
				img[x + y*width] = 2;
				return img;
			}
		}
	}
	img[x + y*width] = 1;
	return img;
}

double maximum(double a, double b, double c, double d)
{
	double par = a;
	double var[3] = { b,c,d };
	for (int i = 0; i < 3; i++)
		par = (par > var[i]) ? par : var[i];
	return par;
}

double minimum(double a, double b, double c, double d)
{
	double par = a;
	double var[3] = { b,c,d };
	for (int i = 0; i < 3; i++)
		par = (par < var[i]) ? par : var[i];
	return par;
}

double r_compute(double var_x, double var_y, BMP* bmp, double cost, double sint)
{
	UCHAR R[4] = { 0,0,0,0 };
	UCHAR G[4] = { 0,0,0,0 };
	UCHAR B[4] = { 0,0,0,0 };

	BMP_GetPixelRGB(bmp, (unsigned long)var_x, (unsigned long)var_y, R, G, B);
	BMP_GetPixelRGB(bmp, (unsigned long)var_x + 1, (unsigned long)var_y, R + 1, G + 1, B + 1);
	BMP_GetPixelRGB(bmp, (unsigned long)var_x, (unsigned long)var_y + 1, R + 2, G + 2, B + 2);
	BMP_GetPixelRGB(bmp, (unsigned long)var_x + 1, (unsigned long)var_y + 1, R + 3, G + 3, B + 3);

	//Area average of the colors
	return ((1 - var_x + (unsigned long)var_x)*(1 - var_y + (unsigned long)var_y)*R[0] + (var_x - (unsigned long)var_x)*(1 - var_y + (unsigned long)var_y)*R[1] + (1 - var_x + (unsigned long)var_x)*(var_y - (unsigned long)var_y)*R[2] + (var_x - (unsigned long)var_x)*(var_y - (unsigned long)var_y)*R[3]);
}

double g_compute(double var_x, double var_y, BMP* bmp, double cost, double sint)
{
	UCHAR R[4] = { 0,0,0,0 };
	UCHAR G[4] = { 0,0,0,0 };
	UCHAR B[4] = { 0,0,0,0 };

	BMP_GetPixelRGB(bmp, (unsigned long)var_x, (unsigned long)var_y, R, G, B);
	BMP_GetPixelRGB(bmp, (unsigned long)var_x + 1, (unsigned long)var_y, R + 1, G + 1, B + 1);
	BMP_GetPixelRGB(bmp, (unsigned long)var_x, (unsigned long)var_y + 1, R + 2, G + 2, B + 2);
	BMP_GetPixelRGB(bmp, (unsigned long)var_x + 1, (unsigned long)var_y + 1, R + 3, G + 3, B + 3);

	//Area average of the colors
	return ((1 - var_x + (unsigned long)var_x)*(1 - var_y + (unsigned long)var_y)*G[0] + (var_x - (unsigned long)var_x)*(1 - var_y + (unsigned long)var_y)*G[1] + (1 - var_x + (unsigned long)var_x)*(var_y - (unsigned long)var_y)*G[2] + (var_x - (unsigned long)var_x)*(var_y - (unsigned long)var_y)*G[3]);
}

double b_compute(double var_x, double var_y, BMP* bmp, double cost, double sint)
{
	UCHAR R[4] = { 0,0,0,0 };
	UCHAR G[4] = { 0,0,0,0 };
	UCHAR B[4] = { 0,0,0,0 };

	BMP_GetPixelRGB(bmp, (unsigned long)var_x, (unsigned long)var_y, R, G, B);
	BMP_GetPixelRGB(bmp, (unsigned long)var_x + 1, (unsigned long)var_y, R + 1, G + 1, B + 1);
	BMP_GetPixelRGB(bmp, (unsigned long)var_x, (unsigned long)var_y + 1, R + 2, G + 2, B + 2);
	BMP_GetPixelRGB(bmp, (unsigned long)var_x + 1, (unsigned long)var_y + 1, R + 3, G + 3, B + 3);

	//Area average of the colors
	return ((1 - var_x + (unsigned long)var_x)*(1 - var_y + (unsigned long)var_y)*B[0] + (var_x - (unsigned long)var_x)*(1 - var_y + (unsigned long)var_y)*B[1] + (1 - var_x + (unsigned long)var_x)*(var_y - (unsigned long)var_y)*B[2] + (var_x - (unsigned long)var_x)*(var_y - (unsigned long)var_y)*B[3]);
}

BMP* AMRotation(BMP* bmp, double theta)
{
	//Compute some important coefficient
	double cost = cos(theta);
	double sint = sin(theta);
	double pre_width = BMP_GetWidth(bmp);
	double pre_height = BMP_GetHeight(bmp);
	USHORT depth = BMP_GetDepth(bmp);

	//Analysing size of new picture;
	double max_x = maximum((pre_width*cost - pre_height*sint), -pre_height*sint, 0.1 * 0, pre_width*cost);
	double min_x = minimum((pre_width*cost - pre_height*sint), -pre_height*sint, 0.1 * 0, pre_width*cost);
	double max_y = maximum((pre_width*sint + pre_height*cost), pre_height*cost, 0.1 * 0, pre_width*sint);
	double min_y = minimum((pre_width*sint + pre_height*cost), pre_height*cost, 0.1 * 0, pre_width*sint);
	double width = max_x - min_x;
	double height = max_y - min_y;

	BMP* newbmp = BMP_Create(width, height, depth);

	//x`=Ax + (-minx,-miny);
	for (UINT x = 0; x < width; x++)
	{
		for (UINT y = 0; y < height; y++)
		{
			double r = 0, g = 0, b = 0;

			//Compute original coordinate in newbmp system
			double pre_x = cost*(x + min_x) + sint*(y + min_y);
			double pre_y = -sint*(x + min_x) + cost*(y + min_y);

			if ((pre_x > 0 && pre_x < pre_width) && (pre_y > 0 && pre_y < pre_height))
			{
				r = r_compute(pre_x, pre_y, bmp, cost, sint);
				g = g_compute(pre_x, pre_y, bmp, cost, sint);
				b = b_compute(pre_x, pre_y, bmp, cost, sint);
			}
			BMP_SetPixelRGB(newbmp, x, y, (UCHAR)r, (UCHAR)g, (UCHAR)b);
		}
	}

	return newbmp;
}

double ABS(double a) {
	if (a >= 0) return a;
	else return -a;
}

double s(double a) {
	if (a >= 0 && a < 1) return 1 - 2 * a*a + a*a*a;
	if (a >= 1 && a < 2) return 4 - 8 * a + 5 * a*a - 3 * a*a*a;
	if (a >= 2) return 0;
}

double S(double a) {
	if (a == 0)return 1;
	else return sin(PI*a) / (PI*a);
}

double huidu(double row, double column, BMP* bmp)
{
	int width = BMP_GetWidth(bmp);
	int height = BMP_GetHeight(bmp);
	int i = row;
	int j = column;
	double u = row - (double)i;
	double v = column - (double)j;
	/*double a[4] = { s(u + 1),s(u),s(ABS(u - 1)),s(ABS(u - 2)) };
	double B[4][4] = { 0 };
	for (int x = 0; x < 4; x++)
	{
	for (int y = 0; y < 4; y++)
	{
	UCHAR r, g, b;
	if (i - 1 + x >= 0 && i - 1 + x < width&&j - 1 + y >= 0 && j - 1 + y < height)
	{
	BMP_GetPixelRGB(bmp, i - 1 + x, j - 1 + y, &r, &g, &b);
	B[x][y] = r;
	}
	}
	}
	double c[4] = { s(v + 1),s(v),s(ABS(v - 1)),s(ABS(v - 2)) };
	double d[4] = { 0 };
	double e = 0;
	for (int i = 0; i < 4; i++) {
	d[i] = a[0] * B[0][i] + a[1] * B[1][i] + a[2] * B[2][i] + a[3] * B[3][i];
	}
	for (int i = 0; i < 4; i++) {
	e += d[i] * c[i];
	}
	return e;*/
	double B[2][2] = { 0 };
	for (int x = 0; x < 2; x++)
	{
		for (int y = 0; y < 2; y++) {
			UCHAR r, g, b;
			if (i - 1 + x >= 0 && i - 1 + x < width&&j - 1 + y >= 0 && j - 1 + y < height)
			{
				BMP_GetPixelRGB(bmp, i - 1 + x, j - 1 + y, &r, &g, &b);
				B[x][y] = r;
			}
		}
	}
	return (1 - u)*(1 - v)*B[0][0] + u*(1 - v)*B[1][0] + (1 - u)*v*B[0][1] + u*v*B[1][1];

}

int huidu2(double row, double column, BMP* bmp)
{
	int width = BMP_GetWidth(bmp);
	int height = BMP_GetHeight(bmp);
	int i = row;
	int j = column;
	double u = row - (double)i;
	double v = column - (double)j;
	double a[4] = { s(ABS(u + 1)),s(ABS(u)),s(ABS(u - 1)),s(ABS(u - 2)) };
	double B[4][4] = { 0 };
	for (int x = 0; x < 4; x++)
	{
		for (int y = 0; y < 4; y++)
		{
			UCHAR r, g, b;
			if (i - 1 + x >= 0 && i - 1 + x < width&&j - 1 + y >= 0 && j - 1 + y < height)
			{
				BMP_GetPixelRGB(bmp, i - 1 + x, j - 1 + y, &r, &g, &b);
				B[x][y] = r;
			}
		}
	}
	double c[4] = { s(ABS(v + 1)),s(ABS(v)),s(ABS(v - 1)),s(ABS(v - 2)) };
	double d[4] = { 0 };
	double e = 0;
	for (int i = 0; i < 4; i++) {
		d[i] = a[0] * B[0][i] + a[1] * B[1][i] + a[2] * B[2][i] + a[3] * B[3][i];
	}
	for (int i = 0; i < 4; i++) {
		e += d[i] * c[i];
	}
	return e/320;
}

BMP* DINTRotation(BMP* bmp, double theta)
{
	RGBtoBW(bmp);
	//Compute some important coefficient
	double cost = cos(theta);
	double sint = sin(theta);
	double pre_width = BMP_GetWidth(bmp);
	double pre_height = BMP_GetHeight(bmp);
	USHORT depth = BMP_GetDepth(bmp);

	//Analysing size of new picture;
	double max_x = maximum((pre_width*cost - pre_height*sint), -pre_height*sint, 0.1 * 0, pre_width*cost);
	double min_x = minimum((pre_width*cost - pre_height*sint), -pre_height*sint, 0.1 * 0, pre_width*cost);
	double max_y = maximum((pre_width*sint + pre_height*cost), pre_height*cost, 0.1 * 0, pre_width*sint);
	double min_y = minimum((pre_width*sint + pre_height*cost), pre_height*cost, 0.1 * 0, pre_width*sint);
	double width = max_x - min_x;
	double height = max_y - min_y;

	BMP* newbmp = BMP_Create(width, height, depth);

	//x`=Ax + (-minx,-miny);
	for (UINT x = 0; x < width; x++)
	{
		for (UINT y = 0; y < height; y++)
		{
			double r = 0, g = 0, b = 0;

			//Compute original coordinate in newbmp system
			double pre_x = cost*(x + min_x) + sint*(y + min_y);
			double pre_y = -sint*(x + min_x) + cost*(y + min_y);
			if ((pre_x > 0 && pre_x < pre_width) && (pre_y > 0 && pre_y < pre_height))
			{
				r = huidu(pre_x, pre_y, bmp);
				BMP_SetPixelRGB(newbmp, x, y, (UCHAR)r, (UCHAR)r, (UCHAR)r);
			}
		}
	}
	return newbmp;
}

descripter* HarrisCorner(BMP* bmp, double thrhld,int* poic) //poic for point count
{
	const double Lambda = 0.2;
	int SobelX[3][3] = { { -1,0,1 },{ -2,0,2 },{ -1,0,1 } };
	int SobelY[3][3] = { { -1,-2,-1 },{ 0,0,0 },{ 1,2,1 } };
	//Convert to Greyscale
	RGBtoBW(bmp);
	//Get basic information
	UINT width, height;
	USHORT depth;
	width = BMP_GetWidth(bmp);
	height = BMP_GetHeight(bmp);
	depth = BMP_GetDepth(bmp);
	int* GradX = malloc(sizeof(int)*height*width); //Gradient in X direction
	int* GradX2 = malloc(sizeof(int)*height*width); //Square of gradient in X direction
	int* GradY = malloc(sizeof(int)*height*width); //Gradient in Y direction
	int* GradY2 = malloc(sizeof(int)*height*width); //Square of gradient in Y direction
	int* GradXY = malloc(sizeof(int)*height*width); //Product of gradient in X and Y direction
	int* Direc = malloc(sizeof(int)*height*width); //Direction of the gradient (Rounded to eight main directions)
	int* GradA = malloc(sizeof(int)*height*width); //Amplitude of gradient
	memset(GradX, 0, sizeof(int)*height*width); 
	memset(GradY, 0, sizeof(int)*height*width); 
	memset(GradX2, 0, sizeof(int)*height*width); 
	memset(GradY2, 0, sizeof(int)*height*width); 
	memset(GradXY, 0, sizeof(int)*height*width); 
	memset(Direc, -1, sizeof(int)*height*width); 
	memset(GradA, 0, sizeof(int)*height*width); 
	//BMP* Harris = BMP_Create(width, height, depth);
	//Convolve with Sobel matrix
	for (int i = 1; i < width - 1; i++)
	{
		for (int j = 1; j < height - 1; j++)
		{
			int sumX = 0, sumY = 0;
			//Calculate gradient
			for (int iX = 0; iX < 3; iX++)
			{
				for (int iY = 0; iY < 3; iY++)
				{
					UCHAR r = 0, g = 0, b = 0;
					BMP_GetPixelRGB(bmp, i + iY - 1, j + iX - 1, &r, &g, &b);
					sumX += r*SobelX[iX][iY];
					sumY += r*SobelY[iX][iY];
				}
			}
			if (abs(sumX) > 0 && abs(sumX) < 1020) GradX[i + j*width] = sumX;
			if (abs(sumY) > 0 && abs(sumY) < 1020) GradY[i + j*width] = sumY;
			GradA[i + j*width] = sqrt(GradX[i + j*width] * GradX[i + j*width] + GradY[i + j*width] * GradY[i + j*width]);
			double theta = 0;
			if (GradA[i + j*width] != 0)
			{
				theta = atan2(GradY[i + j*width], GradX[i + j*width]);
				//Angle rounding ( long code here. Do not expand if not necessary)
				{
					if (theta < PI / 8 && theta >= -PI / 8)
					{
						Direc[i + j*width] = 0; // 0 for East
					}
					else
					{
						if (theta >= PI / 8 && theta < PI * 3 / 8)
						{
							Direc[i + j*width] = 1;// 1 for North-East
						}
						else
						{
							if (theta >= PI * 3 / 8 && theta < PI * 5 / 8)
							{
								Direc[i + j*width] = 2;// 2 for North
							}
							else
							{
								if (theta >= PI * 5 / 8 && theta < PI * 7 / 8)
								{
									Direc[i + j*width] = 3;// 3 for North-West
								}
								else
								{
									if (theta >= PI * 7 / 8 || theta < -PI * 7 / 8)
									{
										Direc[i + j*width] = 4;// 4 for West
									}
									else
									{
										if (theta >= -PI * 7 / 8 && theta < -PI * 5 / 8)
										{
											Direc[i + j*width] = 5;// 5 for South-West
										}
										else
										{
											if (theta >= -PI * 5 / 8 && theta < -PI * 3 / 8)
											{
												Direc[i + j*width] = 6;// 6 for South
											}
											else
											{
												if (theta >= -PI * 3 / 8 && theta < -PI * 1 / 8)
												{
													Direc[i + j*width] = 7;// 7 for South-East
												}
											}
										}
									}
								}
							}
						}
					}
				}
			}
		}
	}
	//Calculating Ix2,Iy2 and Ixy
	for (int i = 0; i < width; i++)
	{
		for (int j = 0; j < height; j++)
		{
			GradX2[i + j*width] = GradX[i + j*width] * GradX[i + j*width];
			GradY2[i + j*width] = GradY[i + j*width] * GradY[i + j*width];
			GradXY[i + j*width] = GradX[i + j*width] * GradY[i + j*width];
			//printf("%d ", GradXY[i + j*width]);
		}
	}
	//Gaussian
	double *Ix2, *Iy2, *Ixy;
	Ix2 = GaussianKernel(GradX2, 1, 3, width, height);
	Iy2 = GaussianKernel(GradY2, 1, 3, width, height);
	Ixy = GaussianKernel(GradXY, 1, 3, width, height);
	GradX2 = NULL;
	GradY2 = NULL;
	GradXY = NULL;
	//Calculating response
	double *CornerStr = malloc(sizeof(double)*width*height);
	memset(CornerStr,0, sizeof(double)*width*height);
	for (int i = 0; i < width; i++)
	{
		for (int j = 0; j < height; j++)
		{
			double det_M = Ix2[i + j*width] * Iy2[i + j*width] - Ixy[i + j*width] * Ixy[i + j*width];
			double tra_M = Ix2[i + j*width] + Iy2[i + j*width];
			CornerStr[i + j*width] = det_M - Lambda*tra_M*tra_M;
			//Thresholding
			if (CornerStr[i + j*width] < thrhld)
			{
				CornerStr[i + j*width] = 0;
			}
		}
	}
	//Non-maximum Supression
	int *HarrisFin = NonMaximumSpr(CornerStr, width, height, 20);
	//Generate Descripter
	int pointcount = 0;
	descripter* EigenPoints = GenerateDescripter(HarrisFin, width, height, GradY, GradX, bmp, GradA, &pointcount);
	*poic = pointcount;
	///////////////////////////////////////////////////////////////////////////////////////////////////////

	/*for (int i = 0; i < pointcount; i++)
	{
		printf("P%d: ", i + 1);
		for (int j = 0; j < 32; j++)
		{
			printf("%-4.3f", EigenPoints[i].vector[j]);
		}
		printf("\n");
	}*/


	///////////////////////////////////////////////////////////////////////////////////////////////////////
	//Freeeeeee
	free(GradX);
	free(GradY);
	free(GradA);
	free(Direc);
	free(Ix2);
	free(Iy2);
	free(Ixy);
	//Write IMG
	/*for (int i = 0; i < width; i++)
	{
		for (int j = 0; j < height; j++)
		{
			UCHAR r, g, b;
			BMP_GetPixelRGB(bmp, i, j, &r, &g, &b);
			BMP_SetPixelRGB(Harris, i, j, r, g, b);
			if (HarrisFin[i + j*width] != 0)
			{
				for (int ki = 0; ki < 5; ki++)
				{
					for (int kj = 0; kj < 5; kj++)
					{
						BMP_SetPixelRGB(Harris, i - 2 + ki, j - 2 + kj, 0, 255, 0);
					}
				}
			}
		}
	}*/
	return EigenPoints;
}

BMP* HarrisCornerDetector(BMP* bmp, double thrhld)
{
	const double Lambda = 0.2;
	int SobelX[3][3] = { { -1,0,1 },{ -2,0,2 },{ -1,0,1 } };
	int SobelY[3][3] = { { -1,-2,-1 },{ 0,0,0 },{ 1,2,1 } };
	//Convert to Greyscale
	RGBtoBW(bmp);
	//Get basic information
	UINT width, height;
	USHORT depth;
	width = BMP_GetWidth(bmp);
	height = BMP_GetHeight(bmp);
	depth = BMP_GetDepth(bmp);
	int* GradX = malloc(sizeof(int)*height*width); //Gradient in X direction
	int* GradX2 = malloc(sizeof(int)*height*width); //Square of gradient in X direction
	int* GradY = malloc(sizeof(int)*height*width); //Gradient in Y direction
	int* GradY2 = malloc(sizeof(int)*height*width); //Square of gradient in Y direction
	int* GradXY = malloc(sizeof(int)*height*width); //Product of gradient in X and Y direction
	int* Direc = malloc(sizeof(int)*height*width); //Direction of the gradient (Rounded to eight main directions)
	int* GradA = malloc(sizeof(int)*height*width); //Amplitude of gradient
	memset(GradX, 0, sizeof(int)*height*width);
	memset(GradY, 0, sizeof(int)*height*width);
	memset(GradX2, 0, sizeof(int)*height*width);
	memset(GradY2, 0, sizeof(int)*height*width);
	memset(GradXY, 0, sizeof(int)*height*width);
	memset(Direc, -1, sizeof(int)*height*width);
	memset(GradA, 0, sizeof(int)*height*width);
	BMP* Harris = BMP_Create(width, height, depth);
	//Convolve with Sobel matrix
	for (int i = 1; i < width - 1; i++)
	{
		for (int j = 1; j < height - 1; j++)
		{
			int sumX = 0, sumY = 0;
			//Calculate gradient
			for (int iX = 0; iX < 3; iX++)
			{
				for (int iY = 0; iY < 3; iY++)
				{
					UCHAR r = 0, g = 0, b = 0;
					BMP_GetPixelRGB(bmp, i + iY - 1, j + iX - 1, &r, &g, &b);
					sumX += r*SobelX[iX][iY];
					sumY += r*SobelY[iX][iY];
				}
			}
			if (abs(sumX) > 0 && abs(sumX) < 1020) GradX[i + j*width] = sumX;
			if (abs(sumY) > 0 && abs(sumY) < 1020) GradY[i + j*width] = sumY;
			GradA[i + j*width] = sqrt(GradX[i + j*width] * GradX[i + j*width] + GradY[i + j*width] * GradY[i + j*width]);
			double theta = 0;
			if (GradA[i + j*width] != 0)
			{
				theta = atan2(GradY[i + j*width], GradX[i + j*width]);
				//Angle rounding ( long code here. Do not expand if not necessary)
				{
					if (theta < PI / 8 && theta >= -PI / 8)
					{
						Direc[i + j*width] = 0; // 0 for East
					}
					else
					{
						if (theta >= PI / 8 && theta < PI * 3 / 8)
						{
							Direc[i + j*width] = 1;// 1 for North-East
						}
						else
						{
							if (theta >= PI * 3 / 8 && theta < PI * 5 / 8)
							{
								Direc[i + j*width] = 2;// 2 for North
							}
							else
							{
								if (theta >= PI * 5 / 8 && theta < PI * 7 / 8)
								{
									Direc[i + j*width] = 3;// 3 for North-West
								}
								else
								{
									if (theta >= PI * 7 / 8 || theta < -PI * 7 / 8)
									{
										Direc[i + j*width] = 4;// 4 for West
									}
									else
									{
										if (theta >= -PI * 7 / 8 && theta < -PI * 5 / 8)
										{
											Direc[i + j*width] = 5;// 5 for South-West
										}
										else
										{
											if (theta >= -PI * 5 / 8 && theta < -PI * 3 / 8)
											{
												Direc[i + j*width] = 6;// 6 for South
											}
											else
											{
												if (theta >= -PI * 3 / 8 && theta < -PI * 1 / 8)
												{
													Direc[i + j*width] = 7;// 7 for South-East
												}
											}
										}
									}
								}
							}
						}
					}
				}
			}
		}
	}
	//Calculating Ix2,Iy2 and Ixy
	for (int i = 0; i < width; i++)
	{
		for (int j = 0; j < height; j++)
		{
			GradX2[i + j*width] = GradX[i + j*width] * GradX[i + j*width];
			GradY2[i + j*width] = GradY[i + j*width] * GradY[i + j*width];
			GradXY[i + j*width] = GradX[i + j*width] * GradY[i + j*width];
			//printf("%d ", GradXY[i + j*width]);
		}
	}
	//Gaussian
	double *Ix2, *Iy2, *Ixy;
	Ix2 = GaussianKernel(GradX2, 1, 3, width, height);
	Iy2 = GaussianKernel(GradY2, 1, 3, width, height);
	Ixy = GaussianKernel(GradXY, 1, 3, width, height);
	GradX2 = NULL;
	GradY2 = NULL;
	GradXY = NULL;
	//Calculating response
	double *CornerStr = malloc(sizeof(double)*width*height);
	memset(CornerStr, 0, sizeof(double)*width*height);
	for (int i = 0; i < width; i++)
	{
		for (int j = 0; j < height; j++)
		{
			double det_M = Ix2[i + j*width] * Iy2[i + j*width] - Ixy[i + j*width] * Ixy[i + j*width];
			double tra_M = Ix2[i + j*width] + Iy2[i + j*width];
			CornerStr[i + j*width] = det_M - Lambda*tra_M*tra_M;
			//Thresholding
			if (CornerStr[i + j*width] < thrhld)
			{
				CornerStr[i + j*width] = 0;
			}
		}
	}
	//Non-maximum Supression
	int *HarrisFin = NonMaximumSpr(CornerStr, width, height, MINSECTION);
	//Freeeeeee
	free(GradX);
	free(GradY);
	free(GradA);
	free(Direc);
	free(Ix2);
	free(Iy2);
	free(Ixy);
	//Write IMG
	for (int i = 0; i < width; i++)
	{
		for (int j = 0; j < height; j++)
		{
			UCHAR r, g, b;
			BMP_GetPixelRGB(bmp, i, j, &r, &g, &b);
			BMP_SetPixelRGB(Harris, i, j, r, g, b);
			if (HarrisFin[i + j*width] != 0)
			{
				for (int ki = 0; ki < 5; ki++)
				{
					for (int kj = 0; kj < 5; kj++)
					{
						BMP_SetPixelRGB(Harris, i - 2 + ki, j - 2 + kj, 0, 255, 0);
					}
				}
			}
		}
	}
	return Harris;
}

BMP* AMRotationParts(BMP* bmp,int x, int y, int radius, double theta,BMP** bmpt)
{
	int height = BMP_GetHeight(bmp);
	int width = BMP_GetWidth(bmp);
	int depth = BMP_GetDepth(bmp);
	BMP* newbmp = BMP_Create(radius * 2, radius * 2, depth);
	for (int i = 0; i < radius*2; i++)
	{
		for (int j = 0; j < radius*2; j++)
		{
			int r, g, b;
			BMP_GetPixelRGB(bmp, x - radius + i, y - radius + j, &r, &g, &b);
			BMP_SetPixelRGB(newbmp, i, j, r, g, b);
		}
	}
	*bmpt = newbmp;
	newbmp = AMRotation(newbmp, theta);
	//BMP_WriteFile(*bmpt, "C:\\Users\\HP\\Documents\\Visual Studio 2015\\Projects\\ImageProcessing\\Debug\\Tst\\DOWS\\Before.bmp");
	//BMP_WriteFile(newbmp, "C:\\Users\\HP\\Documents\\Visual Studio 2015\\Projects\\ImageProcessing\\Debug\\Tst\\DOWS\\After.bmp");
	return newbmp;
}

descripter* GenerateDescripter(int* Harris, int width, int height,int* GradY,int* GradX, BMP* bmp, int* GrA, int* point)
{
	const int area = 16;
	int PointCount = 0;
	//Malloc Descriptor
	descripter* KeyPoints = malloc(sizeof(descripter)*(height / MINSECTION * width / MINSECTION + 30));
	memset(KeyPoints, 0, sizeof(descripter)*(height / MINSECTION * width / MINSECTION + 30));
	for (int i = 0; i < width; i++)
	{
		for (int j = 0; j < height; j++)
		{
			int direc = 0;
			//Check eigenpoint
			if (Harris[i + j*width] != 0)
			{
				//Calculating histogram
				double histo1[36] = { 0 };
				for (int x = 0; x < 16; x++)
				{
					for (int y = 0; y < 16; y++)
					{
						double theta = atan2(GradY[i + j*width], GradX[i + j*width]);
						if (theta < 0) theta += 2 * PI;
						direc = (int)(theta / (PI / 18));
						if ((i - 8 + x + (j - 8 + y)*width) < width*height && (i - 8 + x + (j - 8 + y)*width) >= 0)
						{
							histo1[direc] += GrA[i - 8 + x + (j - 8 + y)*width];
						}
					}
				}
				//Find its peak
				int max = 0;
				int maxi = 0;
				for (int t = 0; t < 36; t++)
				{
					if (histo1[t] > max)
					{
						max = histo1[t];
						maxi = t;
					}
				}
				direc = maxi;
				KeyPoints[PointCount].x = i;
				KeyPoints[PointCount].y = j;
				//Rotating image
				KeyPoints[PointCount].direc = direc;
				double theta = direc*PI / 18 + PI / 36;
				BMP* bmpt;
				BMP* tmpbmp = AMRotationParts(bmp, i, j, area, -theta, &bmpt);
				//Convert to Greyscale
				RGBtoBW(tmpbmp);
				/*
				//Compute some important coefficient
				double cost = cos(theta);
				double sint = sin(theta);
				double pre_width = BMP_GetWidth(bmpt);
				double pre_height = BMP_GetHeight(bmpt);
				//Analysing size of new picture;
				double max_x = maximum((pre_width*cost - pre_height*sint), -pre_height*sint, 0.1 * 0, pre_width*cost);
				double min_x = minimum((pre_width*cost - pre_height*sint), -pre_height*sint, 0.1 * 0, pre_width*cost);
				double max_y = maximum((pre_width*sint + pre_height*cost), pre_height*cost, 0.1 * 0, pre_width*sint);
				double min_y = minimum((pre_width*sint + pre_height*cost), pre_height*cost, 0.1 * 0, pre_width*sint);
				*/
				//Calculating new coordinates
				//int newx = (int)(cost*area - sint*area - min_x);
				//int newy = (int)(sint*area + cost*area - min_y);
				int newwid = BMP_GetWidth(tmpbmp);
				int newhei = BMP_GetHeight(tmpbmp);
				int newx = (int)(newwid / 2);
				int newy = (int)(newwid / 2);
				//Calculating Gradients
				double *GrX = malloc(sizeof(double) * area * area);
				double *GrY = malloc(sizeof(double) * area * area);
				int *Dre = malloc(sizeof(int) * area * area);
				////////////////////////////////////////////////////////////////////////////////////////////////////////////
				/*BMP* tmbmp=BMP_Create(area, area, 24);
				for (int ix = 0; ix < area; ix++)
				{
					for (int jy = 0; jy < area; jy++)
					{
						int r, g, b;
						BMP_GetPixelRGB(tmpbmp, newx + ix - area / 2, newy + jy - area / 2, &r, &g, &b);
						BMP_SetPixelRGB(tmbmp, ix, jy, r, g, b);
					}
				}
				BMP_WriteFile(tmbmp, "C:\\Users\\HP\\Documents\\Visual Studio 2015\\Projects\\ImageProcessing\\Debug\\Tst\\DOWS\\tmbm.bmp");
				BMP_Free(tmbmp);*/
				////////////////////////////////////////////////////////////////////////////////////////////////////////////////
				for (int tx = 0; tx < area; tx++)
				{
					for (int ty = 0; ty < area; ty++)
					{
						int SobelY[3][3] = { { -1,0,1 },{ -2,0,2 },{ -1,0,1 } };
						int SobelX[3][3] = { { -1,-2,-1 },{ 0,0,0 },{ 1,2,1 } };
						int sumX = 0, sumY = 0;
						for (int iX = 0; iX < 3; iX++)
						{
							for (int iY = 0; iY < 3; iY++)
							{
								UCHAR r = 0, g = 0, b = 0;
								BMP_GetPixelRGB(tmpbmp, newx + tx + iX - area / 2, newy + ty + iY - area / 2, &r, &g, &b);
								sumX += r*SobelX[iX][iY];
								sumY += r*SobelY[iX][iY];
							}
						}
						GrX[tx + ty*area] = sumX;
						GrY[tx + ty*area] = sumY;
						double theta = 0;
						theta = atan2(GrY[tx + ty*area], GrX[tx + ty*area]);
						//Angle rounding( long code here )
						{
							if (theta < PI / 8 && theta >= -PI / 8)
							{
								Dre[tx + ty*area] = 0; // 0 for East
							}
							else
							{
								if (theta >= PI / 8 && theta < PI * 3 / 8)
								{
									Dre[tx + ty*area] = 1;// 1 for North-East
								}
								else
								{
									if (theta >= PI * 3 / 8 && theta < PI * 5 / 8)
									{
										Dre[tx + ty*area] = 2;// 2 for North
									}
									else
									{
										if (theta >= PI * 5 / 8 && theta < PI * 7 / 8)
										{
											Dre[tx + ty*area] = 3;// 3 for North-West
										}
										else
										{
											if (theta >= PI * 7 / 8 || theta < -PI * 7 / 8)
											{
												Dre[tx + ty*area] = 4;// 4 for West
											}
											else
											{
												if (theta >= -PI * 7 / 8 && theta < -PI * 5 / 8)
												{
													Dre[tx + ty*area] = 5;// 5 for South-West
												}
												else
												{
													if (theta >= -PI * 5 / 8 && theta < -PI * 3 / 8)
													{
														Dre[tx + ty*area] = 6;// 6 for South
													}
													else
													{
														if (theta >= -PI * 3 / 8 && theta < -PI * 1 / 8)
														{
															Dre[tx + ty*area] = 7;// 7 for South-East
														}
													}
												}
											}
										}
									}
								}
							}
						}
					}
				}
				double* gaussmatrix = malloc(sizeof(double) * area * area);
				double sum = GaussFunction(1.0*area / 4, area, gaussmatrix);
				for (int t = 0; t < area; t++)
				{
					for (int l = 0; l < area; l++)
					{
						GrX[t + l * area] *= gaussmatrix[t + l * area];
						//GrX[t + l * area] /= sum;
						GrY[t + l * area] *= gaussmatrix[t + l * area];
						//GrY[t + l * area] /= sum;
					}
				}
				//Generate 128-d vector
				//Calculate Histogram
				for (int lx = 0; lx < 4; lx++)
				{
					for (int ly = 0; ly < 4; ly++)
					{
						double histo[8] = { 0 };
						for (int tk = lx*(area / 4); tk < (lx + 1)*(area / 4); tk++)
						{
							for (int tl = ly*(area / 4); tl < (ly + 1)*(area / 4); tl++)
							{
								histo[Dre[tk + tl * area]] += sqrt(GrX[tk + tl * area] * GrX[tk + tl * area] + GrY[tk + tl * area] * GrY[tk + tl * area]);
							}
						}
						double isum = 0;
						for (int tmp = 0; tmp < 8; tmp++)
						{
							isum += histo[tmp] * histo[tmp];
						}
						for (int tmp = 0; tmp < 8; tmp++)
						{
							if (isum != 0) histo[tmp] /= sqrt(isum);
							KeyPoints[PointCount].vector[tmp + 8 * (lx + ly * 4)] = histo[tmp];
						}
					}
				}

				//////////////////////////////////////////////////////////
				/*printf("P%d: ", PointCount + 1);
				for (int j = 0; j < 128; j++)
				{
					printf("%-6.3f", KeyPoints[PointCount].vector[j]);
				}
				printf("\n");*/
				//////////////////////////////////////////////////////////
				PointCount++;
				//free
				free(GrX);
				free(GrY);
				free(Dre);
				BMP_Free(tmpbmp);
				free(gaussmatrix);
			}

		}
	}
	printf("\n\n");
	*point = PointCount;
	return KeyPoints;
}

pointpair* GeneratePointpair(BMP* img1, BMP* img2, double thrhld,int *paircount)
{
	int img1count, img2count;
	int flipped = 0;
	pointpair relation;
	descripter* img1des = HarrisCorner(img1, thrhld, &img1count);
	descripter* img2des = HarrisCorner(img2, thrhld, &img2count);
	if (img1count > img2count) //swap img1 and img2
	{
		int tmp = img1count;
		img1count = img2count;
		img2count = tmp;
		BMP* tmpp = img1;
		img1 = img2;
		img2 = tmpp;
		descripter* tmpdes = img1des;
		img1des = img2des;
		img2des = tmpdes;
		flipped = 1;
	}
	*paircount = 0;
	pointpair* points = malloc(sizeof(pointpair)*img1count);
	memset(points, 0, sizeof(pointpair)*img1count);
	int* flagimg2 = malloc(sizeof(int)*img2count);
	memset(flagimg2, 0, sizeof(int)*img2count);
	for (int i = 0; i < img1count; i++)
	{
		double mindist = 1e50;
		double semindist = 1e50;
		int minj = -1;
		for (int j = 0; j < img2count; j++)
		{
			if (flagimg2[j] == 0)
			{
				double dist = distof128(img1des[i], img2des[j]);
				if (dist < mindist)
				{
					semindist = mindist;
					mindist = dist;
					minj = j;
				}
			}
		}
		if (mindist/semindist <=POINTPAIR_THRHLD)
		{
			flagimg2[minj] = 1;
			if (flipped == 0)
			{
				(points[*paircount].p1).x = img1des[i].x;
				(points[*paircount].p1).y = img1des[i].y;
				(points[*paircount].p2).x = img2des[minj].x;
				(points[*paircount].p2).y = img2des[minj].y;
			}
			else
			{
				(points[*paircount].p2).x = img1des[i].x;
				(points[*paircount].p2).y = img1des[i].y;
				(points[*paircount].p1).x = img2des[minj].x;
				(points[*paircount].p1).y = img2des[minj].y;
			}
			(*paircount)++;
		}
	}
	//Free
	free(img1des);
	free(img2des);
	return points;
}

BMP* ImageStitching(BMP* img1, BMP* img2, double thrhld)
{
	int wid1 = BMP_GetWidth(img1);
	int wid2 = BMP_GetWidth(img2);
	int hei1 = BMP_GetHeight(img1);
	int hei2 = BMP_GetHeight(img2);
	int depth = BMP_GetDepth(img1);
	int paircount = 0;
	pointpair* points = RANSACMatch(img1, img2, thrhld, &paircount);
	//pointpair* points = GeneratePointpair(img1, img2, thrhld, &paircount);
	BMP* newbmp = BMP_Create(wid1 + wid2, max(hei1, hei2), depth);
	//Draw img1
	for (int i = 0; i < wid1; i++)
	{
		for (int j = 0; j < hei1; j++)
		{
			UCHAR r=0, g=0, b=0;
			BMP_GetPixelRGB(img1, i, j, &r, &g, &b);
			BMP_SetPixelRGB(newbmp, i, j, r, g, b);
		}
	}
	//Draw img2
	for (int i = wid1; i < wid1+wid2; i++)
	{
		for (int j = 0; j < hei2; j++)
		{
			UCHAR r = 0, g = 0, b = 0;
			BMP_GetPixelRGB(img2, i-wid1, j, &r, &g, &b);
			BMP_SetPixelRGB(newbmp, i, j, r, g, b);
		}
	}
	//Draw line
	for (int i = 0; i < paircount; i++)
	{
		DrawLine(newbmp, (points[i].p1).x, (points[i].p1).y, (points[i].p2).x + wid1, (points[i].p2).y, 0, 255, 0);
	}
	//Free
	free(points);
	return newbmp;
}

void DrawLine(BMP* bmp, int inix, int iniy, int desx, int desy, int r, int g, int b)
{
	double disy = 1.0*(desy - iniy) / abs(desx - inix);
	int dire = -1;
	if (desx > inix) dire = 1;
	double y = iniy;
	for (int i = inix;i <=desx; i+=dire)
	{
		BMP_SetPixelRGB(bmp, i, (int)y, r, g, b);
		y += disy;
	}
}

double distof128(descripter desc1, descripter desc2) //distance of 36-d vector
{
	double dist = 0;
	for (int i = 0; i < 128; i++)
	{
		dist += (desc1.vector[i] - desc2.vector[i])*(desc1.vector[i] - desc2.vector[i]);
	}
	return sqrt(dist);
}

int* NonMaximumSpr(double* img, int width, int height, int windowSize)
{
	int *ans = malloc(sizeof(int)*width*height);
	memset(ans, 0, sizeof(int)*width*height);
	for (int i = 0; i < width/windowSize; i++)
	{
		for (int j = 0; j < height / windowSize; j++)
		{
			double max = 0;
			int maxpos = -1;
			for (int ti = 0; ti < windowSize; ti++)
			{
				for (int tj = 0; tj < windowSize; tj++)
				{
					int pos = i*windowSize + ti + (j*windowSize + tj)*width;
					if (pos >= 0 && pos < width*height)
					{
						if (img[i*windowSize + ti + (j*windowSize + tj)*width] > max)
						{
							max = img[i*windowSize + ti + (j*windowSize + tj)*width];
							maxpos = pos;
						}
					}
				}
			}
			if (max > 10) ans[maxpos] = 1;
		}
	}
	return ans;
}


double determinant(double a11, double a12, double a13, double a21, double a22, double a23, double a31, double a32, double a33)
{
	double result = a11*a22*a33 + a12*a23*a31 + a21*a32*a13 - a13*a22*a31 - a12*a21*a33 - a11*a23*a32;
	return result;
}

int matchjudgement(affinematrix matrix, pointpair*naivepair, int paircount, int* match)
{
	int result = 0;
	for (int i = 0; i < paircount; i++)
	{
		//Affine Transformation;
		point p;
		p.x = (matrix.A1)*(naivepair[i].p1.x) + (matrix.B1)*(naivepair[i].p1.y) + matrix.e;
		p.y = matrix.A2*naivepair[i].p1.x + matrix.B2*naivepair[i].p1.y + matrix.f;
		if (((naivepair[i].p2).x - p.x)*((naivepair[i].p2).x - p.x) + ((naivepair[i].p2).y - p.y)*((naivepair[i].p2).y - p.y) < RANSAC_RADIUS*RANSAC_RADIUS && ((naivepair[i].p2).x - p.x)*((naivepair[i].p2).x - p.x) + ((naivepair[i].p2).y - p.y)*((naivepair[i].p2).y - p.y) >= 0)
		{
			result++;
			match[i] = 1;
		}
		else
		{
			match[i] = 0;
		}
	}

	return result;
}

affinematrix RANSACAffm(BMP* img1, BMP* img2, double thrhld, int times, double radius)
{
	affinematrix best_affm;
	best_affm.A1 = 0;
	best_affm.A2 = 0;
	best_affm.B1 = 0;
	best_affm.B2 = 0;
	best_affm.e = 0;
	best_affm.f = 0;
	int num_of_match = 0;
	int count = 0;
	int paircount = 0;
	affinematrix result = best_affm;
	pointpair* naivepair = GeneratePointpair(img1, img2, thrhld, paircount);
	//Initial 3 point pairs;
	for (int i = 0; i < times; i++)
	{
		int initialpair[3] = { -1,-1,-1 };

		initialpair[0] = rand() % paircount;
		initialpair[1] = rand() % paircount;
		while (initialpair[1] == initialpair[0] && initialpair[1] == -1)
			initialpair[1] = rand() % paircount;
		initialpair[2] = rand() % paircount;
		while (initialpair[2] == initialpair[0] || initialpair[2] == initialpair[1] && initialpair[2] == -1)
			initialpair[2] = rand() % paircount;

		pointpair pair1 = naivepair[initialpair[0]];
		pointpair pair2 = naivepair[initialpair[1]];
		pointpair pair3 = naivepair[initialpair[2]];

		//Compute the model affine matrix;
		/*double D1 = determinant(pair1.p1.x, pair1.p1.y, 1, pair2.p1.x, pair2.p1.y, 1, pair3.p1.x, pair3.p1.y, 1);
		if (D1 != 0)
		{
			best_affm.A1 = determinant(pair1.p2.x, pair1.p1.y, 1, pair2.p2.x, pair2.p1.y, 1, pair3.p2.x, pair3.p1.y, 1) / D1;
			best_affm.B1 = determinant(pair1.p1.x, pair1.p2.x, 1, pair2.p1.x, pair2.p2.x, 1, pair3.p1.x, pair3.p2.x, 1) / D1;
			best_affm.e = determinant(pair1.p1.x, pair1.p2.y, pair1.p2.x, pair2.p1.x, pair2.p2.y, pair2.p2.x, pair3.p1.x, pair3.p2.y, pair3.p2.x) / D1;
			best_affm.A2 = determinant(pair1.p2.y, pair1.p1.y, 1, pair2.p2.y, pair2.p1.y, 1, pair3.p2.y, pair3.p1.y, 1) / D1;
			best_affm.B2 = determinant(pair1.p1.x, pair1.p2.y, 1, pair2.p1.x, pair2.p2.y, 1, pair3.p1.x, pair3.p2.y, 1) / D1;
			best_affm.f = determinant(pair1.p1.x, pair1.p1.y, pair1.p2.y, pair2.p1.x, pair2.p1.y, pair2.p2.y, pair3.p1.x, pair3.p1.y, pair3.p2.y) / D1;
		}*/

		double D1 = determinant(pair1.p1.x, pair1.p1.y, 1, pair2.p1.x, pair2.p1.y, 1, pair3.p1.x, pair3.p1.y, 1);
		if (D1 != 0)
		{
		best_affm.A1 = determinant(pair1.p2.x, pair1.p1.y, 1, pair2.p2.x, pair2.p1.y, 1, pair3.p2.x, pair3.p1.y, 1) / D1;
		best_affm.B1 = determinant(pair1.p1.x, pair1.p2.x, 1, pair2.p1.x, pair2.p2.x, 1, pair3.p1.x, pair3.p2.x, 1) / D1;
		best_affm.e = determinant(pair1.p1.x, pair1.p1.y, pair1.p2.x, pair2.p1.x, pair2.p1.y, pair2.p2.x, pair3.p1.x, pair3.p1.y, pair3.p2.x) / D1;
		best_affm.A2 = determinant(pair1.p2.y, pair1.p1.y, 1, pair2.p2.y, pair2.p1.y, 1, pair3.p2.y, pair3.p1.y, 1) / D1;
		best_affm.B2 = determinant(pair1.p1.x, pair1.p2.y, 1, pair2.p1.x, pair2.p2.y, 1, pair3.p1.x, pair3.p2.y, 1) / D1;
		best_affm.f = determinant(pair1.p1.x, pair1.p1.y, pair1.p2.y, pair2.p1.x, pair2.p1.y, pair2.p2.y, pair3.p1.x, pair3.p1.y, pair3.p2.y) / D1;
		}

		//Compute the number of points successfully matched;
		if ((best_affm.A1)*(best_affm.A1) + (best_affm.A2) *(best_affm.A2) + (best_affm.B1) *(best_affm.B1) + (best_affm.B2) *(best_affm.B2) + (best_affm.e) *(best_affm.e) + (best_affm.f) *(best_affm.f) == 0)
			continue;
		else
		{
			int* match = malloc(sizeof(int)*(paircount));
			num_of_match = matchjudgement(best_affm, naivepair, paircount, match);
			free(match);
		}

		//Update the affm;
		if (num_of_match > count)
		{
			count = num_of_match;
			result = best_affm;
		}
	}

	free(naivepair);
	return result;
}

pointpair* RANSACMatch(BMP* img1, BMP* img2, double thrhld, int* RANSACcount)
{
	affinematrix best_affm;
	best_affm.A1 = 0;
	best_affm.A2 = 0;
	best_affm.B1 = 0;
	best_affm.B2 = 0;
	best_affm.e = 0;
	best_affm.f = 0;
	int num_of_match = 0;
	int count = 0;
	int paircount = 0;
	affinematrix result;
	pointpair* naivepair = GeneratePointpair(img1, img2, thrhld, &paircount);
	int* match = malloc(sizeof(int)*(paircount));

	//Initial 3 point pairs;
	for (int i = 0; i < RANSAC_TIMES; i++)
	{
		int initialpair[3] = { -1,-1,-1 };

		initialpair[0] = rand() % paircount;
		initialpair[1] = rand() % paircount;
		while (initialpair[1] == initialpair[0] && initialpair[1] == -1)
			initialpair[1] = rand() % paircount;
		initialpair[2] = rand() % paircount;
		while (initialpair[2] == initialpair[0] || initialpair[2] == initialpair[1] && initialpair[2] == -1)
			initialpair[2] = rand() % paircount;

		pointpair pair1 = naivepair[initialpair[0]];
		pointpair pair2 = naivepair[initialpair[1]];
		pointpair pair3 = naivepair[initialpair[2]];

		//Compute the model affine matrix;
		/*double D1 = determinant(pair1.p1.x, pair1.p1.y, 1, pair2.p1.x, pair2.p1.y, 1, pair3.p1.x, pair3.p1.y, 1);
		if (D1 != 0)
		{
			best_affm.A1 = determinant(pair1.p2.x, pair1.p1.y, 1, pair2.p2.x, pair2.p1.y, 1, pair3.p2.x, pair3.p1.y, 1) / D1;
			best_affm.B1 = determinant(pair1.p1.x, pair1.p2.x, 1, pair2.p1.x, pair2.p2.x, 1, pair3.p1.x, pair3.p2.x, 1) / D1;
			best_affm.e = determinant(pair1.p1.x, pair1.p2.y, pair1.p2.x, pair2.p1.x, pair2.p2.y, pair2.p2.x, pair3.p1.x, pair3.p2.y, pair3.p2.x) / D1;
			best_affm.A2 = determinant(pair1.p2.y, pair1.p1.y, 1, pair2.p2.y, pair2.p1.y, 1, pair3.p2.y, pair3.p1.y, 1) / D1;
			best_affm.B2 = determinant(pair1.p1.x, pair1.p2.y, 1, pair2.p1.x, pair2.p2.y, 1, pair3.p1.x, pair3.p2.y, 1) / D1;
			best_affm.f = determinant(pair1.p1.x, pair1.p1.y, pair1.p2.y, pair2.p1.x, pair2.p1.y, pair2.p2.y, pair3.p1.x, pair3.p1.y, pair3.p2.y) / D1;
		}*/

		double D1 = determinant(pair1.p1.x, pair1.p1.y, 1, pair2.p1.x, pair2.p1.y, 1, pair3.p1.x, pair3.p1.y, 1);
		if (D1 != 0)
		{
			best_affm.A1 = determinant(pair1.p2.x, pair1.p1.y, 1, pair2.p2.x, pair2.p1.y, 1, pair3.p2.x, pair3.p1.y, 1) / D1;
			best_affm.B1 = determinant(pair1.p1.x, pair1.p2.x, 1, pair2.p1.x, pair2.p2.x, 1, pair3.p1.x, pair3.p2.x, 1) / D1;
			best_affm.e = determinant(pair1.p1.x, pair1.p1.y, pair1.p2.x, pair2.p1.x, pair2.p1.y, pair2.p2.x, pair3.p1.x, pair3.p1.y, pair3.p2.x) / D1;
			best_affm.A2 = determinant(pair1.p2.y, pair1.p1.y, 1, pair2.p2.y, pair2.p1.y, 1, pair3.p2.y, pair3.p1.y, 1) / D1;
			best_affm.B2 = determinant(pair1.p1.x, pair1.p2.y, 1, pair2.p1.x, pair2.p2.y, 1, pair3.p1.x, pair3.p2.y, 1) / D1;
			best_affm.f = determinant(pair1.p1.x, pair1.p1.y, pair1.p2.y, pair2.p1.x, pair2.p1.y, pair2.p2.y, pair3.p1.x, pair3.p1.y, pair3.p2.y) / D1;
		}

		//Compute the number of points successfully matched;
		if ((best_affm.A1) == 0 && (best_affm.A2) == 0 && (best_affm.B1) == 0 && (best_affm.B2) == 0 && (best_affm.e) == 0 && (best_affm.f) == 0)
			continue;
		else
			num_of_match = matchjudgement(best_affm, naivepair, paircount, match);
		//Update the affm;
		if (num_of_match > count)
		{
			count = num_of_match;
			result = best_affm;
		}

	}

	count = matchjudgement(result, naivepair, paircount, match);
	pointpair* Matchedpairs = malloc(sizeof(pointpair)*count);
	*RANSACcount = 0;
	for (int k = 0; k < paircount; k++)
	{
		if (match[k] == 1)
		{
			Matchedpairs[*RANSACcount] = naivepair[k];
			(*RANSACcount)++;
		}
	}

	free(match);
	free(naivepair);
	return Matchedpairs;
}
