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
#define RANSAC_RADIUS 10
#define RANSAC_TIMES 1000000
#define USM_GAUSSIAN_RADIUS 3//Gaussian radius in USM
#define LAMDA 3 //lamda is a coefficient controling intensity of USM

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
BMP* ImageComparing(BMP* img1, BMP* img2, double thrhld);
pointpair* GeneratePointpair(BMP* img1, BMP* img2, double thrhld, int *paircount);
double determinant(double a11, double a12, double a13, double a21, double a22, double a23, double a31, double a32, double a33);
affinematrix RANSACAffm(BMP* img1, BMP* img2, double thrhld);
pointpair* RANSACMatch(BMP* img1, BMP* img2, double thrhld, int* RANSACcount, affinematrix* matrix);
BMP* ImageStitching(BMP* img1, BMP* img2, double thrhld);
pointpair* RANSACHomo_match(BMP* img1, BMP* img2, double thrhld, int* RANSACcount, homographicmatrix* transform);
double Exponentiation(double base, int exponential);
int matchjudgement(affinematrix matrix, pointpair *naivepair, int paircount, int* match);//Judge how many point pairs has been matched under certain affine matrix;
BMP* NaiveUSM(BMP* bmp);
double GaussianElimination_determinant(double**matrix, int size);
double Laplacian_expansion(double** matrix, int initial_size, int col, int** seize);
int matchjudgement_homo(homographicmatrix matrix, pointpair* naivepair, int paircount, int* match);
void* SeparateFrequency(int* img, int radius, int hei, int wid, int** highfre, int** lowfre);
BMP* ImageBlending(int* img1, BMP* img2, int img1wid, int img1hei, int img1x0, int img1y0, int img2x0, int img2y0, BMP* newimg);
point homographicTransform(point a, homographicmatrix matrix);
void AreaInterpolation(homographicmatrix matrix, BMP* img1, int x, int y, double* R, double* G, double* B);
BMP* ImageBlending_rgb(BMP* Originalbmp, BMP* img1, BMP* img2, int img1x0, int img1y0, int img2x0, int img2y0, BMP* newimg, homographicmatrix mat);

int main()
{
	BMP *bmp, *bmp1;
	char filename[] = "C:\\Users\\HP\\Documents\\Visual Studio 2015\\Projects\\ImageProcessing\\Debug\\swmp1.bmp";
	char filename1[] = "C:\\Users\\HP\\Documents\\Visual Studio 2015\\Projects\\ImageProcessing\\Debug\\swmp2.bmp";
	bmp = BMP_ReadFile(filename);
	BMP_CHECK_ERROR(stderr, -1);
	bmp1 = BMP_ReadFile(filename1);
	BMP_CHECK_ERROR(stderr, -1);
	/////////////////////////////////////////////////////////////////////////
	//Your code in between
	//srand(time(NULL));
	srand(1);
	//BMP* newbmp = ImageStitching(bmp, bmp1, 100000);
	//char filepath[100] = "C:\\Users\\HP\\Documents\\Visual Studio 2015\\Projects\\ImageProcessing\\Debug\\Tst\\DOWS\\NewStitch.bmp";
	//BMP_WriteFile(newbmp, filepath);
	//BMP_Free(newbmp);
	BMP* newbmp = HarrisCornerDetector(bmp, 100000);
	BMP_WriteFile(newbmp, "C:\\Users\\HP\\Documents\\Visual Studio 2015\\Projects\\ImageProcessing\\Debug\\Tst\\DOWS\\HarrisSWMP1.bmp");
	BMP_Free(newbmp);
	newbmp = HarrisCornerDetector(bmp1, 100000);
	BMP_WriteFile(newbmp, "C:\\Users\\HP\\Documents\\Visual Studio 2015\\Projects\\ImageProcessing\\Debug\\Tst\\DOWS\\HarrisSWMP2.bmp");
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
	for (int i = 0; i < width / k; i++)
	{
		for (int j = 0; j < height / k; j++)
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

double* GaussianKernelDouble(double* img, double theta, int radius, int width, int height)
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
	free(img1);
	return ans;
}

int* intGaussianKernel(int* img, double theta, int radius, int width, int height)
{
	double* gaussmatrix = malloc(sizeof(double)*radius*radius);
	double sum = GaussFunction(theta, radius, gaussmatrix);
	sum = 0;
	for (int i = 0; i < radius; i++)
	{
		sum += gaussmatrix[i + radius*(radius / 2)];
	}
	int *ans = malloc(sizeof(int)*width*height);
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
	int SobelY[3][3] = { { -1,0,1 },{ -2,0,2 },{ -1,0,1 } };
	int SobelX[3][3] = { { -1,-2,-1 },{ 0,0,0 },{ 1,2,1 } };
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
	return e / 320;
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

descripter* HarrisCorner(BMP* bmp, double thrhld, int* poic) //poic for point count
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

BMP* AMRotationParts(BMP* bmp, int x, int y, int radius, double theta, BMP** bmpt)
{
	int height = BMP_GetHeight(bmp);
	int width = BMP_GetWidth(bmp);
	int depth = BMP_GetDepth(bmp);
	BMP* newbmp = BMP_Create(radius * 2, radius * 2, depth);
	for (int i = 0; i < radius * 2; i++)
	{
		for (int j = 0; j < radius * 2; j++)
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

descripter* GenerateDescripter(int* Harris, int width, int height, int* GradY, int* GradX, BMP* bmp, int* GrA, int* point)
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

pointpair* GeneratePointpair(BMP* img1, BMP* img2, double thrhld, int *paircount)
{
	int img1count, img2count;
	int flipped = 0;
	pointpair relation;
	printf("Processing Picture 1 ...\n");
	descripter* img1des = HarrisCorner(img1, thrhld, &img1count);
	printf("Processing Picture 2 ...\n");
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
		if (mindist / semindist <= POINTPAIR_THRHLD)
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

BMP* ImageComparing(BMP* img1, BMP* img2, double thrhld)
{
	int wid1 = BMP_GetWidth(img1);
	int wid2 = BMP_GetWidth(img2);
	int hei1 = BMP_GetHeight(img1);
	int hei2 = BMP_GetHeight(img2);
	int depth = BMP_GetDepth(img1);
	int paircount = 0;
	affinematrix transform;
	pointpair* points = RANSACMatch(img1, img2, thrhld, &paircount, &transform);
	//pointpair* points = GeneratePointpair(img1, img2, thrhld, &paircount);
	BMP* newbmp = BMP_Create(wid1 + wid2, max(hei1, hei2), depth);
	//Draw img1
	for (int i = 0; i < wid1; i++)
	{
		for (int j = 0; j < hei1; j++)
		{
			UCHAR r = 0, g = 0, b = 0;
			BMP_GetPixelRGB(img1, i, j, &r, &g, &b);
			BMP_SetPixelRGB(newbmp, i, j, r, g, b);
		}
	}
	//Draw img2
	for (int i = wid1; i < wid1 + wid2; i++)
	{
		for (int j = 0; j < hei2; j++)
		{
			UCHAR r = 0, g = 0, b = 0;
			BMP_GetPixelRGB(img2, i - wid1, j, &r, &g, &b);
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
	for (int i = inix; i <= desx; i += dire)
	{
		BMP_SetPixelRGB(bmp, i, (int)y, r, g, b);
		y += disy;
	}
}

double distof128(descripter desc1, descripter desc2) //distance of 128-d vector
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
	for (int i = 0; i < width / windowSize; i++)
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

point affineTransform(point a, affinematrix matrix)
{
	point p;
	//Affine Transformation;
	p.x = (matrix.A1)*(a.x) + (matrix.B1)*(a.y) + matrix.e;
	p.y = matrix.A2*(a.x) + matrix.B2*(a.y) + matrix.f;
	return p;
}

point homographicTransform(point a, homographicmatrix matrix)
{
	point p;
	//homographic Transformation;
	double var = matrix.H31*a.x + matrix.H32*a.y + 1;
	p.x = ((matrix.H11)*(a.x) + (matrix.H12)*(a.y) + matrix.H13) / var;
	p.y = (matrix.H21*a.x + matrix.H22*a.y + matrix.H23) / var;
	return p;
}

affinematrix RANSACAffm(BMP* img1, BMP* img2, double thrhld)
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
	printf("Generating pointpair ...\n");
	pointpair* naivepair = GeneratePointpair(img1, img2, thrhld, &paircount);
	printf("Calculating RANSAC match ...\n");
	//Initial 3 point pairs;
	for (int i = 0; i < RANSAC_TIMES; i++)
	{
		int initialpair[3] = { -1,-1,-1 };

		initialpair[0] = rand() % paircount;
		initialpair[1] = rand() % paircount;
		while (initialpair[1] == initialpair[0] || initialpair[1] == -1)
			initialpair[1] = rand() % paircount;
		initialpair[2] = rand() % paircount;
		while (initialpair[2] == initialpair[0] || initialpair[2] == initialpair[1] || initialpair[2] == -1)
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

pointpair* RANSACMatch(BMP* img1, BMP* img2, double thrhld, int* RANSACcount, affinematrix* matrix)
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
	printf("Generating pointpair ...\n");
	pointpair* naivepair = GeneratePointpair(img1, img2, thrhld, &paircount);
	int* match = malloc(sizeof(int)*(paircount));
	printf("Calculating RANSAC match ...\n");
	//Initial 3 point pairs;
	for (int i = 0; i < RANSAC_TIMES; i++)
	{
		int initialpair[3] = { -1,-1,-1 };

		initialpair[0] = rand() % paircount;
		initialpair[1] = rand() % paircount;
		while (initialpair[1] == initialpair[0] || initialpair[1] == -1)
			initialpair[1] = rand() % paircount;
		initialpair[2] = rand() % paircount;
		while (initialpair[2] == initialpair[0] || initialpair[2] == initialpair[1] || initialpair[2] == -1)
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

	*matrix = result;
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

int maximumint(int a, int b, int c, int d, int e)
{
	int par = a;
	int var[4] = { b,c,d,e };
	for (int i = 0; i < 4; i++)
		par = (par > var[i]) ? par : var[i];
	return par;
}

int minimumint(int a, int b, int c, int d, int e)
{
	int par = a;
	int var[4] = { b,c,d,e };
	for (int i = 0; i < 4; i++)
		par = (par < var[i]) ? par : var[i];
	return par;
}

BMP* ImageStitching(BMP* bmp1, BMP* bmp2, double thrhld)
{
	int wid1 = BMP_GetWidth(bmp1);
	int wid2 = BMP_GetWidth(bmp2);
	int hei1 = BMP_GetHeight(bmp1);
	int hei2 = BMP_GetHeight(bmp2);
	int depth = BMP_GetDepth(bmp1);
	int paircount = 0;
	BMP* backbmp1 = BMP_Create(wid1, hei1, depth);
	//Copy BMP1
	for (int i = 0; i < wid1; i++)
	{
		for (int j = 0; j < hei1; j++)
		{
			int r = 0, g = 0, b = 0;
			BMP_GetPixelRGB(bmp1, i, j, &r, &g, &b);
			BMP_SetPixelRGB(backbmp1, i, j, r, g, b);
		}
	}
	BMP* backbmp2 = BMP_Create(wid2, hei2, depth);
	//Copy BMP2
	for (int i = 0; i < wid2; i++)
	{
		for (int j = 0; j < hei2; j++)
		{
			int r = 0, g = 0, b = 0;
			BMP_GetPixelRGB(bmp2, i, j, &r, &g, &b);
			BMP_SetPixelRGB(backbmp2, i, j, r, g, b);
		}
	}
	homographicmatrix transform;
	int count;
	pointpair* points = RANSACHomo_match(bmp1, bmp2, thrhld, &count, &transform);
	//Transform Boundary
	point zero_zero_tr = { 0,0 }, width_zero_tr = { wid1,0 }, zero_height_tr = { 0,hei1 }, width_height_tr = { wid1,hei1 };
	zero_zero_tr = homographicTransform(zero_zero_tr, transform);
	width_zero_tr = homographicTransform(width_zero_tr, transform);
	zero_height_tr = homographicTransform(zero_height_tr, transform);
	width_height_tr = homographicTransform(width_height_tr, transform);
	//Calculate Boundary
	int minx = minimumint(zero_zero_tr.x, width_zero_tr.x, width_height_tr.x, zero_height_tr.x, 0);
	int maxx = maximumint(zero_zero_tr.x, width_zero_tr.x, width_height_tr.x, zero_height_tr.x, wid2);
	int miny = minimumint(zero_zero_tr.y, width_zero_tr.y, width_height_tr.y, zero_height_tr.y, 0);
	int maxy = maximumint(zero_zero_tr.y, width_zero_tr.y, width_height_tr.y, zero_height_tr.y, hei2);
	int newwid = maxx - minx;
	int newhei = maxy - miny;
	//Create New Image
	BMP* newimg = BMP_Create(newwid, newhei, depth);
	//Transform img1
	int img1x0 = (int)minimum(zero_zero_tr.x, width_zero_tr.x, width_height_tr.x, zero_height_tr.x) - minx;
	int img1wid = (int)maximum(zero_zero_tr.x, width_zero_tr.x, width_height_tr.x, zero_height_tr.x) - minx - img1x0;
	int img1y0 = (int)minimum(zero_zero_tr.y, width_zero_tr.y, width_height_tr.y, zero_height_tr.y) - miny;
	int img1hei = (int)maximum(zero_zero_tr.y, width_zero_tr.y, width_height_tr.y, zero_height_tr.y) - miny - img1y0;
	int* img1_r = malloc((img1hei + 1)*(img1wid + 1) * sizeof(int));
	int* img1_g = malloc((img1hei + 1)*(img1wid + 1) * sizeof(int));
	int* img1_b = malloc((img1hei + 1)*(img1wid + 1) * sizeof(int));
	///////////////////////////////
	BMP* imgbmp1 = BMP_Create(img1wid, img1hei, depth);
	///////////////////////////////
	memset(img1_r, -1, img1hei*img1wid * sizeof(int));
	memset(img1_g, -1, img1hei*img1wid * sizeof(int));
	memset(img1_b, -1, img1hei*img1wid * sizeof(int));
	for (int i = 0; i < img1wid; i++)
	{
		for (int j = 0; j < img1hei; j++)
		{
			double r = 0, g = 0, b = 0;
			if (judge_edge(bmp1, transform, i, j, img1wid, img1hei) == 1)
			{
				int i1 = i + (int)minimum(zero_zero_tr.x, width_zero_tr.x, width_height_tr.x, zero_height_tr.x);
				int j1 = j + (int)minimum(zero_zero_tr.y, width_zero_tr.y, width_height_tr.y, zero_height_tr.y);
				AreaInterpolation(transform, backbmp1, i1, j1, &r, &g, &b);
				BMP_SetPixelRGB(imgbmp1, i, j, r, g, b);
				if (i + j*img1wid >= 0 && i + j*img1wid < img1wid*img1hei)
				{
					img1_r[i + j*img1wid] = r;
					img1_g[i + j*img1wid] = g;
					img1_b[i + j*img1wid] = b;
				}
			}
			else
				BMP_SetPixelRGB(imgbmp1, i, j, (UCHAR)r, (UCHAR)g, (UCHAR)b);
		}
	}
	//////////////////////////////////////
	printf("Writing img1.bmp ...\n");
	BMP_WriteFile(imgbmp1, "C:\\Users\\HP\\Documents\\Visual Studio 2015\\Projects\\ImageProcessing\\Debug\\Tst\\DOWS\\IMG1.bmp");
	//////////////////////////////////////
	//Transform img2
	int img2x0 = -minx, img2y0 = -miny;
	//BLEND
	ImageBlending_rgb(backbmp1, imgbmp1, backbmp2, img1x0, img1y0, img2x0, img2y0, newimg, transform, img1wid, img1hei);
	BMP_Free(imgbmp1);
	///////////////////////////////////////////////////////////////////////////////////
	BMP* newbmp = BMP_Create(wid1 + wid2, max(hei1, hei2), depth);
	//Draw img1
	for (int i = 0; i < wid1; i++)
	{
		for (int j = 0; j < hei1; j++)
		{
			int r = 0, g = 0, b = 0;
			BMP_GetPixelRGB(bmp1, i, j, &r, &g, &b);
			BMP_SetPixelRGB(newbmp, i, j, r, g, b);
		}
	}
	//Draw img2
	for (int i = wid1; i < wid1 + wid2; i++)
	{
		for (int j = 0; j < hei2; j++)
		{
			int r = 0, g = 0, b = 0;
			BMP_GetPixelRGB(bmp2, i - wid1, j, &r, &g, &b);
			BMP_SetPixelRGB(newbmp, i, j, r, g, b);
		}
	}
	//Draw line
	for (int i = 0; i < count; i++)
	{
		DrawLine(newbmp, (points[i].p1).x, (points[i].p1).y, (points[i].p2).x + wid1, (points[i].p2).y, 0, 255, 0);
	}
	printf("Writing Compare.bmp ...\n");
	BMP_WriteFile(newbmp, "C:\\Users\\HP\\Documents\\Visual Studio 2015\\Projects\\ImageProcessing\\Debug\\Tst\\DOWS\\COMP.bmp");
	//Free
	free(points);
	BMP_Free(newbmp);
	free(img1_r);
	free(img1_g);
	free(img1_b);
	/////////////////////////////////////////////////////////////////////////////////////////////
	return newimg;
}

BMP* NaiveUSM(BMP* bmp)
{
	double width = BMP_GetWidth(bmp);
	double height = BMP_GetHeight(bmp);
	double depth = BMP_GetDepth(bmp);
	BMP* newbmp = BMP_Create(width, height, depth);

	/* Memorize Blurred image's information*/
	BMP* Gaussianbmp = BMP_Create(width, height, depth);
	Gaussianbmp = FastGaussianBlur(bmp, 2 * PI, USM_GAUSSIAN_RADIUS);
	int** pre_R = malloc(sizeof(int*)*width);
	int** pre_G = malloc(sizeof(int*)*width);
	int** pre_B = malloc(sizeof(int*)*width);

	for (int i = 0; i < width; i++)
	{
		pre_R[i] = malloc(sizeof(int)*height);
		pre_G[i] = malloc(sizeof(int)*height);
		pre_B[i] = malloc(sizeof(int)*height);
	}

	for (int i = 0; i < width; i++)
	{
		for (int j = 0; j < height; j++)
			BMP_GetPixelRGB(bmp, i, j, &pre_R[i][j], &pre_G[i][j], &pre_B[i][j]);
	}

	for (int i = 0; i < width; i++)
	{
		for (int j = 0; j < height; j++)
		{
			int grad_r = 0;
			int grad_g = 0;
			int grad_b = 0;

			int R[9] = { 0,0,0,0,0,0,0,0,0 };
			int G[9] = { 0,0,0,0,0,0,0,0,0 };
			int B[9] = { 0,0,0,0,0,0,0,0,0 };

			BMP_GetPixelRGB(bmp, i, j, &grad_r, &grad_g, &grad_b);

			if (i == 0 || j == 0 || i == width - 1 || j == height - 1)
				BMP_SetPixelRGB(newbmp, i, j, grad_r, grad_g, grad_b);
			else
			{
				BMP_GetPixelRGB(Gaussianbmp, i - 1, j - 1, R, G, B);
				BMP_GetPixelRGB(Gaussianbmp, i, j - 1, R + 1, G + 1, B + 1);
				BMP_GetPixelRGB(Gaussianbmp, i + 1, j - 1, R + 2, G + 2, B + 2);
				BMP_GetPixelRGB(Gaussianbmp, i - 1, j, R + 3, G + 3, B + 3);
				BMP_GetPixelRGB(Gaussianbmp, i, j, R + 4, G + 4, B + 4);
				BMP_GetPixelRGB(Gaussianbmp, i + 1, j, R + 5, G + 5, B + 5);
				BMP_GetPixelRGB(Gaussianbmp, i - 1, j + 1, R + 6, G + 6, B + 6);
				BMP_GetPixelRGB(Gaussianbmp, i, j + 1, R + 7, G + 7, B + 7);
				BMP_GetPixelRGB(Gaussianbmp, i + 1, j + 1, R + 8, G + 8, B + 8);

				int model[9] = { -1,-1,-1,-1,8,-1,-1,-1,-1 };
				for (int var = 0; var < 9; var++)
				{
					grad_r += (int)(LAMDA*model[var] * R[var] + 0.5);
					grad_g += (int)(LAMDA* model[var] * G[var] + 0.5);
					grad_b += (int)(LAMDA*model[var] * B[var] + 0.5);
				}

				/*RGB's gradient computation*/
				if (grad_r < 0)
					grad_r = -grad_r;
				if (grad_r > 255)
					grad_r = 255;
				if (grad_g < 0)
					grad_g = -grad_g;
				if (grad_g > 255)
					grad_g = 255;
				if (grad_b < 0)
					grad_b = -grad_b;
				if (grad_b > 255)
					grad_b = 255;

				BMP_SetPixelRGB(newbmp, i, j, grad_r, grad_g, grad_b);
			}
		}
	}

	/*Free allocate memories*/
	for (int i = 0; i < height; i++)
	{
		free(pre_R[i]);
		free(pre_G[i]);
		free(pre_B[i]);
	}

	free(pre_R);
	free(pre_G);
	free(pre_B);

	BMP_Free(Gaussianbmp);
	return newbmp;
}

double Exponentiation(double base, int exponential)//Calculation a double's integer exponentiation;
{
	double result = 1;
	for (int i = 0; i < exponential; i++)
		result *= base;
	return result;
}

double GaussianElimination_determinant(double**matrix, int size)
{
	int sign = 1;
	for (int i = 0; i < size; i++)//(i+1)-th pivot;
	{
		///////////////////////////////////////////////////////
		/*for (int i = 0; i < 8; i++)
		{
		for (int j = 0; j < 8; j++)
		{
		printf("%16.4f", matrix[i][j]);
		}
		printf("\n\n");
		}
		printf("\n\n");*/
		//////////////////////////////////////////////////////
		//Judge whether (i+1)-th pivot equals to 0;
		if (matrix[i][i] == 0)//Interchange 2 rows;
		{
			if (i == size - 1)
				return 0;
			else
			{
				int change_row = i + 1;
				while (change_row < size)
				{
					if (matrix[change_row][1] == 0) change_row++;
					else
					{
						break;
					}
				}

				if (change_row == size)//Loop to the last row...
					return 0;
				else
				{
					for (int j = 0; j < size; j++)
					{
						double mid_row = matrix[i][j];
						matrix[i][j] = matrix[change_row][j];
						matrix[change_row][j] = mid_row;
					}
					sign = -sign;
				}
			}
			//Finish row's interchange;
		}
		for (int j = i + 1; j < size; j++)
		{

			double coefficient = -matrix[j][i] / matrix[i][i];
			for (int k = 0; k < size; k++)
				matrix[j][k] = coefficient*matrix[i][k] + matrix[j][k];
		}
	}

	double result = 1;
	for (int i = 0; i < size; i++)
		result *= matrix[i][i];
	result *= sign;
	return result;
}

double Laplacian_expansion(double** matrix, int initial_size, int col, int** seize)
//row and col begin at 0, size begin with 1;
//Array "seize": Choosen row/col: 1; else: 0.
{
	if (col == initial_size - 1)
	{
		for (int i = 0; i < initial_size; i++)
		{
			for (int j = 0; j < initial_size; j++)
			{
				if (seize[i][j] == 0)
					return matrix[i][j];
			}
		}
	}
	else
	{
		double result = 0;
		for (int i = 0; i < initial_size; i++)
		{
			if (seize[i][col] == 0)
			{
				for (int j = 0; j < initial_size; j++)
				{
					seize[i][j] = 1;
					seize[j][col] = 1;
				}

				result += Exponentiation(-1, i + col)*matrix[i][col] * Laplacian_expansion(matrix, initial_size, col + 1, seize);

				for (int j = 0; j < initial_size; j++)
				{
					seize[i][j] = 0;
					seize[j][col] = 0;
				}
			}

		}
		return result;
	}
}


homographicmatrix RANSACHomo(BMP* img1, BMP* img2, double thrhld)
{
	homographicmatrix best;
	best.H11 = 0;
	best.H12 = 0;
	best.H13 = 0;
	best.H21 = 0;
	best.H22 = 0;
	best.H23 = 0;
	best.H31 = 0;
	best.H32 = 0;

	int num_of_match = 0;
	int count = 0;
	int paircount = 0;
	homographicmatrix result = best;
	pointpair* naivepair = GeneratePointpair(img1, img2, thrhld, &paircount);

	for (int i = 0; i < RANSAC_TIMES; i++)
	{

		//Initial 4 point pairs;
		int initialpair[4] = { -1,-1,-1,-1 };

		initialpair[0] = rand() % paircount;
		initialpair[1] = rand() % paircount;
		while (initialpair[1] == initialpair[0] || initialpair[1] == -1)
			initialpair[1] = rand() % paircount;
		initialpair[2] = rand() % paircount;
		while (initialpair[2] == initialpair[0] || initialpair[2] == initialpair[1] || initialpair[2] == -1)
			initialpair[2] = rand() % paircount;
		initialpair[3] = rand() % paircount;
		while (initialpair[3] == initialpair[0] || initialpair[3] == initialpair[1] || initialpair[3] == initialpair[2] || initialpair[2] == -1)
			initialpair[3] = rand() % paircount;

		pointpair pairs[4];
		for (int k = 0; k < 4; k++)
			pairs[k] = naivepair[initialpair[k]];

		//From computation...
		//Au = v;
		//u = [ H11, H12, H13, H21, H22, H23, H31, H32]^T
		//v = [ x12, y12, x22, y22, x32, y32, x42, y42]^T
		//A: [ x11, y11, 1, 0, 0, 0, -x11x12, -x12y11 ]
		//   [ 0, 0, 0, x11, y11, 1, -x11y12, -y11y12 ]
		//	 [ x21, y21, 1, 0, 0, 0, -x21x22, -x22y21 ]
		//   [ 0, 0, 0, x21, y21, 1, -x21y22, -y21y22 ]
		//	 [ x31, y31, 1, 0, 0, 0, -x31x32, -x32y31 ]
		//   [ 0, 0, 0, x31, y31, 1, -x31y32, -y31y32 ]
		//	 [ x41, y41, 1, 0, 0, 0, -x41x42, -x42y41 ]
		//   [ 0, 0, 0, x41, y41, 1, -x41y42, -y41y42 ]

		double** com_mat = malloc(sizeof(double*) * 8);
		for (int k = 0; k < 8; k++)
			com_mat[k] = malloc(sizeof(double) * 8);

		//Generate matrix A;
		for (int k = 0; k < 4; k++)
		{
			com_mat[2 * k][0] = pairs[k].p1.x;
			com_mat[2 * k][1] = pairs[k].p1.y;
			com_mat[2 * k][2] = 1;
			com_mat[2 * k][3] = 0;
			com_mat[2 * k][4] = 0;
			com_mat[2 * k][5] = 0;
			com_mat[2 * k][6] = -pairs[k].p1.x*pairs[k].p2.x;
			com_mat[2 * k][7] = -pairs[k].p2.x*pairs[k].p1.y;
			com_mat[2 * k + 1][0] = 0;
			com_mat[2 * k + 1][1] = 0;
			com_mat[2 * k + 1][2] = 0;
			com_mat[2 * k + 1][3] = pairs[k].p1.x;
			com_mat[2 * k + 1][4] = pairs[k].p1.y;
			com_mat[2 * k + 1][5] = 1;
			com_mat[2 * k + 1][6] = -pairs[k].p1.x*pairs[k].p2.y;
			com_mat[2 * k + 1][7] = -pairs[k].p1.y*pairs[k].p2.y;
		}

		//Solution vector;
		double vec[8] = { pairs[0].p2.x , pairs[0].p2.y , pairs[1].p2.x , pairs[1].p2.y , pairs[2].p2.x , pairs[2].p2.y , pairs[3].p2.x , pairs[3].p2.y };
		//Result matrix...
		double entry_H[8];
		for (int k = 0; k < 8; k++)
			entry_H[k] = 0;

		//First step of Cramer: computing the determinant of matrix;
		double D1 = GaussianElimination_determinant(com_mat, 8);

		for (int k = 0; k < 8; k++)
		{
			double mid[8];
			for (int var = 0; var < 8; var++)//...I can't find a good name for my variable...
			{
				mid[var] = com_mat[var][k];
				com_mat[var][k] = vec[k];
			}

			entry_H[k] = GaussianElimination_determinant(com_mat, 8) / D1;


			for (int var = 0; var < 8; var++)//Back!
				com_mat[var][k] = mid[var];
		}

		//free generated matrix;
		for (int k = 0; k < 8; k++)
			free(com_mat[k]);
		free(com_mat);

		//Assignment results to matrix;
		best.H11 = entry_H[0];
		best.H12 = entry_H[1];
		best.H13 = entry_H[2];
		best.H21 = entry_H[3];
		best.H22 = entry_H[4];
		best.H23 = entry_H[5];
		best.H31 = entry_H[6];
		best.H32 = entry_H[7];

		//Compute the number of points successfully matched;
		if (best.H11 == 0 && best.H12 == 0 && best.H13 == 0 && best.H21 == 0 && best.H22 == 0 && best.H23 == 0 && best.H31 == 0 && best.H32 == 0)
			continue;
		else
		{
			int* match = malloc(sizeof(int)*(paircount));
			num_of_match = matchjudgement_homo(best, naivepair, paircount, match);
			free(match);
		}

		//Update the affm;
		if (num_of_match >= count)
		{
			count = num_of_match;
			result = best;
		}
	}

	free(naivepair);
	return result;

}


pointpair* RANSACHomo_match(BMP* img1, BMP* img2, double thrhld, int* RANSACcount, homographicmatrix* transform)
{
	homographicmatrix best;
	best.H11 = 0;
	best.H12 = 0;
	best.H13 = 0;
	best.H21 = 0;
	best.H22 = 0;
	best.H23 = 0;
	best.H31 = 0;
	best.H32 = 0;

	int num_of_match = 0;
	int count = 0;
	int paircount = 0;
	homographicmatrix result = best;
	printf("Generating point pair ... \n\n");
	pointpair* naivepair = GeneratePointpair(img1, img2, thrhld, &paircount);
	printf("Finding homography using RANSAC ...\n");
	for (int i = 0; i < RANSAC_TIMES; i++)
	{

		//Initial 4 point pairs;
		int* initialpair = malloc(sizeof(int) * 4);
		for (int k = 0; k < 4; k++)
			initialpair[k] = -1;

		initialpair[0] = rand() % paircount;
		initialpair[1] = rand() % paircount;
		while (initialpair[1] == initialpair[0] || initialpair[1] == -1)
			initialpair[1] = rand() % paircount;
		initialpair[2] = rand() % paircount;
		while (initialpair[2] == initialpair[0] || initialpair[2] == initialpair[1] || initialpair[2] == -1)
			initialpair[2] = rand() % paircount;
		initialpair[3] = rand() % paircount;
		while (initialpair[3] == initialpair[0] || initialpair[3] == initialpair[1] || initialpair[3] == initialpair[2] || initialpair[2] == -1)
			initialpair[3] = rand() % paircount;

		pointpair pairs[4];
		for (int k = 0; k < 4; k++)
			pairs[k] = naivepair[initialpair[k]];



		//From computation...
		//Au = v;
		//u = [ H11, H12, H13, H21, H22, H23, H31, H32]^T
		//v = [ x12, y12, x22, y22, x32, y32, x42, y42]^T
		//A: [ pairs[k].p1.x,  pairs[k].p1.y, 1, 0, 0, 0, - pairs[k].p1.x*pairs[k].p2.x,- pairs[k].p1.y*pairs[k].p2.x ]
		//   [ 0, 0, 0,  pairs[k].p1.x,  pairs[k].p1.y, 1,- pairs[k].p1.x*pairs[k].p2.y, - pairs[k].p1.y*pairs[k].p2.y ]
		//	 [ x21, y21, 1, 0, 0, 0, -x21x22, -x22y21 ]
		//   [ 0, 0, 0, x21, y21, 1, -x21y22, -y21y22 ]
		//	 [ x31, y31, 1, 0, 0, 0, -x31x32, -x32y31 ]
		//   [ 0, 0, 0, x31, y31, 1, -x31y32, -y31y32 ]
		//	 [ x41, y41, 1, 0, 0, 0, -x41x42, -x42y41 ]
		//   [ 0, 0, 0, x41, y41, 1, -x41y42, -y41y42 ]

		double** com_mat = malloc(sizeof(double*) * 8);
		for (int k = 0; k < 8; k++)
			com_mat[k] = malloc(sizeof(double) * 8);

		//Generate matrix A;
		for (int k = 0; k < 4; k++)
		{
			com_mat[2 * k][0] = pairs[k].p1.x;
			com_mat[2 * k][1] = pairs[k].p1.y;
			com_mat[2 * k][2] = 1;
			com_mat[2 * k][3] = 0;
			com_mat[2 * k][4] = 0;
			com_mat[2 * k][5] = 0;
			com_mat[2 * k][6] = -pairs[k].p1.x*pairs[k].p2.x;
			com_mat[2 * k][7] = -pairs[k].p1.y*pairs[k].p2.x;
			com_mat[2 * k + 1][0] = 0;
			com_mat[2 * k + 1][1] = 0;
			com_mat[2 * k + 1][2] = 0;
			com_mat[2 * k + 1][3] = pairs[k].p1.x;
			com_mat[2 * k + 1][4] = pairs[k].p1.y;
			com_mat[2 * k + 1][5] = 1;
			com_mat[2 * k + 1][6] = -pairs[k].p1.x*pairs[k].p2.y;
			com_mat[2 * k + 1][7] = -pairs[k].p1.y*pairs[k].p2.y;
		}

		//Solution vector;
		double vec[8] = { pairs[0].p2.x , pairs[0].p2.y , pairs[1].p2.x , pairs[1].p2.y , pairs[2].p2.x , pairs[2].p2.y , pairs[3].p2.x , pairs[3].p2.y };
		//Result matrix...
		double entry_H[8] = { 0 };

		//First step of Cramer: computing the determinant of matrix;
		double D1 = GaussianElimination_determinant(com_mat, 8);
		if (D1 == 0)
		{
			//free generated matrix;
			for (int k = 0; k < 8; k++)
				free(com_mat[k]);
			free(com_mat);
			free(initialpair);
			continue;
		}

		//Regenerate com_mat
		for (int k = 0; k < 4; k++)
		{
			com_mat[2 * k][0] = pairs[k].p1.x;
			com_mat[2 * k][1] = pairs[k].p1.y;
			com_mat[2 * k][2] = 1;
			com_mat[2 * k][3] = 0;
			com_mat[2 * k][4] = 0;
			com_mat[2 * k][5] = 0;
			com_mat[2 * k][6] = -pairs[k].p1.x*pairs[k].p2.x;
			com_mat[2 * k][7] = -pairs[k].p1.y*pairs[k].p2.x;
			com_mat[2 * k + 1][0] = 0;
			com_mat[2 * k + 1][1] = 0;
			com_mat[2 * k + 1][2] = 0;
			com_mat[2 * k + 1][3] = pairs[k].p1.x;
			com_mat[2 * k + 1][4] = pairs[k].p1.y;
			com_mat[2 * k + 1][5] = 1;
			com_mat[2 * k + 1][6] = -pairs[k].p1.x*pairs[k].p2.y;
			com_mat[2 * k + 1][7] = -pairs[k].p1.y*pairs[k].p2.y;
		}

		for (int k = 0; k < 8; k++)
		{
			//Generate tmp_mat
			double** tmp_mat = malloc(sizeof(double*) * 8);
			for (int k = 0; k < 8; k++)
				tmp_mat[k] = malloc(sizeof(double) * 8);

			for (int i = 0; i < 8; i++)
			{
				for (int j = 0; j < 8; j++)
				{
					tmp_mat[i][j] = com_mat[i][j];
				}
			}

			for (int var = 0; var < 8; var++)//...I can't find a good name for my variable...
			{
				tmp_mat[var][k] = vec[var];
			}

			entry_H[k] = GaussianElimination_determinant(tmp_mat, 8) / D1;

			//free generated matrix;
			for (int k = 0; k < 8; k++)
				free(tmp_mat[k]);
			free(tmp_mat);

			/*for (int var = 0; var < 8; var++)//Back!
			com_mat[var][k] = mid[var];*/

		}

		//free generated matrix;
		for (int k = 0; k < 8; k++)
			free(com_mat[k]);
		free(com_mat);


		//Assignment results to matrix;
		best.H11 = entry_H[0];
		best.H12 = entry_H[1];
		best.H13 = entry_H[2];
		best.H21 = entry_H[3];
		best.H22 = entry_H[4];
		best.H23 = entry_H[5];
		best.H31 = entry_H[6];
		best.H32 = entry_H[7];

		//Compute the number of points successfully matched;
		if (best.H11 == 0 && best.H12 == 0 && best.H13 == 0 && best.H21 == 0 && best.H22 == 0 && best.H23 == 0 && best.H31 == 0 && best.H32 == 0)
			continue;
		else
		{
			int* match = malloc(sizeof(int)*(paircount));
			for (int var = 0; var < paircount; var++)
				match[var] = 0;
			int paircount1 = paircount;
			num_of_match = matchjudgement_homo(best, naivepair, paircount1, match);
			free(match);
		}

		//Update the homom;
		if (num_of_match >= count)
		{
			count = num_of_match;
			result = best;
		}
		free(initialpair);
	}

	int* match = malloc(sizeof(int)*(paircount));
	count = matchjudgement_homo(result, naivepair, paircount, match);
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
	*transform = result;
	free(match);
	free(naivepair);
	return Matchedpairs;

}

int matchjudgement_homo(homographicmatrix matrix, pointpair* naivepair, int paircount, int* match)
{
	int result = 0;
	for (int i = 0; i < paircount; i++)
	{
		//homographic Transformation;
		point p;
		double var = matrix.H31*naivepair[i].p1.x + matrix.H32*naivepair[i].p1.y + 1;
		p.x = ((matrix.H11)*(naivepair[i].p1.x) + (matrix.H12)*(naivepair[i].p1.y) + matrix.H13) / var;
		p.y = (matrix.H21*naivepair[i].p1.x + matrix.H22*naivepair[i].p1.y + matrix.H23) / var;
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

BMP* ImageBlending(int* img1, BMP* img2, int img1wid, int img1hei, int img1x0, int img1y0, int img2x0, int img2y0, BMP* newimg)
{
	int img2wid = BMP_GetWidth(img2);
	int img2hei = BMP_GetHeight(img2);
	int depth = BMP_GetDepth(img2);
	int regx0 = max(img1x0, img2x0);
	int regwid = min(img1x0 + img1wid, img2x0 + img2wid) - regx0;
	int regy0 = max(img1y0, img2y0);
	int reghei = min(img1y0 + img1hei, img2y0 + img2hei) - regy0;
	//Create the overlap part of img1
	int* img1part = malloc(regwid*reghei * sizeof(int));
	memset(img1part, -1, sizeof(int)*regwid*reghei); // -1 means transparent
	for (int i = 0; i < regwid; i++)
	{
		for (int j = 0; j < reghei; j++)
		{
			int intensity = 0;
			intensity = img1[(regx0 - img1x0 + i) + (regy0 - img1y0 + j)*img1wid];
			img1part[i + j*regwid] = intensity;
		}
	}
	//Create the overlap part of img2
	int* img2part = malloc(regwid*reghei * sizeof(int));
	for (int i = 0; i < regwid; i++)
	{
		for (int j = 0; j < reghei; j++)
		{
			int r = 0, g = 0, b = 0;
			BMP_GetPixelRGB(img2, regx0 - img2x0 + i, regy0 - img2y0 + j, &r, &g, &b);
			img2part[i + j*regwid] = r;
		}
	}
	/*//Separating high frequency and low frequency
	int* img1high = NULL;
	int* img1low = NULL;
	SeparateFrequency(img1part, 3, reghei, regwid, &img1high, &img1low);
	int* img2high = NULL;
	int* img2low = NULL;
	SeparateFrequency(img2part, 3, reghei, regwid, &img2high, &img2low);
	//Blending high frequency
	int* blendhigh = malloc(sizeof(int)*reghei*regwid);
	memset(blendhigh, -1, sizeof(int)*reghei*regwid);
	for (int i = 0; i < regwid; i++)
	{
	for (int j = 0; j < reghei; j++)
	{
	int intensity1 = img1high[i + j*regwid];
	int intensity2 = img2high[i + j*regwid];
	blendhigh[i + j*regwid] = (intensity1 + intensity2) / 2;//binary blending
	}
	}
	//Blending Low frequency
	int* blendlow = malloc(sizeof(int)*reghei*regwid);
	for (int i = 0; i < regwid; i++)
	{
	double img1weight = 1 - (1.0*i / regwid);
	double img2weight = 1.0*i / regwid;
	for (int j = 0; j < reghei; j++)
	{
	int flag = img1part[i + j*regwid];
	int i1 = img1low[i + j*regwid];
	int i2 = img2low[i + j*regwid];
	if (flag >= 0)	blendlow[i + j*regwid] = i1*img1weight + i2*img2weight; //Linear blending
	else
	{
	blendlow[i + j*regwid] = i2;
	}
	}
	}*/
	/*
	//Reconstruct
	int* blendreg = malloc(sizeof(int)*reghei*regwid);
	for (int i = 0; i < regwid; i++)
	{
	for (int j = 0; j < reghei; j++)
	{
	int i1 = blendhigh[i + j*regwid];
	int i2 = blendlow[i + j*regwid];
	blendreg[i + j*regwid] = blendhigh[i + j*regwid] + blendlow[i + j*regwid];
	}
	}
	*/
	//Reconstruct
	int* blendreg = malloc(sizeof(int)*reghei*regwid);
	for (int i = 0; i < regwid; i++)
	{
		double img1weight = 1 - (1.0*i / regwid);
		double img2weight = 1.0*i / regwid;
		for (int j = 0; j < reghei; j++)
		{
			int i1 = img1part[i + j*regwid];
			int i2 = img2part[i + j*regwid];
			if (i1 >= 0)	blendreg[i + j*regwid] = i1*img1weight + i2*img2weight; //Linear blending
			else
			{
				blendreg[i + j*regwid] = i2;
			}
		}
	}
	//Write Newimg
	int newwid = BMP_GetWidth(newimg);
	int newhei = BMP_GetHeight(newimg);
	for (int i = 0; i < newwid; i++)
	{
		for (int j = 0; j < newhei; j++)
		{
			if (i >= regx0&&i < (regx0 + regwid) && j >= regy0&&j < (regy0 + reghei)) // if in overlap region
			{
				int i1 = blendreg[(i - regx0) + (j - regy0)*regwid];
				BMP_SetPixelRGB(newimg, i, j, i1, i1, i1);
			}
			else
			{
				if (i >= img1x0&&i < (img1x0 + img1wid) && j >= img1y0&&j < (img1y0 + img1hei)) // if in img1
				{
					int i1 = img1[(i - img1x0) + (j - img1y0)*img1wid];
					if (i1 >= 0)	BMP_SetPixelRGB(newimg, i, j, i1, i1, i1);
				}
				else
				{
					if (i >= img2x0&&i < (img2x0 + img2wid) && j >= img2y0&&j < (img2y0 + img2hei)) // if in img2
					{
						int r1 = 0, g1 = 0, b1 = 0;
						BMP_GetPixelRGB(img2, i - img2x0, j - img2y0, &r1, &g1, &b1);
						BMP_SetPixelRGB(newimg, i, j, r1, g1, b1);
					}
				}
			}
		}
	}
	//Freeeeeeee!!!
	free(img1part);
	free(img2part);
	/*free(img1high);
	free(img1low);
	free(img2high);
	free(img2low);
	free(blendhigh);
	free(blendlow);*/
	free(blendreg);
}

void* SeparateFrequency(int* img, int radius, int hei, int wid, int** highfre, int** lowfre) //DO NOT MALLOC HIGHFRE AND LOWFRE
{
	(*highfre) = malloc(sizeof(int)*hei*wid);
	memset(*highfre, 0, sizeof(int)*hei*wid);
	(*lowfre) = intGaussianKernel(img, 1.0*radius / 3, radius, wid, hei);
	for (int i = 0; i < wid; i++)
	{
		for (int j = 0; j < hei; j++)
		{
			int intensity = 0;
			int lowintensity = 0;
			intensity = img[i + j*wid];
			lowintensity = (*lowfre)[i + j*wid];
			(*highfre)[i + j*wid] = intensity - lowintensity;
		}
	}

}

void AreaInterpolation(homographicmatrix matrix, BMP* img1, int x, int y, double* R, double* G, double* B)
{
	//Determinant 3*3;
	double determ_of_matrix = matrix.H11*matrix.H22 + matrix.H12*matrix.H23*matrix.H31 + matrix.H21*matrix.H32*matrix.H13 - matrix.H13*matrix.H22*matrix.H31 - matrix.H32*matrix.H23*matrix.H11 - matrix.H12*matrix.H21;
	double inverse_matrix[3][3];
	//Adjugate matrix...
	inverse_matrix[0][0] = (matrix.H22 - matrix.H23*matrix.H32) / determ_of_matrix;
	inverse_matrix[1][0] = -(matrix.H21 - matrix.H23*matrix.H31) / determ_of_matrix;
	inverse_matrix[2][0] = (matrix.H21*matrix.H32 - matrix.H22*matrix.H31) / determ_of_matrix;
	inverse_matrix[0][1] = -(matrix.H12 - matrix.H32*matrix.H13) / determ_of_matrix;
	inverse_matrix[1][1] = (matrix.H11 - matrix.H13*matrix.H31) / determ_of_matrix;
	inverse_matrix[2][1] = -(matrix.H11*matrix.H32 - matrix.H12*matrix.H31) / determ_of_matrix;
	inverse_matrix[0][2] = (matrix.H12*matrix.H23 - matrix.H22*matrix.H13) / determ_of_matrix;
	inverse_matrix[1][2] = -(matrix.H11*matrix.H23 - matrix.H13*matrix.H21) / determ_of_matrix;
	inverse_matrix[2][2] = (matrix.H11*matrix.H22 - matrix.H12*matrix.H21) / determ_of_matrix;

	//First mapping...
	/*double homo_z = (matrix.H31*x + matrix.H32*y + 1);
	*homo_x = (int)((matrix.H11*x + matrix.H12*y + matrix.H13) / homo_z);
	*homo_y = (int)((matrix.H21*x + matrix.H22*y + matrix.H23) / homo_z);
	int var1 = (int)((matrix.H11*x + matrix.H12*y + matrix.H13));
	int var2 = (int)((matrix.H21*x + matrix.H22*y + matrix.H23));*/

	double coefficient = (inverse_matrix[2][0] * x + inverse_matrix[2][1] * y + inverse_matrix[2][2]);
	//Second re-mapping;
	double interpolate_x = (inverse_matrix[0][0] * x + inverse_matrix[0][1] * y + inverse_matrix[0][2]) / coefficient;
	double interpolate_y = (inverse_matrix[1][0] * x + inverse_matrix[1][1] * y + inverse_matrix[1][2]) / coefficient;

	//Get color value of 4 coners;
	UCHAR pre_r[4];
	UCHAR pre_g[4];
	UCHAR pre_b[4];
	BMP_GetPixelRGB(img1, (unsigned long)interpolate_x, (unsigned long)interpolate_y, pre_r, pre_g, pre_b);//LU
	BMP_GetPixelRGB(img1, (unsigned long)interpolate_x + 1, (unsigned long)interpolate_y, pre_r + 1, pre_g + 1, pre_b + 1);//RU
	BMP_GetPixelRGB(img1, (unsigned long)interpolate_x, (unsigned long)interpolate_y + 1, pre_r + 2, pre_g + 2, pre_b + 2);//LB
	BMP_GetPixelRGB(img1, (unsigned long)interpolate_x + 1, (unsigned long)interpolate_y + 1, pre_r + 3, pre_g + 3, pre_b + 3);//RB

	*R = (1 - interpolate_x + (unsigned long)interpolate_x)*(1 - interpolate_y + (unsigned long)interpolate_y)*pre_r[0] + (interpolate_x - (unsigned long)interpolate_x)*(1 - interpolate_y + (unsigned long)interpolate_y)*pre_r[1] + (1 - interpolate_x + (unsigned long)interpolate_x)*(interpolate_y - (unsigned long)interpolate_y)*pre_r[2] + (interpolate_x - (unsigned long)interpolate_x)*(interpolate_y - (unsigned long)interpolate_y)*pre_r[3];
	*G = (1 - interpolate_x + (unsigned long)interpolate_x)*(1 - interpolate_y + (unsigned long)interpolate_y)*pre_g[0] + (interpolate_x - (unsigned long)interpolate_x)*(1 - interpolate_y + (unsigned long)interpolate_y)*pre_g[1] + (1 - interpolate_x + (unsigned long)interpolate_x)*(interpolate_y - (unsigned long)interpolate_y)*pre_g[2] + (interpolate_x - (unsigned long)interpolate_x)*(interpolate_y - (unsigned long)interpolate_y)*pre_g[3];
	*B = (1 - interpolate_x + (unsigned long)interpolate_x)*(1 - interpolate_y + (unsigned long)interpolate_y)*pre_b[0] + (interpolate_x - (unsigned long)interpolate_x)*(1 - interpolate_y + (unsigned long)interpolate_y)*pre_b[1] + (1 - interpolate_x + (unsigned long)interpolate_x)*(interpolate_y - (unsigned long)interpolate_y)*pre_b[2] + (interpolate_x - (unsigned long)interpolate_x)*(interpolate_y - (unsigned long)interpolate_y)*pre_b[3];

}

int judge_edge(BMP* img, homographicmatrix matrix, int x, int y, int img1wid, int img1hei)
{
	double pre_width = BMP_GetWidth(img);
	double pre_height = BMP_GetHeight(img);
	//Coorsinate transformation;
	point LL, RL, LU, RU;
	LL.x = 0;
	LL.y = 0;
	RL.x = pre_width;
	RL.y = 0;
	LU.x = 0;
	LU.y = pre_height;
	RU.x = pre_width;
	RU.y = pre_height;
	LL = homographicTransform(LL, matrix);
	RL = homographicTransform(RL, matrix);
	LU = homographicTransform(LU, matrix);
	RU = homographicTransform(RU, matrix);
	double min_x = minimum(LL.x, RL.x, LU.x, RU.x);
	double min_y = minimum(LL.y, RL.y, LU.y, RU.y);

	//Transformation of testing point;
	point test_p;
	test_p.x = x + min_x;
	test_p.y = y + min_y;

	//Determinant 3*3;
	double determ_of_matrix = matrix.H11*matrix.H22 + matrix.H12*matrix.H23*matrix.H31 + matrix.H21*matrix.H32*matrix.H13 - matrix.H13*matrix.H22*matrix.H31 - matrix.H32*matrix.H23*matrix.H11 - matrix.H12*matrix.H21;
	double inverse_matrix[3][3];
	//Adjugate matrix...
	inverse_matrix[0][0] = (matrix.H22 - matrix.H23*matrix.H32) / determ_of_matrix;
	inverse_matrix[1][0] = -(matrix.H21 - matrix.H23*matrix.H31) / determ_of_matrix;
	inverse_matrix[2][0] = (matrix.H21*matrix.H32 - matrix.H22*matrix.H31) / determ_of_matrix;
	inverse_matrix[0][1] = -(matrix.H12 - matrix.H32*matrix.H13) / determ_of_matrix;
	inverse_matrix[1][1] = (matrix.H11 - matrix.H13*matrix.H31) / determ_of_matrix;
	inverse_matrix[2][1] = -(matrix.H11*matrix.H32 - matrix.H12*matrix.H31) / determ_of_matrix;
	inverse_matrix[0][2] = (matrix.H12*matrix.H23 - matrix.H22*matrix.H13) / determ_of_matrix;
	inverse_matrix[1][2] = -(matrix.H11*matrix.H23 - matrix.H13*matrix.H21) / determ_of_matrix;
	inverse_matrix[2][2] = (matrix.H11*matrix.H22 - matrix.H12*matrix.H21) / determ_of_matrix;

	double remap_x = inverse_matrix[0][0] * test_p.x + inverse_matrix[0][1] * test_p.y + inverse_matrix[0][2];
	double remap_y = inverse_matrix[1][0] * test_p.x + inverse_matrix[1][1] * test_p.y + inverse_matrix[1][2];
	remap_x /= inverse_matrix[2][0] * test_p.x + inverse_matrix[2][1] * test_p.y + inverse_matrix[2][2];
	remap_y /= inverse_matrix[2][0] * test_p.x + inverse_matrix[2][1] * test_p.y + inverse_matrix[2][2];

	if (remap_x >= 0 && remap_x <= pre_width && remap_y >= 0 && remap_y <= pre_height)
		return 1;
	else
		return 0;
}

BMP* ImageBlending_rgb(BMP* Originalbmp, BMP* img1, BMP* img2, int img1x0, int img1y0, int img2x0, int img2y0, BMP* newimg, homographicmatrix mat)
{

	int img1wid = BMP_GetWidth(img1);
	int img1hei = BMP_GetHeight(img1);
	int img2wid = BMP_GetWidth(img2);
	int img2hei = BMP_GetHeight(img2);
	int depth = BMP_GetDepth(img2);
	int regx0 = max(img1x0, img2x0);
	int regwid = min(img1x0 + img1wid, img2x0 + img2wid) - regx0;
	int regy0 = max(img1y0, img2y0);
	int reghei = min(img1y0 + img1hei, img2y0 + img2hei) - regy0;
	//Create the overlap part of img1
	double* img1part_r = malloc(regwid*reghei * sizeof(double));
	double* img1part_g = malloc(regwid*reghei * sizeof(double));
	double* img1part_b = malloc(regwid*reghei * sizeof(double));
	memset(img1part_r, -1, sizeof(int)*regwid*reghei); // -1 means transparent
	memset(img1part_g, -1, sizeof(int)*regwid*reghei);
	memset(img1part_b, -1, sizeof(int)*regwid*reghei);
	for (int i = 0; i < regwid; i++)
	{
		for (int j = 0; j < reghei; j++)
		{
			UCHAR r = 0;
			UCHAR g = 0;
			UCHAR b = 0;
			point p = { 0,0 };
			if (judge_edge(Originalbmp, mat, i + regx0 - img1x0, j + regy0 - img1y0, img1wid, img1hei) == 1)
			{
				BMP_GetPixelRGB(img1, i + regx0 - img1x0, j + regy0 - img1y0, &r, &g, &b);
				img1part_r[i + j*regwid] = r;
				img1part_g[i + j*regwid] = g;
				img1part_b[i + j*regwid] = b;
			}
		}
	}
	//Create the overlap part of img2
	int* img2part_r = malloc(regwid*reghei * sizeof(int));
	int* img2part_g = malloc(regwid*reghei * sizeof(int));
	int* img2part_b = malloc(regwid*reghei * sizeof(int));
	for (int i = 0; i < regwid; i++)
	{
		for (int j = 0; j < reghei; j++)
		{
			UCHAR r = 0, g = 0, b = 0;
			BMP_GetPixelRGB(img2, regx0 - img2x0 + i, regy0 - img2y0 + j, &r, &g, &b);
			img2part_r[i + j*regwid] = r;
			img2part_g[i + j*regwid] = g;
			img2part_b[i + j*regwid] = b;
		}
	}

	//Reconstruct
	int* blendreg_r = malloc(sizeof(int)*reghei*regwid);
	int* blendreg_g = malloc(sizeof(int)*reghei*regwid);
	int* blendreg_b = malloc(sizeof(int)*reghei*regwid);
	for (int i = 0; i < regwid; i++)
	{
		double img1weight = 1 - (1.0*i / regwid);
		double img2weight = 1.0*i / regwid;
		for (int j = 0; j < reghei; j++)
		{
			//R;
			int i1 = img1part_r[i + j *regwid];
			int i2 = img2part_r[i + j*regwid];
			if (i1 >= 0)	blendreg_r[i + j*regwid] = i1*img1weight + i2*img2weight; //Linear blending
			else
				blendreg_r[i + j*regwid] = i2;
			//G;
			i1 = img1part_g[i + j*regwid];
			i2 = img2part_g[i + j*regwid];
			if (i1 >= 0)	blendreg_g[i + j*regwid] = i1*img1weight + i2*img2weight; //Linear blending
			else
				blendreg_g[i + j*regwid] = i2;
			//B;
			i1 = img1part_b[i + j*regwid];
			i2 = img2part_b[i + j*regwid];
			if (i1 >= 0)	blendreg_b[i + j*regwid] = i1*img1weight + i2*img2weight; //Linear blending
			else
				blendreg_b[i + j*regwid] = i2;
		}
	}
	//Write Newimg
	int newwid = BMP_GetWidth(newimg);
	int newhei = BMP_GetHeight(newimg);
	for (int i = 0; i < newwid; i++)
	{
		for (int j = 0; j < newhei; j++)
		{
			if (i >= regx0&&i < (regx0 + regwid) && j >= regy0&&j < (regy0 + reghei)) // if in overlap region
			{
				UCHAR i1 = blendreg_r[(i - regx0) + (j - regy0)*regwid];
				UCHAR i2 = blendreg_g[(i - regx0) + (j - regy0)*regwid];
				UCHAR i3 = blendreg_b[(i - regx0) + (j - regy0)*regwid];
				BMP_SetPixelRGB(newimg, (i), (j), i1, i2, i3);
			}
			else
			{
				if (i >= img1x0&&i < (img1x0 + img1wid) && j >= img1y0&&j < (img1y0 + img1hei)) // if in img1
				{
					UCHAR i1 = 0;
					UCHAR i2 = 0;
					UCHAR i3 = 0;
					BMP_GetPixelRGB(img1, i-img1x0, j-img1y0, &i1, &i2, &i3);
					BMP_SetPixelRGB(newimg, i, j, i1, i2, i3);
				}
				else
				{
					if (i >= img2x0&&i < (img2x0 + img2wid) && j >= img2y0&&j < (img2y0 + img2hei)) // if in img2
					{
						UCHAR r1 = 0, g1 = 0, b1 = 0;
						BMP_GetPixelRGB(img2, i - img2x0, j - img2y0, &r1, &g1, &b1);
						BMP_SetPixelRGB(newimg, i, j, r1, g1, b1);
					}
				}
			}
		}
	}
	//Freeeeeeee!!!
	free(img1part_r);
	free(img1part_g);
	free(img1part_b);
	free(img2part_r);
	free(img2part_g);
	free(img2part_b);
	/*free(img1high);
	free(img1low);
	free(img2high);
	free(img2low);
	free(blendhigh);
	free(blendlow);*/
	free(blendreg_r);
	free(blendreg_g);
	free(blendreg_b);
}