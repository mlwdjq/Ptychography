//-------------------------------------------------------------------------------
// 文件：unwrap-sr.c
//
// 描述：相位解包 -- 基于可靠性排序的二维相位解包 SR (Sorting Reliability)
//       参考文献 "Fast two-dimensional phase-unwrapping algorithm based on 
//       sorting by reliability following a noncontinuous path"
//
// 作者：李建欣
//
// 修改：2009年03月21日 第一版
//-------------------------------------------------------------------------------
//#include "stdafx.h"

#include "unwrap-sr.h"
#include <math.h>
#include "mex.h"
 
/* 输入变量 */
#define Phase	prhs[0]
#define Mask	prhs[1]

/* 输出变量 */
#define UP     plhs[0]

/* 宏定义 */
#define PI 3.14159265


/* 主函数接口 */
//-------------------------------------------------------------------------------
void mexFunction( int nlhs, mxArray *plhs[], 
          int nrhs, const mxArray*prhs[] )
{ 
    double *phase; 
	double *phase2;
    double *unwrapped_phase; 
	double *mask_double; 
	unsigned char * mask;
    int m,n,height,width; 
    
    /* 检查输入输出格式 */    
    if (nrhs != 2) { 
    mexErrMsgTxt("请输入一个参数自变量"); 
	} 
    else if (nlhs > 1) {
    mexErrMsgTxt("输出变量应该只有一个"); 
    } 

    /* 计算图片尺寸 */ 
    height = mxGetM(Phase); 
    width = mxGetN(Phase);

	
    /* 创建矩阵 */ 
    UP = mxCreateDoubleMatrix(height, width, mxREAL); 

    /* 给指针赋值 */ 
    phase = mxGetPr(Phase);
	mask_double = mxGetPr(Mask);
	unwrapped_phase = mxGetPr(UP);
//	phase = mxGetPr(UP);

	mask=new unsigned char [width*height]; 
	for(m=0;m<height;m++)
		for(n=0;n<width;n++)
			mask[m*width+n]=(unsigned char)mask_double[m*width+n];

	phase2=new double [width*height]; 
	memcpy(phase2,phase,width*height*sizeof(double));
        
    /* 计算 */
//	songle(i,iout,m,n);
	UnwrapPhaseBySortingReliabilityWithMask(phase2, -PI, PI, mask, 255, width, height);

	memcpy(unwrapped_phase,phase2,width*height*sizeof(double));
	delete [] mask;
	delete [] phase2;

    return;
}

//-------------------------------------------------------------------------------


static void UnwrapSRErrorShow(char *fmt, ...)
{
#if (!defined(__UNWRAP_SR_DLL_ENABLE__) && defined(__UNWRAP_SR_DEBUG_ENABLE__))
	va_list argptr;
	va_start(argptr, fmt);
	vprintf(fmt, argptr);
	va_end(argptr);
#endif
}


static int EdgeNo = 0;

static UNWRAP_SR_YES_NO FindPivot(UNWRAP_SR_EDGE *left, UNWRAP_SR_EDGE *right, double *PivotPtr)
{
	UNWRAP_SR_EDGE a, b, c, *p;
	
	a = *left;
	b = *(left + (right - left) / 2);
	c = *right;
	UNWRAP_SR_O3(a,b,c);
	
	if (a.reliab < b.reliab)
	{
		*PivotPtr = b.reliab;
		return UNWRAP_SR_YES;
	}
	
	if (b.reliab < c.reliab)
	{
		*PivotPtr = c.reliab;
		return UNWRAP_SR_YES;
	}
	
	for (p = left + 1; p <= right; ++p)
	{
		if (p->reliab != left->reliab)
		{
			*PivotPtr = (p->reliab < left->reliab) ? left->reliab : p->reliab;
			return UNWRAP_SR_YES;
		}

		return UNWRAP_SR_NO;
	}

	return UNWRAP_SR_NO;
}


static UNWRAP_SR_EDGE *partition(UNWRAP_SR_EDGE *left, UNWRAP_SR_EDGE *right, double pivot)
{
	while (left <= right)
	{
		while (left->reliab < pivot) ++left;
		while (right->reliab >= pivot) --right;

		if (left < right)
		{
			UNWRAP_SR_SWAP (*left, *right);
			++left;
			--right;
		}
	}

	return left;
}


static void QuickerSort(UNWRAP_SR_EDGE *left, UNWRAP_SR_EDGE *right)
{
	double pivot;
	UNWRAP_SR_EDGE *p;
	
	if (FindPivot(left, right, &pivot) == UNWRAP_SR_YES)
	{
		p = partition(left, right, pivot);
		QuickerSort(left, p - 1);
		QuickerSort(p, right);
	}
}


// initialse pixels. See the explanation of the pixel class above.
// initially every pixel is assumed to belong to a group consisiting of only itself
static void InitPIXEL(double *WrappedImage, unsigned char *InputMask, unsigned char *ExtendedMask, 
					  UNWRAP_SR_PIXELM *pixel, int ImageWidth, int ImageHeight)
{
	int i, j;
	double *WrappedImagePointer = WrappedImage;
	unsigned char *InputMaskPointer = InputMask;
	unsigned char *ExtendedMaskPointer = ExtendedMask;
	UNWRAP_SR_PIXELM *PixelPointer = pixel;
	
    for (i = 0; i < ImageHeight; i++)
	{
		for (j = 0; j < ImageWidth; j++)
        {
			PixelPointer->increment = 0;
			PixelPointer->PixelNumInGroup = 1;		
			PixelPointer->value = *WrappedImagePointer;
			PixelPointer->reliability = 9999999 + rand();
			PixelPointer->InputMask = *InputMaskPointer;
			PixelPointer->ExtendedMask = *ExtendedMaskPointer;
			PixelPointer->head = PixelPointer;
			PixelPointer->last = PixelPointer;
			PixelPointer->next = NULL;			
            PixelPointer->NewGroup = 0;
            PixelPointer->group = -1;
            PixelPointer++;
            WrappedImagePointer++;
			InputMaskPointer++;
			ExtendedMaskPointer++;
		}
	}
}


// gamma finction in the paper
static double wrap(double PixelValue)
{
	double WrappedPixelValue;
	if (PixelValue > UNWRAP_SR_PI)	WrappedPixelValue = PixelValue - UNWRAP_SR_TWO_PI;
	else if (PixelValue < -UNWRAP_SR_PI)	WrappedPixelValue = PixelValue + UNWRAP_SR_TWO_PI;
	else WrappedPixelValue = PixelValue;
	return WrappedPixelValue;
}


// PixelLValue is the left pixel,	PixelRValue is the right pixel
static int FindWrap(double PixelLValue, double PixelRValue)
{
	int WrapValue;
	double difference; 
	difference = PixelLValue - PixelRValue;	
	if (difference > UNWRAP_SR_PI) WrapValue = -1;
	else if (difference < -UNWRAP_SR_PI) WrapValue = 1;
	else WrapValue = 0;	
	return WrapValue;
} 


static void ExtendMask(unsigned char *InputMask, unsigned char *ExtendedMask, int ImageWidth, int ImageHeight)
{
	int i,j;
	int ImageWidthPlusOne = ImageWidth + 1;
	int ImageWidthMinusOne = ImageWidth - 1;
	unsigned char *IMP = InputMask    + ImageWidth + 1;	//input mask pointer
	unsigned char *EMP = ExtendedMask + ImageWidth + 1;	//extended mask pointer
	
	// extend the mask for the image except borders
	for (i = 1; i < ImageHeight - 1; ++i)
	{
		for (j = 1; j < ImageWidth - 1; ++j)
		{
			if ((*IMP) == 255 && (*(IMP + 1) == 255) && (*(IMP - 1) == 255) && 
				(*(IMP + ImageWidth) == 255) && (*(IMP - ImageWidth) == 255) &&
				(*(IMP - ImageWidthMinusOne) == 255) && (*(IMP - ImageWidthPlusOne) == 255) &&
				(*(IMP + ImageWidthMinusOne) == 255) && (*(IMP + ImageWidthPlusOne) == 255))
			{		
				*EMP = 255;
			}

			++EMP;
			++IMP;
		}

		EMP += 2;
		IMP += 2;
	}
}


static void CalcReliability(double *WrappedImage, UNWRAP_SR_PIXELM *pixel, int ImageWidth, int ImageHeight)
{
	int i, j;
	int ImageWidthPlusOne = ImageWidth + 1;
	int ImageWidthMinusOne = ImageWidth - 1;
	double H, V, D1, D2;
	double *WIP = WrappedImage + ImageWidthPlusOne; // WIP is the wrapped image pointer
	UNWRAP_SR_PIXELM *PixelPointer = pixel + ImageWidthPlusOne;
	
	for (i = 1; i < ImageHeight -1; ++i)
	{
		for (j = 1; j < ImageWidth - 1; ++j)
		{
			if (PixelPointer->InputMask == 255)
			{
				H = wrap(*(WIP - 1) - *WIP) - wrap(*WIP - *(WIP + 1));
				V = wrap(*(WIP - ImageWidth) - *WIP) - wrap(*WIP - *(WIP + ImageWidth));
				D1 = wrap(*(WIP - ImageWidthPlusOne) - *WIP) - wrap(*WIP - *(WIP + ImageWidthPlusOne));
				D2 = wrap(*(WIP - ImageWidthMinusOne) - *WIP) - wrap(*WIP - *(WIP + ImageWidthMinusOne));
				PixelPointer->reliability = H*H + V*V + D1*D1 + D2*D2;
			}

			PixelPointer++;
			WIP++;
		}

		PixelPointer += 2;
		WIP += 2;
	}
}


// calculate the reliability of the horizontal edges of the image
// it is calculated by adding the reliability of pixel and the relibility of 
// its right-hand neighbour
// edge is calculated between a pixel and its next neighbour
static void HorizontalEDGE(UNWRAP_SR_PIXELM *pixel, UNWRAP_SR_EDGE *edge, int ImageWidth, int ImageHeight)
{
	int i, j;
	UNWRAP_SR_EDGE *EdgePointer = edge;
	UNWRAP_SR_PIXELM *PixelPointer = pixel;
	
	for (i = 0; i < ImageHeight; i++)
	{
		for (j = 0; j < ImageWidth - 1; j++) 
		{
			if (PixelPointer->InputMask == 255 && (PixelPointer + 1)->InputMask == 255)
			{
				EdgePointer->pointer1 = PixelPointer;
				EdgePointer->pointer2 = (PixelPointer+1);
				EdgePointer->reliab = PixelPointer->reliability + (PixelPointer + 1)->reliability;
				EdgePointer->increment = FindWrap(PixelPointer->value, (PixelPointer + 1)->value);
				EdgePointer++;
				EdgeNo++;
			}

			PixelPointer++;
		}

		PixelPointer++;
	}
}


// calculate the reliability of the vertical edges of the image
// it is calculated by adding the reliability of pixel and the relibility of 
// its lower neighbour in the image.
static void VerticalEDGE(UNWRAP_SR_PIXELM *pixel, UNWRAP_SR_EDGE *edge, int ImageWidth, int ImageHeight)
{
	int i, j;
	UNWRAP_SR_PIXELM *PixelPointer = pixel;
	UNWRAP_SR_EDGE *EdgePointer = edge + EdgeNo; 
	
	for (i = 0; i < ImageHeight - 1; i++)
	{
		for (j = 0; j < ImageWidth; j++) 
		{
			if (PixelPointer->InputMask == 255 && (PixelPointer + ImageWidth)->InputMask == 255)
			{
				EdgePointer->pointer1 = PixelPointer;
				EdgePointer->pointer2 = (PixelPointer + ImageWidth);
				EdgePointer->reliab = PixelPointer->reliability + (PixelPointer + ImageWidth)->reliability;
				EdgePointer->increment = FindWrap(PixelPointer->value, (PixelPointer + ImageWidth)->value);
				EdgePointer++;
				EdgeNo++;
			}

			PixelPointer++;
		} 
	} 
}


// gather the pixels of the image into groups 
static void GatherPIXEL(UNWRAP_SR_EDGE *edge, int ImageWidth, int ImageHeight)
{
	int k, incremento;
	UNWRAP_SR_PIXELM *PIXEL1, *PIXEL2, *group1, *group2;   
	UNWRAP_SR_EDGE *PointerEdge = edge;
	
	for (k = 0; k < EdgeNo; k++)
	{
		PIXEL1 = PointerEdge->pointer1;
		PIXEL2 = PointerEdge->pointer2;
		
		// UNWRAP_SR_PIXELM 1 and UNWRAP_SR_PIXELM 2 belong to different groups
		// initially each pixel is a group by it self and one pixel can construct a group
		// UNWRAP_SR_NO else or else if to this if
		if (PIXEL2->head != PIXEL1->head)
		{
			// UNWRAP_SR_PIXELM 2 is alone in its group
			// merge this pixel with UNWRAP_SR_PIXELM 1 group and find the number of 2 pi to add 
			// to or subtract to unwrap it
			if ((PIXEL2->next == NULL) && (PIXEL2->head == PIXEL2))
			{
				PIXEL1->head->last->next = PIXEL2;
				PIXEL1->head->last = PIXEL2;
				(PIXEL1->head->PixelNumInGroup)++;
				PIXEL2->head=PIXEL1->head;
				PIXEL2->increment = PIXEL1->increment - PointerEdge->increment;
			}			
			// UNWRAP_SR_PIXELM 1 is alone in its group
			// merge this pixel with UNWRAP_SR_PIXELM 2 group and find the number of 2 pi to add 
			// to or subtract to unwrap it
			else if ((PIXEL1->next == NULL) && (PIXEL1->head == PIXEL1))
			{
				PIXEL2->head->last->next = PIXEL1;
				PIXEL2->head->last = PIXEL1;
				(PIXEL2->head->PixelNumInGroup)++;
				PIXEL1->head = PIXEL2->head;
				PIXEL1->increment = PIXEL2->increment + PointerEdge->increment;
			} 			
			// UNWRAP_SR_PIXELM 1 and UNWRAP_SR_PIXELM 2 both have groups
			else
            {
				group1 = PIXEL1->head;
                group2 = PIXEL2->head;

				// if the UNWRAP_SR_NO. of pixels in UNWRAP_SR_PIXELM 1 group is larger than the UNWRAP_SR_NO. of pixels
				// in UNWRAP_SR_PIXELM 2 group.   Merge UNWRAP_SR_PIXELM 2 group to UNWRAP_SR_PIXELM 1 group
				// and find the number of wraps between UNWRAP_SR_PIXELM 2 group and UNWRAP_SR_PIXELM 1 group
				// to unwrap UNWRAP_SR_PIXELM 2 group with respect to UNWRAP_SR_PIXELM 1 group.
				// the UNWRAP_SR_NO. of wraps will be added to UNWRAP_SR_PIXELM 2 group in the future
				if (group1->PixelNumInGroup > group2->PixelNumInGroup)
				{
					// merge UNWRAP_SR_PIXELM 2 with UNWRAP_SR_PIXELM 1 group
					group1->last->next = group2;
					group1->last = group2->last;
					group1->PixelNumInGroup = group1->PixelNumInGroup + group2->PixelNumInGroup;
					incremento = PIXEL1->increment-PointerEdge->increment - PIXEL2->increment;

					// merge the other pixels in UNWRAP_SR_PIXELM 2 group to UNWRAP_SR_PIXELM 1 group
					while (group2 != NULL)
					{
						group2->head = group1;
						group2->increment += incremento;
						group2 = group2->next;
					}
				} 
				// if the UNWRAP_SR_NO. of pixels in UNWRAP_SR_PIXELM 2 group is larger than the UNWRAP_SR_NO. of pixels
				// in UNWRAP_SR_PIXELM 1 group.   Merge UNWRAP_SR_PIXELM 1 group to UNWRAP_SR_PIXELM 2 group
				// and find the number of wraps between UNWRAP_SR_PIXELM 2 group and UNWRAP_SR_PIXELM 1 group
				// to unwrap UNWRAP_SR_PIXELM 1 group with respect to UNWRAP_SR_PIXELM 2 group.
				// the UNWRAP_SR_NO. of wraps will be added to UNWRAP_SR_PIXELM 1 grop in the future
				else
                {
					// merge UNWRAP_SR_PIXELM 1 with UNWRAP_SR_PIXELM 2 group
					group2->last->next = group1;
					group2->last = group1->last;
					group2->PixelNumInGroup = group2->PixelNumInGroup + group1->PixelNumInGroup;
					incremento = PIXEL2->increment + PointerEdge->increment - PIXEL1->increment;

					// merge the other pixels in UNWRAP_SR_PIXELM 2 group to UNWRAP_SR_PIXELM 1 group
					while (group1 != NULL)
					{
						group1->head = group2;
						group1->increment += incremento;
						group1 = group1->next;
					} // while
                } // else
            } // else
        } // if

        PointerEdge++;
	}
}


// unwrap the image 
static void UnwrapImage(UNWRAP_SR_PIXELM *pixel, int ImageWidth, int ImageHeight)
{
	int i, ImageSize = ImageWidth * ImageHeight;
	UNWRAP_SR_PIXELM *PixelPointer = pixel;
	
	for (i = 0; i < ImageSize; i++)
	{
		PixelPointer->value += UNWRAP_SR_TWO_PI * (double)(PixelPointer->increment);
        PixelPointer++;
    }
}


// set the masked pixels (mask = 0) to the minimum of the unwrapper phase
static void MaskImage(UNWRAP_SR_PIXELM *pixel, unsigned char *InputMask, int ImageWidth, int ImageHeight)
{
	int i;
	int ImageWidthPlusOne  = ImageWidth + 1;
	int ImageHeightPlusOne  = ImageHeight + 1;
	int ImageWidthMinusOne = ImageWidth - 1;
	int ImageHeightMinusOne = ImageHeight - 1;
	int ImageSize = ImageWidth * ImageHeight;
	double min = 99999999.;
	unsigned char *IMP = InputMask;	//input mask pointer
	UNWRAP_SR_PIXELM *PointerPixel = pixel;
	
	//find the minimum of the unwrapped phase
	for (i = 0; i < ImageSize; i++)
	{
		if ((PointerPixel->value < min) && (*IMP == 255)) 
		{
			min = PointerPixel->value;
		}
		
		PointerPixel++;
		IMP++;
	}
	
	PointerPixel = pixel;
	IMP = InputMask;	
	
	//set the masked pixels to minimum
	for (i = 0; i < ImageSize; i++)
	{
		if ((*IMP) == 0)
		{
			PointerPixel->value = min;
		}

		PointerPixel++;
		IMP++;
	}
}


// the input to this unwrapper is an array that contains the wrapped phase map. 
// copy the image on the buffer passed to this unwrapper to over-write the unwrapped 
// phase map on the buffer of the wrapped phase map.
static void ReturnImage(UNWRAP_SR_PIXELM *pixel, double *UnwrappedImage, int ImageWidth, int ImageHeight)
{
	int i, ImageSize = ImageWidth * ImageHeight;
    double *UnwrappedImagePointer = UnwrappedImage;
    UNWRAP_SR_PIXELM *PixelPointer = pixel;
	
    for (i = 0; i < ImageSize; i++) 
	{
        *UnwrappedImagePointer = PixelPointer->value;
        PixelPointer++;
		UnwrappedImagePointer++;
	}
}


int UnwrapPhaseBySortingReliabilityWithMask(double *phase, double LowVal, double HighVal, 
											unsigned char *mask, int ValidPoint, int height, 
											int width)
{
	int i, j, ExectResult;
	int ImageSize = height * width;
	int InitEdgeNum = width * (height - 1) + (width - 1) * height;
	unsigned char *ExtendedMask = NULL;
	UNWRAP_SR_PIXELM *pixel = NULL;
	UNWRAP_SR_EDGE *edge = NULL;

	EdgeNo = 0;
	
	if ((ExtendedMask = (unsigned char *) calloc (ImageSize, sizeof(unsigned char))) == NULL)
	{
		UnwrapSRErrorShow("UnwrapPhaseBySortingReliabilityWithMask(): fail to allocate ExtendedMask.\n");
		ExectResult = UNWRAP_SR_ERROR_MEMORY;
		goto CleanReturn;
	}

	if ((pixel = (UNWRAP_SR_PIXELM *) calloc (ImageSize, sizeof(UNWRAP_SR_PIXELM))) == NULL)
	{
		UnwrapSRErrorShow("UnwrapPhaseBySortingReliabilityWithMask(): fail to allocate pixel.\n");
		ExectResult = UNWRAP_SR_ERROR_MEMORY;
		goto CleanReturn;
	}

	if ((edge = (UNWRAP_SR_EDGE *) calloc (InitEdgeNum, sizeof(UNWRAP_SR_EDGE))) == NULL)
	{
		UnwrapSRErrorShow("UnwrapPhaseBySortingReliabilityWithMask(): fail to allocate edge.\n");
		ExectResult = UNWRAP_SR_ERROR_MEMORY;
		goto CleanReturn;
	}

	// phase must be normalized to range: -PI - +PI
	// mask must be valid: 255, invalid: 0
	// When the mask is 255 this means that the pixel is valid 
	// When the mask is 0 this means that the pixel is invalid (noisy or corrupted pixel)
	for (j = 0; j < height; j++) 
	{
		for (i = 0; i < width; i++)
		{
			phase[j * width + i] = -UNWRAP_SR_PI + UNWRAP_SR_TWO_PI *
				((phase[j * width + i] - LowVal) / (HighVal - LowVal));
			mask[j * width + i] = mask[j * width + i] == ValidPoint ? 255 : 0;
		}
	}

	ExtendMask(mask, ExtendedMask, width, height);	
	InitPIXEL(phase, mask, ExtendedMask, pixel, width, height);	
	CalcReliability(phase, pixel, width, height);	
	HorizontalEDGE(pixel, edge, width, height);	
	VerticalEDGE(pixel, edge, width, height);	

	// sort the EDGEs depending on their reiability. The PIXELs with higher relibility (small value) first
	QuickerSort(edge, edge + EdgeNo - 1);	

	// gather PIXELs into groups
	GatherPIXEL(edge, width, height);	
	UnwrapImage(pixel, width, height);	
	MaskImage(pixel, mask, width, height);
	
	// copy the image from UNWRAP_SR_PIXELM structure to the unwrapped phase array passed to this function
	ReturnImage(pixel, phase, width, height);
	ExectResult = UNWRAP_SR_SUCCESS;
	
CleanReturn:
	if (edge) free(edge);
	if (pixel) free(pixel);
	if (ExtendedMask) free(ExtendedMask);
	return ExectResult;
}


#if (defined(__UNWRAP_SR_DLL_ENABLE__) && defined(WIN32))
#include <windows.h>

BOOL APIENTRY DllMain(HINSTANCE hInstance, DWORD dwReason, LPVOID lpReserved)
{
    switch (dwReason)
	{
	case DLL_PROCESS_ATTACH: break;
	case DLL_THREAD_ATTACH:
	case DLL_THREAD_DETACH: break;
	case DLL_PROCESS_DETACH: break;
	}
	
    return TRUE;
}
#endif

