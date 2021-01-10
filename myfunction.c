#include <stdbool.h>

typedef struct {
    unsigned char red;
    unsigned char green;
    unsigned char blue;
} pixel;

typedef struct {
    int red;
    int green;
    int blue;
    // int num;
} pixel_sum;


/* Compute min and max of two integers, respectively */
int min(int a, int b) { return (a < b ? a : b); }

int max(int a, int b) { return (a > b ? a : b); }

int calcIndex(int i, int j, int n) {
    return ((i) * (n) + (j));
}

/*
 * initialize_pixel_sum - Initializes all fields of sum to 0
 */
void initialize_pixel_sum(pixel_sum *sum) {
    sum->red = sum->green = sum->blue = 0;
    // sum->num = 0;
    return;
}

/*
 * assign_sum_to_pixel - Truncates pixel's new value to match the range [0,255]
 */
static void assign_sum_to_pixel(pixel *current_pixel, pixel_sum sum, int kernelScale) {

    // divide by kernel's weight
    sum.red = sum.red / kernelScale;
    sum.green = sum.green / kernelScale;
    sum.blue = sum.blue / kernelScale;

    // truncate each pixel's color values to match the range [0,255]
    current_pixel->red = (unsigned char) (min(max(sum.red, 0), 255));
    current_pixel->green = (unsigned char) (min(max(sum.green, 0), 255));
    current_pixel->blue = (unsigned char) (min(max(sum.blue, 0), 255));
    return;
}

/*
* sum_pixels_by_weight - Sums pixel values, scaled by given weight
*/
static void sum_pixels_by_weight(pixel_sum *sum, pixel p, int weight) {
    sum->red += ((int) p.red) * weight;
    sum->green += ((int) p.green) * weight;
    sum->blue += ((int) p.blue) * weight;
    // sum->num++;
    return;
}

/*
 *  Applies kernel for pixel at (i,j)
 */
static pixel
applyKernel(int dim, int i, int j, pixel *src, int kernelSize, int kernel[kernelSize][kernelSize], int kernelScale,
            bool filter) {

    int ii, jj;
    pixel_sum sum;
    pixel current_pixel;
    int min_intensity = 766; // arbitrary value that is higher than maximum possible intensity, which is 255*3=765
    int max_intensity = -1; // arbitrary value that is lower than minimum possible intensity, which is 0
    int min_row, min_col, max_row, max_col;
    pixel loop_pixel;
    int pixel_loopSum;

    //initialize_pixel_sum(&sum);
    sum.red = sum.green = sum.blue = 0;
    int minI = min(i + 1, dim - 1), maxI = max(i - 1, 0), minJ = min(j + 1, dim - 1), maxJ = max(j - 1, 0);
    ii = maxI;
    int iin = ii * dim;
    for (; ii <= minI; ii++) {
        int kRow = ii - i + 1;
        // compute row index in kernel
        int kCol = maxJ - j + 1;
        int weight = kernel[kRow][kCol], currentIdx = iin + maxJ;
        // apply kernel on pixel at [ii,jj-jj+2]
        sum.blue += src[currentIdx].blue * weight;
        sum.red += src[currentIdx].red * weight;
        sum.green += src[currentIdx].green * weight;

        weight = kernel[kRow][kCol + 1];
        sum.blue += src[currentIdx + 1].blue * weight;
        sum.red += src[currentIdx + 1].red * weight;
        sum.green += src[currentIdx + 1].green * weight;

        weight = kernel[kRow][kCol + 2];
        sum.blue += src[currentIdx + 2].blue * weight;
        sum.red += src[currentIdx + 2].red * weight;
        sum.green += src[currentIdx + 2].green * weight;
        iin+=dim;
    }

    if (filter) {
        ii = maxI, iin = ii* dim, jj = maxJ;
        // find min and max coordinates
        for (; ii <= minI; ii++) {
                // check if smaller than min or higher than max and update
                loop_pixel = src[iin + jj];
                pixel_loopSum = ((int) loop_pixel.red) + ((int) loop_pixel.green) + ((int) loop_pixel.blue);
                if (pixel_loopSum <= min_intensity) {
                    min_intensity = pixel_loopSum;
                    min_row = ii;
                    min_col = jj;
                }
                if (pixel_loopSum > max_intensity) {
                    max_intensity = pixel_loopSum;
                    max_row = ii;
                    max_col = jj;
                }

                loop_pixel = src[iin + jj+1];
                pixel_loopSum = ((int) loop_pixel.red) + ((int) loop_pixel.green) + ((int) loop_pixel.blue);
                if (pixel_loopSum <= min_intensity) {
                    min_intensity = pixel_loopSum;
                    min_row = ii;
                    min_col = jj+1;
                }
                if (pixel_loopSum > max_intensity) {
                    max_intensity = pixel_loopSum;
                    max_row = ii;
                    max_col = jj+1;
                }

                loop_pixel = src[iin + jj+2];
                pixel_loopSum = ((int) loop_pixel.red) + ((int) loop_pixel.green) + ((int) loop_pixel.blue);
                if (pixel_loopSum <= min_intensity) {
                    min_intensity = pixel_loopSum;
                    min_row = ii;
                    min_col = jj+2;
                }
                if (pixel_loopSum > max_intensity) {
                    max_intensity = pixel_loopSum;
                    max_row = ii;
                    max_col = jj+2;
                }
            iin+=dim;
        }

        // filter out min and max
        int minIdx = min_row * dim + min_col, maxIdx = max_row * dim + max_col;
        sum.blue -= src[minIdx].blue;
        sum.red -= src[minIdx].red;
        sum.green -= src[minIdx].green;

        sum.blue -= src[maxIdx].blue;
        sum.red -= src[maxIdx].red;
        sum.green -= src[maxIdx].green;

    }

    // assign kernel's result to pixel at [i,j]
    assign_sum_to_pixel(&current_pixel, sum, kernelScale);
    return current_pixel;
}

/*
* Apply the kernel over each pixel.
* Ignore pixels where the kernel exceeds bounds. These are pixels with row index smaller than kernelSize/2 and/or
* column index smaller than kernelSize/2
*/
void smooth(int dim, pixel *src, pixel *dst, int kernelSize, int kernel[kernelSize][kernelSize], int kernelScale,
            bool filter) {
    int i, j;
    for (i = kernelSize / 2; i < dim - kernelSize / 2; i++) {
        int iDim = i * dim, limit = dim - kernelSize / 2 - 2;
        for (j = kernelSize / 2; j < limit; j += 3) {
            dst[iDim + j] = applyKernel(dim, i, j, src, kernelSize, kernel, kernelScale, filter);
            dst[iDim + j + 1] = applyKernel(dim, i, (j + 1), src, kernelSize, kernel, kernelScale, filter);
            dst[iDim + j + 2] = applyKernel(dim, i, (j + 2), src, kernelSize, kernel, kernelScale, filter);
        }
        for (; j < dim - kernelSize / 2; j++) {
            dst[iDim + j] = applyKernel(dim, i, j, src, kernelSize, kernel, kernelScale, filter);
        }
    }
}

void charsToPixels(pixel *pixels) {

    int row, col, limit = n - 2;
    for (row = 0; row < m; row++) {
        int rowN = row * n;
        for (col = 0; col < limit; col += 3) {
            int currentIdx = rowN + col, threeIdx = 3 * rowN + 3 * col;
            pixels[currentIdx].red = image->data[threeIdx];
            pixels[currentIdx].green = image->data[threeIdx + 1];
            pixels[currentIdx].blue = image->data[threeIdx + 2];

            pixels[currentIdx + 1].red = image->data[threeIdx + 3];
            pixels[currentIdx + 1].green = image->data[threeIdx + 4];
            pixels[currentIdx + 1].blue = image->data[threeIdx + 5];

            pixels[currentIdx + 2].red = image->data[threeIdx + 6];
            pixels[currentIdx + 2].green = image->data[threeIdx + 7];
            pixels[currentIdx + 2].blue = image->data[threeIdx + 8];
        }
        for (; col < n; col++) {
            int currentIdx = rowN + col, threeIdx = 3 * rowN + 3 * col;
            pixels[currentIdx].red = image->data[threeIdx];
            pixels[currentIdx].green = image->data[threeIdx + 1];
            pixels[currentIdx].blue = image->data[threeIdx + 2];
        }
    }
}

void pixelsToChars(pixel *pixels) {

    int row, col, limit = n - 2;
    for (row = 0; row < m; row++) {
        int rowN = row * n;
        for (col = 0; col < limit; col += 3) {
            int currentIdx = rowN + col, threeIdx = 3 * rowN + 3 * col;
            image->data[threeIdx] = pixels[currentIdx].red;
            image->data[threeIdx + 1] = pixels[currentIdx].green;
            image->data[threeIdx + 2] = pixels[currentIdx].blue;

            image->data[threeIdx + 3] = pixels[currentIdx + 1].red;
            image->data[threeIdx + 4] = pixels[currentIdx + 1].green;
            image->data[threeIdx + 5] = pixels[currentIdx + 1].blue;

            image->data[threeIdx + 6] = pixels[currentIdx + 2].red;
            image->data[threeIdx + 7] = pixels[currentIdx + 2].green;
            image->data[threeIdx + 8] = pixels[currentIdx + 2].blue;
        }
        for (; col < n; col++) {
            int currentIdx = rowN + col, threeIdx = 3 * rowN + 3 * col;
            image->data[threeIdx] = pixels[currentIdx].red;
            image->data[threeIdx + 1] = pixels[currentIdx].green;
            image->data[threeIdx + 2] = pixels[currentIdx].blue;

        }
    }
}

void copyPixels(pixel *src, pixel *dst) {

    int row, col;
    for (row = 0; row < m; row++) {
        int rowN = row * n;
        for (col = 0; col < n; col++) {
            int currentIdx = rowN + col;
            dst[currentIdx].red = src[currentIdx].red;
            dst[currentIdx].green = src[currentIdx].green;
            dst[currentIdx].blue = src[currentIdx].blue;
        }
    }
}

void doConvolution(Image *image, int kernelSize, int kernel[kernelSize][kernelSize], int kernelScale, bool filter) {
    int mallocSize = m * n * sizeof(pixel);
    pixel *pixelsImg = malloc(mallocSize);
    pixel *backupOrg = malloc(mallocSize);

    charsToPixels(pixelsImg);
    copyPixels(pixelsImg, backupOrg);

    smooth(m, backupOrg, pixelsImg, kernelSize, kernel, kernelScale, filter);

    pixelsToChars(pixelsImg);

    free(pixelsImg);
    free(backupOrg);
}

void myfunction(Image *image, char *srcImgpName, char *blurRsltImgName, char *sharpRsltImgName,
                char *filteredBlurRsltImgName, char *filteredSharpRsltImgName, char flag) {

    /*
    * [1, 1, 1]
    * [1, 1, 1]
    * [1, 1, 1]
    */
    int blurKernel[3][3] = {{1, 1, 1},
                            {1, 1, 1},
                            {1, 1, 1}};

    /*
    * [-1, -1, -1]
    * [-1, 9, -1]
    * [-1, -1, -1]
    */
    int sharpKernel[3][3] = {{-1, -1, -1},
                             {-1, 9,  -1},
                             {-1, -1, -1}};

    if (flag == '1') {
        // blur image
        doConvolution(image, 3, blurKernel, 9, false);

        // write result image to file
        writeBMP(image, srcImgpName, blurRsltImgName);

        // sharpen the resulting image
        doConvolution(image, 3, sharpKernel, 1, false);

        // write result image to file
        writeBMP(image, srcImgpName, sharpRsltImgName);
    } else {
        // apply extermum filtered kernel to blur image
        doConvolution(image, 3, blurKernel, 7, true);

        // write result image to file
        writeBMP(image, srcImgpName, filteredBlurRsltImgName);

        // sharpen the resulting image
        doConvolution(image, 3, sharpKernel, 1, false);

        // write result image to file
        writeBMP(image, srcImgpName, filteredSharpRsltImgName);
    }
}

