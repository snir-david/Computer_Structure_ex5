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
    int currRow, currCol;
    pixel_sum sum;
    pixel current_pixel;
    int min_intensity = 766; // arbitrary value that is higher than maximum possible intensity, which is 255*3=765
    int max_intensity = -1; // arbitrary value that is lower than minimum possible intensity, which is 0
    int min_row, min_col, max_row, max_col;
    pixel loop_pixel;

    initialize_pixel_sum(&sum);
    int maxI = max(i - 1, 0), maxJ = max(j - 1, 0), minI = min(i + 1, dim - 1), minJ = min(j + 1, dim - 1);
    for (ii = maxI; ii <= minI; ii++) {
        int iin = ii * dim;
        for (jj = maxJ; jj <= minJ; jj++) {

            int kRow, kCol;

            // compute row index in kernel
            if (ii < i) {
                kRow = 0;
            } else if (ii > i) {
                kRow = 2;
            } else {
                kRow = 1;
            }

            // compute column index in kernel
            if (jj < j) {
                kCol = 0;
            } else if (jj > j) {
                kCol = 2;
            } else {
                kCol = 1;
            }

            // apply kernel on pixel at [ii,jj]
            sum_pixels_by_weight(&sum, src[iin + jj], kernel[kRow][kCol]);
        }
    }

    if (filter) {
        // find min and max coordinates
        for (ii = maxI; ii <= minI; ii++) {
            int iin = ii * dim;
            for (jj = maxJ; jj <= minJ; jj++) {
                // check if smaller than min or higher than max and update
                loop_pixel = src[iin + jj];
                if ((((int) loop_pixel.red) + ((int) loop_pixel.green) + ((int) loop_pixel.blue)) <= min_intensity) {
                    min_intensity = (((int) loop_pixel.red) + ((int) loop_pixel.green) + ((int) loop_pixel.blue));
                    min_row = ii;
                    min_col = jj;
                }
                if ((((int) loop_pixel.red) + ((int) loop_pixel.green) + ((int) loop_pixel.blue)) > max_intensity) {
                    max_intensity = (((int) loop_pixel.red) + ((int) loop_pixel.green) + ((int) loop_pixel.blue));
                    max_row = ii;
                    max_col = jj;
                }
            }
        }
        // filter out min and max
        int currentIdx = min_row * dim + min_col;
        sum_pixels_by_weight(&sum, src[currentIdx], -1);
        sum_pixels_by_weight(&sum, src[currentIdx], -1);
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
    int halfKernel = kernelSize / 2;
    int dimMinKer = dim - halfKernel;
    for (i = halfKernel; i < dimMinKer; i++) {
        int in = i * dim;
        //TODO loop unrolling
        for (j = halfKernel; j < dimMinKer; j++) {
            dst[in + j] = applyKernel(dim, i, j, src, kernelSize, kernel, kernelScale, filter);
        }
    }
}

void charsToPixels(Image *charsImg, pixel *pixels) {
    int row, col;

    for (row = 0; row < m; row++) {
        int rowN = row * n, threeRowN = 3 * rowN;
        for (col = 0; col < n; col++) {
            int currentIdx = rowN + col, threeCol = 3 * col;
            pixels[currentIdx].red = image->data[threeRowN + threeCol];
            pixels[currentIdx].green = image->data[threeRowN + threeCol + 1];
            pixels[currentIdx].blue = image->data[threeRowN + threeCol + 2];
        }
    }
}

void pixelsToChars(pixel *pixels, Image *charsImg) {

    int row, col;

    for (row = 0; row < m; row++) {
        int rowN = row * n, threeRowN = 3 * rowN;
        for (col = 0; col < n; col++) {
            int currentIdx = rowN + col, threeCol = 3 * col;
            image->data[threeRowN + threeCol] = pixels[currentIdx].red;
            image->data[threeRowN + threeCol + 1] = pixels[currentIdx].green;
            image->data[threeRowN + threeCol + 2] = pixels[currentIdx].blue;
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

    pixel *pixelsImg = malloc(m * n * sizeof(pixel));
    pixel *backupOrg = malloc(m * n * sizeof(pixel));

    charsToPixels(image, pixelsImg);
    copyPixels(pixelsImg, backupOrg);

    smooth(m, backupOrg, pixelsImg, kernelSize, kernel, kernelScale, filter);

    pixelsToChars(pixelsImg, image);

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

