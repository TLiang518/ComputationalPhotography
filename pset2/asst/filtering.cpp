/* -----------------------------------------------------------------
 * File:    filtering.cpp
 * Created: 2015-09-22
 * -----------------------------------------------------------------
 * 
 * Image convolution and filtering
 * 
 * ---------------------------------------------------------------*/


#include "filtering.h"
#include <cmath>
#include <cassert>

using namespace std;

Image boxBlur(const Image &im, int k, bool clamp) {
    // --------- HANDOUT  PS02 ------------------------------
    // convolve an image with a box filter of size k by k
    
    /*
    In the following problems, you will implement several functions
    in which you will convolve an image with a kernel. This will
    require that you index out of the bounds of the image. Handle
    these boundary effects with the smart accessor from Section 2
    Also, process each of the three color channels independently.

    We have provided you with a function Image impulseImg(int k) 
    that generates a k × k × 1 grayscale image that is black everywhere
    except for one pixel in the center that is completely white. If you
    convolve a kernel with this image, you should get a copy of the kernel
    in the center of the image. An example of this can be seen in in a3 main.cpp.
    
    Implement the box filter Image boxBlur(const Image &im, int k, bool clamp=true) 
    in filtering.cpp . Each pixel of the output is the average of its k × k neighbors 
    in the input image, where k is an integer. Make sure the average is centered. 
    We will only test you on odd k.
    */
    Image output_image(im.width(), im.height(), im.channels());
    for (int i = 0; i < im.width(); i ++){
        for (int j = 0; j < im.height(); j++){
            for (int c = 0; c < im.channels(); c++){
                float avg = 0;
                for (int k_i_counter = 0; k_i_counter < k; k++){
                    for (int k_j_counter = 0; k_j_counter < k; k++){
                        avg += im.smartAccessor(i+k_i_counter, j+k_j_counter, c, clamp);
                    }
                }
                avg = avg / (k*k);
                output_image(i,j,c) = avg;
            }
        }
    }

}

Image Filter::convolve(const Image &im, bool clamp){
    // --------- HANDOUT  PS02 ------------------------------
    // Write a convolution function for the filter class
    return im; // change this
}



Image boxBlur_filterClass(const Image &im, int k, bool clamp) {
    // --------- HANDOUT  PS02 ------------------------------
    // Reimplement the box filter using the filter class.
    // check that your results match those in the previous function "boxBlur"
    return im; // change this
}


Image gradientMagnitude(const Image &im, bool clamp){
    // --------- HANDOUT  PS02 ------------------------------
    // Uses a Sobel kernel to compute the horizontal and vertical
    // components of the gradient of an image and returns the gradient magnitude.
    return im; // changeme
    
}

vector<float> gauss1DFilterValues(float sigma, float truncate){
    // --------- HANDOUT  PS02 ------------------------------
    // Create a vector containing the normalized values in a 1D Gaussian filter
    // Truncate the gaussian at truncate*sigma.
    return vector<float>();
}

Image gaussianBlur_horizontal(const Image &im, float sigma, float truncate, bool clamp){
    // --------- HANDOUT  PS02 ------------------------------
    // Gaussian blur across the rows of an image
    return im;
}

vector<float> gauss2DFilterValues(float sigma, float truncate){
    // --------- HANDOUT  PS02 ------------------------------
    // create a vector containing the normalized values in a 2D Gaussian
    // filter. Truncate the gaussian at truncate*sigma.
    return vector<float>();
}


Image gaussianBlur_2D(const Image &im, float sigma, float truncate, bool clamp){
    // --------- HANDOUT  PS02 ------------------------------
    //  Blur an image with a full  full 2D rotationally symmetric Gaussian kernel
    return im;
}

Image gaussianBlur_separable(const Image &im, float sigma, float truncate, bool clamp){
    // --------- HANDOUT  PS02 ------------------------------
    // Use principles of seperabiltity to blur an image using 2 1D Gaussian Filters
    return im;
}


Image unsharpMask(const Image &im, float sigma, float truncate, float strength, bool clamp){
    // --------- HANDOUT  PS02 ------------------------------
    // sharpen an image
    return im;
}


Image bilateral(const Image &im, float sigmaRange, float sigmaDomain, float truncateDomain, bool clamp){
    // --------- HANDOUT  PS02 ------------------------------
    // Denoise an image using the bilateral filter
    return im;
}


Image bilaYUV(const Image &im, float sigmaRange, float sigmaY, float sigmaUV, float truncateDomain, bool clamp){
    // --------- HANDOUT  PS02 ------------------------------
    // 6.865 only
    // Bilaterial Filter an image seperatly for
    // the Y and UV components of an image
    return im;
}




/**************************************************************
 //               DON'T EDIT BELOW THIS LINE                //
 *************************************************************/

// Create an image of 0's with a value of 1 in the middle. This function
// can be used to test that you have properly set the kernel values in your
// Filter object. Make sure to set k to be larger than the size of your kernel
Image impulseImg(int k){
    // initlize a kxkx1 image of all 0's
    Image impulse(k, k, 1);
    
    // set the center pixel to have intensity 1
    int center = floor(k/2);
    impulse(center,center,0) = 1.0f;
    
    return impulse;
}


// ------------- FILTER CLASS -----------------------
Filter::Filter(const vector<float> &fData, int fWidth, int fHeight) 
    : kernel(fData), width(fWidth), height(fHeight) 
{
        assert(fWidth*fHeight == fData.size());
}


Filter::Filter(int fWidth, int fHeight) 
    : kernel(std::vector<float>(fWidth*fHeight,0)), width(fWidth), height(fHeight) {} 


Filter::~Filter() {}


const float & Filter::operator()(int x, int y) const {
    if (x < 0 || x >= width)
        throw OutOfBoundsException();
    if ( y < 0 || y >= height)
        throw OutOfBoundsException();
    
    return kernel[x + y*width];
}


float & Filter::operator()(int x, int y) {
    if (x < 0 || x >= width)
        throw OutOfBoundsException();
    if ( y < 0 || y >= height)
        throw OutOfBoundsException();
    
    return kernel[x +y*width];
}
// --------- END FILTER CLASS -----------------------
