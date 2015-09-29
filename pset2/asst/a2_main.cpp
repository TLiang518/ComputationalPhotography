/* --------------------------------------------------------------------------
 * File:    a2_main.cpp
 * Created: 2015-09-23
 * --------------------------------------------------------------------------
 * 
 * 
 * 
 * ------------------------------------------------------------------------*/


#include "Image.h"
#include "filtering.h"
#include <ctime>
#include <iostream>
#include <vector>

using namespace std;


// This is a way for you to test your functions. 
// We will only grade the contents of filter.cpp and Image.cpp
int main() {

    // ------- Example tests, change them ! --------------
    Image im = impulseImg(10);
    //Go back and check this. I'm not sure what the behavior should be, even.
    
    Image white_square("./Input/white_square.png");
    //300x300
    Image cambridge1("./Input/Cambridge1.png");
    
    Image blurred = boxBlur(cambridge1, 3, true);
    blurred.write("./Output/boxBlurCambridge.png");    
    Image filter_blurred = boxBlur_filterClass(cambridge1, 3, true);
    filter_blurred.write("./Output/boxBlurFilterCambridge.png");

    cout << "blurred impulse image" << endl;
    cout << "keep testing..." << endl;
    
    // ---------------------------------------------------

    
    // ---------------------------------------------------
    // Test the filter class on an impulse image
    
    Image dirac = impulseImg(31);
    

    Image lounge_view("./Input/lounge_view.png");

    Image sobel = gradientMagnitude(lounge_view, true);
    sobel.write("./Output/sobel.png");

    // Test kernel
    vector<float> kernel{0,0,1,
                         0,1,0,
                         1,0,0}; // C++11 syntax
    Filter testFilter(kernel, 3, 3);
    Image testOutput = testFilter.convolve(dirac);
    // The output should be an exact copy of the kernel in the center of the
    // image
    testOutput.write("./Output/testKernel.png");
    // ---------------------------------------------------
    

    // ---------------------------------------------------
    // E.g. test the sobelKernel
    // create Sobel Filter that extracts horizontal gradients
    // [ -1 0 1 ]
    // [ -2 0 2 ]
    // [ -1 0 1 ]
    float fDataXArray[] = { -1.0, 0.0, 1.0, -2.0, 0.0, 2.0, -1.0, 0.0, 1.0 };
    vector<float> fDataX (fDataXArray, fDataXArray + sizeof(fDataXArray) / sizeof(float) );
    Filter sobelX(fDataX, 3, 3);

    // verify that your filter is correct by using it to filter an impulse image
    Image impulse = impulseImg(11); //create an image containing an impulse
    // convolve the impulse image with the Sobel kernel. We divide the output by 4 and
    // add 0.5 to make the range of the image between 0 and 1
    Image verifyKernel = sobelX.convolve(impulse)/4 + 0.5;
    verifyKernel.write("./Output/verifySobelKernel.png");

    // filter an image using the sobel kernel
    Image im2("./Input/lounge_view.png");
    Image sobelFiltered = sobelX.convolve(im2);

    // make the range of the output image from 0 to 1 for visualization
    // since the Sobel filter changes the range of a (0,1) image to (-2,2)
    Image sobelOut = sobelFiltered/4 + 0.5;
    sobelOut.write("./Output/sobelFiltered.png");
    // ---------------------------------------------------
    float sigma = 2.0f;
    
    Image lens_img("./Input/lens.png");

    Image sharpened_image = unsharpMask(im2, sigma);
    sharpened_image.write("./Output/sharpened_image.png");
        
    Image gaussianBlur_horizontal_img= gaussianBlur_horizontal(im2, sigma, 3.0, true);
    gaussianBlur_horizontal_img.write("./Output/gaussianBlur_horizontal_img.png");

    Image gaussianBlur_separable_img = gaussianBlur_separable(im2, sigma, 3.0, true);
    gaussianBlur_separable_img.write("./Output/gaussianBlur_separable_img.png");
    
    // --- Timer example ---------------------------------
    clock_t start = clock();
    Image gaussianBlur_2D_img = gaussianBlur_2D(im2, sigma);
    gaussianBlur_2D_img.write("./Output/gaussianBlur_2D.png");

    clock_t end = clock();
    double duration = (end-start)*1.0f/CLOCKS_PER_SEC;
    cout << "2D gaussian took: " << duration <<"s" << endl;
    // ---------------------------------------------------
    
    
    Image bilaterial_img = bilateral(lens_img);
    bilaterial_img.write("./Output/bilaterial_img.png");
    
    Image bilateral_yuv_img = bilaYUV(lens_img);
    bilateral_yuv_img.write("./Output/bilateral_yuv_img.png");
}

