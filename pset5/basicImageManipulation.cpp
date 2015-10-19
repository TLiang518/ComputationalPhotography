/* --------------------------------------------------------------------------
 * File:    basicImageManipulation.cpp
 * Created: 2015-09-23
 * --------------------------------------------------------------------------
 * 
 * 
 * 
 * ------------------------------------------------------------------------*/


#include "basicImageManipulation.h"
using namespace std;
  
  
// --------- HANDOUT PS03 ------------------------------
// -----------------------------------------------------
//
Image scaleNN(const Image &im, float factor){
    // --------- HANDOUT  PS03 ------------------------------
    // create a new image that is factor times bigger than the input by using
    // nearest neighbor interpolation.

    /*
    The simplest technique to sample the new value is called nearest-neighbor re-sampling: 
    we round-off the real coordinates to the nearest integers and use the input’s color at 
    this new location to be the value of the current output pixel.
    */
    int new_width = im.width()*factor;
    int new_height = im.height()*factor;
    Image output(new_width, new_height, im.channels());
    for (int a = 0; a < new_width; a++){
        for (int b = 0; b < new_height; b++){
            for (int c = 0; c < im.channels(); c++){
                output(a,b,c) = im(a/factor, b/factor, c);
            }
        }
    }
    return output;
}

/*
Nearest-neighbor re-sampling creates blocky artifacts and pixelated results. 
We will address this using a better reconstruction based on bilinear interpolation. 
For this, we consider the four pixels immediately around the computed real coordinates
 and perform two linear interpolations. We first linearly interpolate along x the colors
  of the top and bottom pairs of pixels. Then we interpolate these two values along y to 
  get the final sample. The interpolation weight are driven by the distance from the corners.
*/

float interpolateLin(const Image &im, float x, float y, int z, bool clamp){
    // --------- HANDOUT  PS03 ------------------------------
    // bilinear interpolation samples the value of a non-integral
    // position (x,y) from its four "on-grid" neighboring pixels.
    //  |           | 
    // -1-----------2- 
    //  |           |  *: my coordinates (x,y) are not integral
    //  |  *        |     since I am not on the pixel grid :(
    //  |           |  1: top-left
    //  |           |  2: top-right
    //  |           |  3: bottom-right
    // -4-----------3- 4: bottom-left, what are our coordinates?
    //  |           |    We are willing to share some color 
    //                   information with * ! Of course, the pixel
    //                   closest to * should influence it more.
    
    /*

    • Take 4 nearest neighbors
    • Weight according to x & y fractional
    coordinates
    • Can be done using two 1D linear
    reconstructions along x then y (or y then x)
    
    */
    std::cout << "Checking interpolate lin " << std::endl;
    int lower_x = floor(x);
    float lower_x_weight = abs(x - lower_x);
    int upper_x = ceil(x);
    float upper_x_weight = 1 - abs(x - upper_x);
    int lower_y = floor(y);
    float lower_y_weight = abs(y - lower_y);
    int upper_y = ceil(y);
    float upper_y_weight = 1 - abs(y - upper_y);

    cout << " lower_x " << lower_x << " upper_x " << upper_x << " lower_y " << lower_y << " upper_y " << upper_y << std::endl;
    cout << " lower_x_weight " << lower_x_weight << " upper_x_weight " << upper_x_weight << " lower_y_weight " << lower_y_weight << " upper_y_weight " << upper_y_weight << std::endl;


    float top_x_value = im(lower_x, upper_y, z)*lower_x_weight + im(upper_x, upper_y,z)*upper_x_weight;
    
    float bottom_x_value = im(lower_x, lower_y,z)*lower_x_weight + im(upper_x, lower_y,z)*upper_x_weight;

    cout << "top_x_value " << top_x_value << " bottom_x_value " << bottom_x_value << endl;

    float interpolated_y_value = top_x_value*upper_y_weight + bottom_x_value*lower_y_weight;

    return interpolated_y_value;
}

Image scaleLin(const Image &im, float factor){
    // --------- HANDOUT  PS03 ------------------------------
    // create a new image that is factor times bigger than the input by using
    // bilinear interpolation
    return im;
}


Image rotate(const Image &im, float theta) {
    // --------- HANDOUT  PS03 ------------------------------
    // (6.865 required, 6.815 extra credit)
    // rotate an image around its center by theta

	// center around which to rotate
    float centerX = (im.width()-1.0)/2.0;
    float centerY = (im.height()-1.0)/2.0;

    return im; // changeme
    
}

// -----------------------------------------------------
// --------- END --- PS03 ------------------------------





// --- NO need to edit below, this is the solution to PS01 ---





// --------- HANDOUT PS01 ------------------------------
// -----------------------------------------------------

// Change the brightness of the image
// const Image & means a reference to im will get passed to the function,
// but the compiler won't let you modify it within the function.
// So you will return a new image
Image brightness(const Image &im, float factor) {
    // // --------- HANDOUT  PS01 ------------------------------
	// // Image output(im.width(), im.height(), im.channels());
	// // Modify image brightness
	// // return output;
	// return Image(1,1,1); // Change this
    
    // --------- SOLUTION PS01 ------------------------------
    return im * factor;
}

Image contrast(const Image &im, float factor, float midpoint) {
    // // --------- HANDOUT  PS01 ------------------------------
    // // Image output(im.width(), im.height(), im.channels());
    // // Modify image contrast
    // // return output;
	// return Image(1,1,1); //Change this
    
    // --------- SOLUTION PS01 ------------------------------
    return (im - midpoint)*factor+midpoint;
}

Image color2gray(const Image &im, const std::vector<float> &weights) {
    // // --------- HANDOUT  PS01 ------------------------------
    // // Image output(im.width(), im.height(), 1);
    // // Convert to grayscale
	// return Image(1,1,1); //Change this

    // --------- SOLUTION PS01 ------------------------------
    Image output(im.width(), im.height(), 1);
    for (int i = 0 ; i < im.width(); i++ ) {
        for (int j = 0 ; j < im.height(); j++ ) {
            output(i,j,0) = im(i,j,0) * weights[0] + im(i,j,1) * weights[1] + im(i,j,2) *weights[2];
        }
    }
    return output;
}

// For this function, we want two outputs, a single channel luminance image 
// and a three channel chrominance image. Return them in a vector with luminance first
std::vector<Image> lumiChromi(const Image &im) {
    // // --------- HANDOUT  PS01 ------------------------------
    // // Create the luminance image
    // // Create the chrominance image
    // // Create the output vector as (luminance, chrominance)
	// return std::vector<Image>(); //Change this
    
    // --------- SOLUTION PS01 ------------------------------

    // Create the luminance
    Image im_luminance = color2gray(im);

    // Create chrominance images
    // We copy the input as starting point for the chrominance
    Image im_chrominance = im; 
    for (int c = 0 ; c < im.channels(); c++ )
    for (int y = 0 ; y < im.height(); y++) 
    for (int x = 0 ; x < im.width(); x++)
    {
        im_chrominance(x,y,c) = im_chrominance(x,y,c) / im_luminance(x,y);
    }

    // Stack luminance and chrominance in the output vector, luminance first
    std::vector<Image> output;
    output.push_back(im_luminance);
    output.push_back(im_chrominance);
    return output;
}


// Modify brightness then contrast
Image brightnessContrastLumi(const Image &im, float brightF, float contrastF, float midpoint) {
    // // --------- HANDOUT  PS01 ------------------------------
    // // Modify brightness, then contrast of luminance image
    // return Image(1,1,1); // Change this
    
    // --------- SOLUTION PS01 ------------------------------
    // Separate luminance and chrominance
    std::vector<Image> lumi_chromi = lumiChromi(im);
    Image im_luminance             = lumi_chromi[0];
    Image im_chrominance           = lumi_chromi[1];

    // Process the luminance channel
    im_luminance = brightness(im_luminance, brightF);
    im_luminance = contrast(im_luminance, contrastF, midpoint);

    // Multiply the chrominance with the new luminance to get the final image
    for (int i = 0 ; i < im.width(); i++ ){
        for (int j = 0 ; j < im.height(); j++) {
            for (int c = 0; c < im.channels(); c++) {
                im_chrominance(i,j,c) = im_chrominance(i,j,c) * im_luminance(i,j);
            }
        }
    }
    // At this point, im_chrominance olds the complete processed image
    return im_chrominance;
}


Image rgb2yuv(const Image &im) {
    // // --------- HANDOUT  PS01 ------------------------------
    // // Create output image of appropriate size
    // // Change colorspace
    // return Image(1,1,1); // Change this
    
    // --------- SOLUTION PS01 ------------------------------
    Image output(im.width(), im.height(), im.channels());
    for (int j = 0 ; j < im.height(); j++)
    for (int i = 0 ; i < im.width(); i++) 
    {
        output(i,j,0) =   0.299 * im(i,j,0) + 0.587 * im(i,j,1) + 0.114 * im(i,j,2);
        output(i,j,1) = - 0.147 * im(i,j,0) - 0.289 * im(i,j,1) + 0.436 * im(i,j,2);
        output(i,j,2) =   0.615 * im(i,j,0) - 0.515 * im(i,j,1) - 0.100 * im(i,j,2);
    }
    return output;
}


Image yuv2rgb(const Image &im) {
    // // --------- HANDOUT  PS01 ------------------------------
    // // Create output image of appropriate size
    // // Change colorspace
    // return Image(1,1,1); // Change this

    // --------- SOLUTION PS01 ------------------------------
    Image output(im.width(), im.height(), im.channels());
    for (int j = 0 ; j < im.height(); j++) 
    for (int i = 0; i < im.width(); i++) 
    {
        output(i,j,0) =  im(i,j,0) + 0     * im(i,j,1) + 1.14  * im(i,j,2);
        output(i,j,1) =  im(i,j,0) - 0.395 * im(i,j,1) - 0.581 * im(i,j,2);
        output(i,j,2) =  im(i,j,0) + 2.032 * im(i,j,1) + 0     * im(i,j,2);
    }
    return output;
}


Image saturate(const Image &im, float factor) {
    // // --------- HANDOUT  PS01 ------------------------------
    // // Create output image of appropriate size
    // // Saturate image
    // // return output; 
    // return Image(1,1,1); // Change this
    
    // --------- SOLUTION PS01 ------------------------------
    Image output = rgb2yuv(im); // Change colorspace
    for (int i = 0 ; i < im.width(); i++) {
        for (int j = 0 ; j < im.height(); j++) {
            output(i,j,1) = output(i,j,1) * factor;
            output(i,j,2) = output(i,j,2) * factor;
        }
    }
    output = yuv2rgb(output); // Back to RGB
    return output;  
}


// Return two images in a C++ vector
std::vector<Image> spanish(const Image &im) {
    // // --------- HANDOUT  PS01 ------------------------------
    // // Remember to create the output images and the output vector
    // // Push the images onto the vector
    // // Do all the required processing
    // // Return the vector, color image first
	// return std::vector<Image>(); //Change this
    
    // --------- SOLUTION PS01 ------------------------------
    // Extract the luminance
    Image output_L = color2gray(im);

    // Convert to YUV for manipulation
    Image output_C = rgb2yuv(im);

    for (int j = 0; j < im.height(); j++)
    for (int i = 0; i < im.width(); i++)
    {
        output_C(i,j,0) = 0.5; // constant luminance
        output_C(i,j,1) = -output_C(i,j,1); // opposite chrominance
        output_C(i,j,2) = -output_C(i,j,2); // opposite chrominance
    }
    // Convert back to RGB
    output_C = yuv2rgb(output_C);

    // Location of the black dot
    int bdot_x = floor(im.width()/2);
    int bdot_y = floor(im.height()/2);

    // Add the black dot to Luminance, and Chrominance images
    output_L(bdot_x, bdot_y,0) = 0.0f;
    output_C(bdot_x, bdot_y,0) = 0.0f; // black is 0
    output_C(bdot_x, bdot_y,1) = 0.0f;
    output_C(bdot_x, bdot_y,2) = 0.0f;

    // Pack the images in a vector, chrominance first
    std::vector<Image> output;
    output.push_back(output_C);
    output.push_back(output_L);
    return output; 
}


// White balances an image using the gray world assumption
Image grayworld(const Image & im) {
    // // --------- HANDOUT  PS01 ------------------------------
    // Implement automatic white balance by multiplying each channel
    // of the input by a factor such that the three channel of the output image
    // have the same mean value. The mean value of the green channel
    // is taken as reference.
    // return Image(1,1,1); // Change this
    
    // --------- SOLUTION PS01 ------------------------------
    // Compute the mean per channel
    float mean_r = 0, mean_g = 0, mean_b = 0;
    float N = im.width()*im.height();
    for (int j = 0 ; j < im.height(); j++)
    for (int i = 0 ; i < im.width(); i++) 
    {
        mean_r += im(i,j,0);
        mean_g += im(i,j,1);
        mean_b += im(i,j,2);
    }
    mean_r /= N;
    mean_g /= N;
    mean_b /= N;

    Image output = im;
    for (int j = 0 ; j < im.height();j ++)
    for (int i = 0 ; i < im.width(); i++)
    {
        output(i,j,0) = output(i,j,0)/mean_r*mean_g;
        // dont process output(i,j,1), since the mean of
        // the green channel is already at the right value
        output(i,j,2) = output(i,j,2)/mean_b*mean_g;
    }
    return output;
}

// -----------------------------------------------------
// --------- END --- PS01 ------------------------------
