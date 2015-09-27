/* -----------------------------------------------------------------
 * File:    a1.cpp
 * Created: 2015-09-15
 * -----------------------------------------------------------------
 * 
 * Assignment 01
 * 
 * ---------------------------------------------------------------*/


#include "a1.h"
using namespace std;

// Change the brightness of the image
// const Image & means a reference to im will get passed to the function,
// but the compiler won't let you modify it within the function.
// So you will return a new image
Image brightness(const Image &im, float factor) {
    // --------- HANDOUT  PS01 ------------------------------
    Image output_image(im.width(), im.height(), im.channels());
	for (int i = 0; i < im.width(); i++)
    {
        for (int j = 0; j < im.height(); j++)
        {
            for (int c = 0; c < im.channels(); c++)
            {
                output_image(i,j,c) = im(i,j,c)*factor;
            }
        }
    }
    return output_image;
}


Image contrast(const Image &im, float factor, float midpoint) {
    // --------- HANDOUT  PS01 ------------------------------
    // Image output(im.width(), im.height(), im.channels());
    // Modify image contrast
    // return output;
	Image output_image(im.width(), im.height(), im.channels());
    for (int i = 0; i < im.width(); i++)
    {
        for (int j = 0; j < im.height(); j++)
        {
            for (int c = 0; c < im.channels(); c++)
            {
                //Handles inputs less than 0 and greater than 1 in a sensible way by
                //Taking the maximum of 0 and the contrasted pixel value (handles less than 0)
                //and then the minimum of that and 1 (handling greater than 1 values)
                output_image(i,j,c) = std::min(float(1), std::max(float(0), factor*(im(i,j,c) - midpoint) + midpoint));
            }
        }
    }
    return output_image;
}


Image color2gray(const Image &im, const std::vector<float> &weights) {
    // --------- HANDOUT  PS01 ------------------------------
    
    Image output(im.width(), im.height(), 1);
    
    for (int i = 0; i < im.width(); i++)
    {
        for (int j = 0; j < im.height(); j++)
        {
            output(i,j) = im(i,j,0)*weights.at(0) + im(i,j,1)*weights.at(1) + im(i,j,2)*weights.at(2);
        }
    }

    return output;	
}


// For this function, we want two outputs, a single channel luminance image 
// and a three channel chrominance image. Return them in a vector with luminance first
std::vector<Image> lumiChromi(const Image &im) {
    // --------- HANDOUT  PS01 ------------------------------
    // Create the luminance image
    // Create the chrominance image
    // Create the output vector as (luminance, chrominance)
	
    std::vector<Image> lumiChromiVect;

    //Default weights
    Image lumiImage = color2gray(im);

    lumiChromiVect.push_back(lumiImage);
    
    Image chromiImage = im;
    for (int i = 0; i < im.width(); i++)
    {
        for (int j = 0; j < im.height(); j++)
        {
            float lumiDivisor = lumiImage(i,j);

            for (int c = 0; c < im.channels(); c++)
            {
                if (lumiDivisor == 0)
                {
                    chromiImage(i,j,c) = 0;
                }
                else
                {
                    chromiImage(i,j,c) = chromiImage(i,j,c)/lumiDivisor;
                }
            }
        }
    }
    lumiChromiVect.push_back(chromiImage);
    return lumiChromiVect;
}

// Modify brightness then contrast
Image brightnessContrastLumi(const Image &im, float brightF, float contrastF, float midpoint) {
    // --------- HANDOUT  PS01 ------------------------------
    // Modify brightness, then contrast of luminance image
    std::vector<Image> lumiChromiVect = lumiChromi(im);
    Image lumiImage = lumiChromiVect.at(0);
    Image chromiImage = lumiChromiVect.at(1);
    
    Image blumiImage = brightness(lumiImage, brightF);
    Image bcLumiImage = contrast(blumiImage,contrastF, midpoint);
    
    Image recombinedImage = chromiImage;
    for (int i = 0; i < im.width(); i++){
        for (int j = 0; j < im.height(); j++){
            for (int c = 0; c < im.channels(); c++)
            {
                recombinedImage(i,j,c) = chromiImage(i,j,c)*bcLumiImage(i,j);
            }
        }
    }
    return recombinedImage;
}


Image rgb2yuv(const Image &im) {
    // --------- HANDOUT  PS01 ------------------------------
    // Create output image of appropriate size
    // Change colorspace
    
    Image yuv = im;

    static const float row_1[] = {0.299, 0.587, 0.114};
    std::vector<float> row_1_vect (row_1, row_1 + sizeof(row_1) / sizeof(row_1[0]) );

    static const float row_2[] = {-0.147, -0.289, 0.436};
    std::vector<float> row_2_vect (row_2, row_2 + sizeof(row_2) / sizeof(row_2[0]) );

    static const float row_3[] = {0.615, -0.515, -0.100};
    std::vector<float> row_3_vect (row_3, row_3 + sizeof(row_3) / sizeof(row_3[0]) );

    std::vector<std::vector<float> > rgb2yuv_factors;
    rgb2yuv_factors.push_back(row_1_vect);
    rgb2yuv_factors.push_back(row_2_vect);
    rgb2yuv_factors.push_back(row_3_vect);

    for (int i = 0; i < im.width(); i++){
            for (int j = 0; j < im.height(); j++){                
                for (int c = 0; c < im.channels(); c++){
                    yuv(i,j,c) = rgb2yuv_factors.at(c).at(0)*im(i,j,0) 
                    + rgb2yuv_factors.at(c).at(1)*im(i,j,1) 
                    + rgb2yuv_factors.at(c).at(2)*im(i,j,2); 
            }
        }
    }
    return yuv;
}


Image yuv2rgb(const Image &im) {
    // --------- HANDOUT  PS01 ------------------------------
    // Create output image of appropriate size
    // Change colorspace
    Image rgb = im;

    static const float row_1[] = {1, 0, 1.14};
    std::vector<float> row_1_vect (row_1, row_1 + sizeof(row_1) / sizeof(row_1[0]) );

    static const float row_2[] = {1, -0.395, -00.581};
    std::vector<float> row_2_vect (row_2, row_2 + sizeof(row_2) / sizeof(row_2[0]) );

    static const float row_3[] = {1, 2.032, 0};
    std::vector<float> row_3_vect (row_3, row_3 + sizeof(row_3) / sizeof(row_3[0]) );

    std::vector<std::vector<float> > yuv2rgb_factors;
    yuv2rgb_factors.push_back(row_1_vect);
    yuv2rgb_factors.push_back(row_2_vect);
    yuv2rgb_factors.push_back(row_3_vect);

    for (int i = 0; i < im.width(); i++){
            for (int j = 0; j < im.height(); j++){                
                for (int c = 0; c < im.channels(); c++){
                    rgb(i,j,c) = yuv2rgb_factors.at(c).at(0)*im(i,j,0)
                     + yuv2rgb_factors.at(c).at(1)*im(i,j,1) 
                     + yuv2rgb_factors.at(c).at(2)*im(i,j,2); 
            }
        }
    }
    return rgb;
}


Image saturate(const Image &im, float factor) {
    // --------- HANDOUT  PS01 ------------------------------
    // Create output image of appropriate size
    // Saturate image
    // return output; 

    Image yuv = rgb2yuv(im);
    Image saturated_yuv = yuv;
    for (int i = 0; i < im.width(); i++){
            for (int j = 0; j < im.height(); j++){                
                for (int c = 1; c < im.channels(); c++){
                    saturated_yuv(i,j,c) = yuv(i,j,c) * factor;
            }
        }
    }

    Image saturated_rgb = yuv2rgb(saturated_yuv);
    return saturated_rgb;
}


// Return two images in a C++ vector
std::vector<Image> spanish(const Image &im) {
    // --------- HANDOUT  PS01 ------------------------------
    // Remember to create the output images and the output vector
    // Push the images onto the vector
    // Do all the required processing
    // Return the vector, color image first

    Image rgb_uv_negative = saturate(im, -1.0);
    Image yuv_uv_negative = rgb2yuv(rgb_uv_negative);
    Image saturated_yuv_constant_luminance = yuv_uv_negative;
    for (int i = 0; i < yuv_uv_negative.width(); i++){
        for (int j = 0; j < yuv_uv_negative.height(); j++){                
            saturated_yuv_constant_luminance(i,j,0) = 0.5;
        }
    }
    Image color_image = yuv2rgb(saturated_yuv_constant_luminance);

    Image gray = color2gray(im);
    
    Image gray_with_dot = add_black_dot(gray);
    Image color_with_dot = add_black_dot(color_image);

    std::vector<Image> spanish_vector;
    spanish_vector.push_back(color_with_dot);
    spanish_vector.push_back(gray_with_dot);

    return spanish_vector;
}

Image add_black_dot(const Image &im)
{
    Image black_dot_image = im;
    for (int c = 0; c < im.channels(); c++)
        {
        black_dot_image(std::floor(black_dot_image.width()/2.0), std::floor(black_dot_image.height()/2.0), c) = 0;
    }
    return black_dot_image;
}


// White balances an image using the gray world assumption
Image grayworld(const Image & im) {
    // --------- HANDOUT  PS01 ------------------------------
    // Implement automatic white balance by multiplying each channel
    // of the input by a factor such that the three channel of the output image
    // have the same mean value. The mean value of the green channel
    // is taken as reference.
    Image whitebalanced = im;

    std::vector<float> mean_per_channel = get_mean_per_channel(im);
    float red_value = mean_per_channel.at(0);
    float green_value = mean_per_channel.at(1);
    float blue_value = mean_per_channel.at(2);

    std::vector<float> multiplier_per_channel;
    float red_multiplier = green_value/red_value;
    float green_multiplier = green_value/green_value;
    float blue_multiplier = green_value/blue_value;
    multiplier_per_channel.push_back(red_multiplier);
    multiplier_per_channel.push_back(green_multiplier);
    multiplier_per_channel.push_back(blue_multiplier);

    std::cout << "red_multiplier " << red_multiplier << std::endl;
    std::cout << "green_multiplier " << green_multiplier << std::endl;
    std::cout << "green_multiplier " << blue_multiplier << std::endl;

    for (int i = 0; i < im.width(); i++){
        for (int j = 0; j < im.height(); j++){                
            for (int c = 0;c < im.channels(); c++){
                whitebalanced(i,j,c) = im(i,j,c)*multiplier_per_channel.at(c);        
            }
        }
    }    
    return whitebalanced;
}

std::vector<float> get_mean_per_channel(const Image & im)
{
    
    std::vector<float> channel_mean;
    channel_mean.push_back(float(0));
    channel_mean.push_back(float(0));
    channel_mean.push_back(float(0));

    int image_dimensions = im.width() * im.height();

    for (int i = 0; i < im.width(); i++){
        for (int j = 0; j < im.height(); j++){ 
            for (int c = 0; c < im.channels(); c++){
                channel_mean.at(c) = channel_mean.at(c) + im(i,j,c)/image_dimensions;
            } 
        }
    }

    return channel_mean;
}