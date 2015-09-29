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
                for (int k_i_counter = 0; k_i_counter < k; k_i_counter++){
                    for (int k_j_counter = 0; k_j_counter < k; k_j_counter++){
                        avg += im.smartAccessor(i+k_i_counter, j+k_j_counter, c, clamp)/(k*k);
                    }
                }
                output_image(i,j,c) = avg;
            }
        }
    }
    return output_image;
}

Image Filter::convolve(const Image &im, bool clamp){
    // --------- HANDOUT  PS02 ------------------------------
    // Write a convolution function for the filter class
    //std::cout << "kernel " << kernel << std::endl;
    Image output_image(im.width(), im.height(), im.channels());
    for (int i = 0; i < im.width(); i ++){
        for (int j = 0; j < im.height(); j++){
            for (int c = 0; c < im.channels(); c++){                
                float conv_value = 0; 
                for (int kernel_a = -width/2; kernel_a <= width/2; kernel_a++){
                    for (int kernel_b = -height/2; kernel_b <= height/2; kernel_b++){
                        conv_value += im.smartAccessor(i + kernel_a, j + kernel_b, c, clamp) * operator()(width/2 - kernel_a, height/2 - kernel_b);
                    }
                }
                output_image(i,j,c) = conv_value;
            }
        }
    }
    return output_image;
}



Image boxBlur_filterClass(const Image &im, int k, bool clamp) {
    // --------- HANDOUT  PS02 ------------------------------
    // Reimplement the box filter using the filter class.
    // check that your results match those in the previous function "boxBlur"
    Filter boxBlur(k, k);
    for (int i = 0; i < k; i ++){
        for (int j = 0; j < k; j++){
            boxBlur(i,j) = float(1.0/(k*k));     
        }
    }

    Image blurred_image = boxBlur.convolve(im, clamp);
    return blurred_image;
}


Image gradientMagnitude(const Image &im, bool clamp){
    // --------- HANDOUT  PS02 ------------------------------
    // Uses a Sobel kernel to compute the horizontal and vertical
    // components of the gradient of an image and returns the gradient magnitude.

    static const float horizontal_component_array[] = {-1, 0, 1, -2, 0, 2, -1, 0, 1};
    vector<float> horizontal_component_vec(horizontal_component_array, horizontal_component_array + sizeof(horizontal_component_array) / sizeof(horizontal_component_array[0]) );
    Filter horizontal_component(horizontal_component_vec, 3, 3);
    //−1 0 1 −1 −2 −1 −202and0 0 0 −1 0 1 1 2 1

    static const float flatgauss2DFilterArray[] = {-1, -2, -1, 0, 0, 0, 1, 2, 1};
    vector<float> vertical_component_vec(flatgauss2DFilterArray, flatgauss2DFilterArray + sizeof(flatgauss2DFilterArray) / sizeof(flatgauss2DFilterArray[0]) );
    Filter vertical_component(vertical_component_vec, 3, 3);

    Image horizontal_component_img = horizontal_component.convolve(im, clamp);
    Image vertical_component_img = vertical_component.convolve(im, clamp);

    Image final_image = horizontal_component_img;

    for (int i = 0; i < final_image.width(); i++){
        for (int j = 0; j < final_image.height(); j++){
            for (int c = 0; c < final_image.channels(); c++){
                final_image(i,j,c) = sqrt(pow(horizontal_component_img(i,j,c), 2.0) + pow(vertical_component_img(i,j,c), 2.0));
            }
        }
    }

    return final_image;
}

vector<float> gauss1DFilterValues(float sigma, float truncate){
    // --------- HANDOUT  PS02 ------------------------------
    // Create a vector containing the normalized values in a 1D Gaussian filter
    // Truncate the gaussian at truncate*sigma.
    // 1+2*ceil(sigma * truncate)
    
    vector<float> gauss1DFilterValues; 
    
    float continual_sum = 0;
    for (int x = (-ceil(sigma * truncate)); x < (ceil(sigma * truncate + 1)); x++){
        float new_val = exp(-1*pow(x,2)/(2*pow(sigma,2)));
        continual_sum = continual_sum + new_val;
        gauss1DFilterValues.push_back(new_val);
    }

    for (auto & element : gauss1DFilterValues) {
        element = element/continual_sum;
    }

    return gauss1DFilterValues;
}

Image gaussianBlur_horizontal(const Image &im, float sigma, float truncate, bool clamp){
    // --------- HANDOUT  PS02 ------------------------------
    // Gaussian blur across the rows of an image
    vector<float> gaussian_kernel = gauss1DFilterValues(sigma, truncate);
    Filter gaussian_filter(gaussian_kernel, 1+2*ceil(sigma * truncate), 1);
    Image horizontal_gaussian_blur =  gaussian_filter.convolve(im, clamp);
    return horizontal_gaussian_blur;
}

Image gaussianBlur_vertical(const Image &im, float sigma, float truncate, bool clamp){
    // --------- HANDOUT  PS02 ------------------------------
    // Gaussian blur across the rows of an image
    vector<float> gaussian_kernel = gauss1DFilterValues(sigma, truncate);
    Filter gaussian_filter(gaussian_kernel, 1, 1+2*ceil(sigma * truncate));
    Image vertical_gaussian_blur =  gaussian_filter.convolve(im, clamp);
    return vertical_gaussian_blur;
}

vector<float> gauss2DFilterValues(float sigma, float truncate){
    // --------- HANDOUT  PS02 ------------------------------
    // create a vector containing the normalized values in a 2D Gaussian
    // filter. Truncate the gaussian at truncate*sigma.
    vector<float> gauss2DFilterValuesVect;
    float continual_sum = 0;
    for (int a = -ceil(sigma * truncate); a < ceil(sigma * truncate) + 1; a++){
        for (int b = -ceil(sigma * truncate); b < ceil(sigma * truncate) + 1; b++){
            float x = sqrt(pow(a,2) + pow(b,2));
            float new_val = exp(-1*pow(x,2)/(2*pow(sigma,2)));
            continual_sum = continual_sum + new_val;
            gauss2DFilterValuesVect.push_back(new_val);
        }
    }

    for (auto & element : gauss2DFilterValuesVect) {
        element = element/continual_sum;
    }

    return gauss2DFilterValuesVect;
}


Image gaussianBlur_2D(const Image &im, float sigma, float truncate, bool clamp){
    // --------- HANDOUT  PS02 ------------------------------
    //  Blur an image with a full  full 2D rotationally symmetric Gaussian kernel
    vector<float> gauss2DFilterVector = gauss2DFilterValues(sigma, truncate);
    Filter gauss2DFilter(gauss2DFilterVector, 1+2*ceil(sigma * truncate), 1+2*ceil(sigma * truncate));
    Image gaussianBlur_2D_img = gauss2DFilter.convolve(im, clamp);
    return gaussianBlur_2D_img;
}

Image gaussianBlur_separable(const Image &im, float sigma, float truncate, bool clamp){
    // --------- HANDOUT  PS02 ------------------------------
    // Use principles of seperabiltity to blur an image using 2 1D Gaussian Filters
    Image horizontal_gaussian_blurred_im = gaussianBlur_horizontal(im, sigma, truncate, clamp);
    Image final_gaussian_blurred_im = gaussianBlur_vertical(horizontal_gaussian_blurred_im, sigma, truncate, clamp);
    return final_gaussian_blurred_im;
}


Image unsharpMask(const Image &im, float sigma, float truncate, float strength, bool clamp){
    // --------- HANDOUT  PS02 ------------------------------
    // sharpen an image
    Image lowpassed_extraction = gaussianBlur_separable(im, sigma, truncate, clamp);
    Image highpassed = im - lowpassed_extraction;
    Image final_image = im + highpassed*strength;
    return final_image;
}


Image bilateral(const Image &im, float sigmaRange, float sigmaDomain, float truncateDomain, bool clamp){
    // --------- HANDOUT  PS02 ------------------------------
    // Denoise an image using the bilateral filter
    Image output_image = im;
    float kernel_dimenions = 2*sigmaDomain * truncateDomain + 1;
    vector<float> gaussian_2d_domain = gauss2DFilterValues(sigmaDomain, truncateDomain);   
    Filter gaussian_2d_domain_filter(gaussian_2d_domain, kernel_dimenions, kernel_dimenions);

    float gaussian_range_constant = 1.0/sqrt(2*M_PI*pow(sigmaRange,2));
    // Read http://people.csail.mit.edu/sparis/publi/2009/fntcgv/Paris_09_Bilateral_filtering.pdf for clarification
    for (int i = 0; i < im.width(); i++){
        for (int j = 0; j < im.height(); j++){
            for (int c = 0; c < im.channels(); c++){
                float filter_value = 0; 
                float k_for_i_j = 0;
                
                for (int kernel_a = -kernel_dimenions/2; kernel_a <= kernel_dimenions/2; kernel_a++){
                    for (int kernel_b = -kernel_dimenions/2; kernel_b <= kernel_dimenions/2; kernel_b++){
                        
                        float spacial_weight = gaussian_2d_domain_filter(kernel_dimenions/2 - kernel_a, kernel_dimenions/2 - kernel_b);
                        
                        float x = im.smartAccessor(i,j,c);
                        float x_prime = im.smartAccessor(i + kernel_a, j + kernel_b, c, clamp);
                        float range_weight = gaussian_range_constant*exp(-1*pow((x - x_prime),2)/(2*pow(sigmaRange,2)));
                        
                        float intensity_at_examined_point = x_prime;

                        float final_value = spacial_weight*range_weight*intensity_at_examined_point;
                        filter_value += final_value;

                        k_for_i_j += spacial_weight*range_weight;
                    
                    }   
                }
                
                output_image(i,j,c) = filter_value/k_for_i_j;
            }
        }
    }
    return output_image;
}


Image bilaYUV(const Image &im, float sigmaRange, float sigmaY, float sigmaUV, float truncateDomain, bool clamp){
    // --------- HANDOUT  PS02 ------------------------------
    // 6.865 only
    // Bilaterial Filter an image seperatly for
    // the Y and UV components of an image

    Image yuv_im = rgb2yuv(im); 
    Image output_image = yuv_im; 

    float kernel_Y_dimenions = 2*sigmaY*truncateDomain + 1;
    vector<float> gaussian_2d_domain_Y = gauss2DFilterValues(sigmaY, truncateDomain);   
    Filter gaussian_2d_domain_filter_Y(gaussian_2d_domain_Y, kernel_Y_dimenions, kernel_Y_dimenions);

    float kernel_UV_dimenions = 2*sigmaUV*truncateDomain + 1;
    vector<float> gaussian_2d_domain_UV = gauss2DFilterValues(sigmaUV, truncateDomain);   
    Filter gaussian_2d_domain_filter_UV(gaussian_2d_domain_UV, kernel_UV_dimenions, kernel_UV_dimenions);

    float gaussian_range_constant = 1.0/sqrt(2*M_PI*pow(sigmaRange,2));
    // Read http://people.csail.mit.edu/sparis/publi/2009/fntcgv/Paris_09_Bilateral_filtering.pdf for clarification
    for (int i = 0; i < yuv_im.width(); i++){
        for (int j = 0; j < yuv_im.height(); j++){
            for (int c = 0; c < yuv_im.channels(); c++){
                float filter_value = 0; 
                float k_for_i_j = 0;
                
                float kernel_dimenions;
                Filter gaussian_2d_domain_filter = gaussian_2d_domain_filter_Y;
                if (c == 0){
                    kernel_dimenions = kernel_Y_dimenions;
                    gaussian_2d_domain_filter = gaussian_2d_domain_filter_Y;
                }
                else{
                    kernel_dimenions = kernel_UV_dimenions;
                    gaussian_2d_domain_filter = gaussian_2d_domain_filter_UV;                    
                }

                for (int kernel_a = -kernel_dimenions/2; kernel_a <= kernel_dimenions/2; kernel_a++){
                    for (int kernel_b = -kernel_dimenions/2; kernel_b <= kernel_dimenions/2; kernel_b++){

                        float spacial_weight = gaussian_2d_domain_filter(kernel_dimenions/2 - kernel_a, kernel_dimenions/2 - kernel_b);
                        
                        float x = yuv_im.smartAccessor(i,j,c);
                        float x_prime = yuv_im.smartAccessor(i + kernel_a, j + kernel_b, c, clamp);
                        float range_weight = gaussian_range_constant*exp(-1*pow((x - x_prime),2)/(2*pow(sigmaRange,2)));
                        
                        float intensity_at_examined_point = x_prime;

                        float final_value = spacial_weight*range_weight*intensity_at_examined_point;
                        filter_value += final_value;

                        k_for_i_j += spacial_weight*range_weight;
                    
                    }   
                }
                
                output_image(i,j,c) = filter_value/k_for_i_j;
            }
        }
    }

    std::cout << "made it to bottom!" << std::endl;

    Image rgb_output_image = yuv2rgb(output_image);
    return rgb_output_image;
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
