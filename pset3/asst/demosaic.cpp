/* --------------------------------------------------------------------------
 * File:    demosaic.cpp
 * Created: 2015-10-01
 * --------------------------------------------------------------------------
 * 
 * 
 * 
 * ------------------------------------------------------------------------*/


#include "demosaic.h"
#include <cmath>

using namespace std;


Image basicGreen(const Image &raw, int offset){
    // --------- HANDOUT  PS03 ------------------------------
    // Takes as input a raw image and returns a single-channel
    // 2D image corresponding to the green channel using simple interpolation

    std::cout << "raw.width() " << raw.width() << " raw.height() " << raw.height() << " raw.channels() " << raw.channels() << std::endl;
    Image output(raw.width(), raw.height(), 1);

    bool condition;
    if (offset % 2 == 0){
        condition = true;
    }
    else {
        condition = false;
    }

    for (int a = 1; a < raw.width() - 1; a++){
        for (int b = 1; b < raw.height() - 1; b++){
            if ((a % 2 == b % 2) == condition){
                output(a,b) = raw(a,b);
            }
            else{
                output(a,b) = (raw(a-1,b) + raw(a+1,b) + raw(a,b+1) + raw(a,b-1))/4.0; 
            }
        }
    }
    return output;
}

Image basicRorB(const Image &raw, int offsetX, int offsetY){
    // --------- HANDOUT  PS03 ------------------------------
    //  Takes as input a raw image and returns a single-channel
    //  2D image corresponding to the red or blue channel using simple interpolation
    
    /*
    
    Similar to the green-channel case, copy the values when they are available. 
    For interpolated pixels that have two direct neighbors that are known (left-right or up-down), 
    simply take the linear interpolation between the two values. For the remaining case, interpolate 
    the four diagonal pixels.

    */

    Image output(raw.width(), raw.height(), 1);
    for (int a = 1; a < raw.width() - 1; a++){
        for (int b = 1; b < raw.height() - 1; b++){
            if (a % 2 == offsetX && b % 2 == offsetY) {
                output(a,b) = raw(a,b);
            }
            else if (a % 2 == offsetX && b % 2 != offsetY){
                output(a,b) = (raw(a,b-1) + raw(a,b+1))/2.0;
            }
            else if (a % 2 != offsetX && b % 2 == offsetY){
                output(a,b) = (raw(a-1,b) + raw(a+1,b))/2.0;
            }
            else{
                output(a,b) = (raw(a-1,b-1) + raw(a-1,b+1) + raw(a+1,b-1) + raw(a+1,b+1))/4.0; 
            }
        }
    }
    return output;
}

Image basicDemosaic(const Image &raw, int offsetGreen, int offsetRedX, int offsetRedY, int offsetBlueX, int offsetBlueY){
    // --------- HANDOUT  PS03 ------------------------------
    // takes as input a raw image and returns an rgb image
    // using simple interpolation to demosaic each of the channels
    Image output(raw.width(), raw.height(), 3);
    Image green = basicGreen(raw, offsetGreen);
    Image red = basicRorB(raw, offsetRedX, offsetRedY);
    Image blue = basicRorB(raw, offsetBlueX, offsetBlueY);
    for (int a = 0; a < raw.width(); a++){
        for (int b = 0; b < raw.height(); b++){
            output(a,b,0) = red(a,b);
            output(a,b,1) = green(a,b);
            output(a,b,2) = blue(a,b);
        }
    }
    return output;
}

Image edgeBasedGreen(const Image &raw, int offset){
    // --------- HANDOUT  PS03 ------------------------------
    // Takes a raw image and outputs a single-channel
    // image corresponding to the green channel taking into account edges
    
    Image output(raw.width(), raw.height(), 1);

    bool condition;
    if (offset % 2 == 0){
        condition = true;
    }
    else {
        condition = false;
    }

    for (int a = 1; a < raw.width() - 1; a++){
        for (int b = 1; b < raw.height() - 1; b++){
            if ((a % 2 == b % 2) == condition){
                output(a,b) = raw(a,b);
            }
            else{
                float row_diff = abs(raw(a-1,b) - raw(a+1,b));
                float col_diff = abs(raw(a,b-1) - raw(a,b+1));
                if (col_diff < row_diff){
                    output(a,b) = (raw(a,b-1) + raw(a,b+1))/2.0;
                }
                else {
                    output(a,b) = (raw(a-1, b) + raw(a+1, b))/2.0;

                }
            }
        }
    }
    return output;
}

Image edgeBasedGreenDemosaic(const Image &raw, int offsetGreen, int offsetRedX, int offsetRedY, int offsetBlueX, int offsetBlueY){
    // --------- HANDOUT  PS03 ------------------------------
    // Takes as input a raw image and returns an rgb image
    // using edge-based green demosaicing for the green channel and
    // simple interpolation to demosaic the red and blue channels
    std::cout << "in edge based green demosaic" << std::endl;
    Image output(raw.width(), raw.height(), 3);
    Image green = edgeBasedGreen(raw, offsetGreen);
    Image red = basicRorB(raw, offsetRedX, offsetRedY);
    Image blue = basicRorB(raw, offsetBlueX, offsetBlueY);
    for (int a = 0; a < raw.width(); a++){
        for (int b = 0; b < raw.height(); b++){
            output(a,b,0) = red(a,b);
            output(a,b,1) = green(a,b);
            output(a,b,2) = blue(a,b);
        }
    }
    return output;
}


Image greenBasedRorB(const Image &raw, Image &green, int offsetX, int offsetY){
    // --------- HANDOUT  PS03 ------------------------------
    // Takes as input a raw image and returns a single-channel
    // 2D image corresponding to the red or blue channel using green based interpolation
    std::cout << "in green based r or b" << std::endl;
    Image three_channel_green = raw;
    for (int a = 0; a < raw.width(); a ++) {
        for (int b = 0; b < raw.height(); b ++) {
            for (int c = 0; c < raw.channels(); c++){
                three_channel_green(a,b,c) = green(a,b);
            }
        }
    }
    Image raw_minus_green = raw - three_channel_green;
    Image basic_r_or_b_minus_green = basicRorB(raw_minus_green, offsetX, offsetY);
    Image greenBasedRorB_im = basic_r_or_b_minus_green + green;
    return greenBasedRorB_im;
}

Image improvedDemosaic(const Image &raw, int offsetGreen, int offsetRedX, int offsetRedY, int offsetBlueX, int offsetBlueY){
    // --------- HANDOUT  PS03 ------------------------------
    // Takes as input a raw image and returns an rgb image
    // using edge-based green demosaicing for the green channel and
    // simple green based demosaicing of the red and blue channels
    std::cout << "in edge based green demosaic" << std::endl;
    Image output(raw.width(), raw.height(), 3);
    Image green = edgeBasedGreen(raw, offsetGreen);
    Image red = greenBasedRorB(raw, green, offsetRedX, offsetRedY);
    Image blue = greenBasedRorB(raw, green, offsetBlueX, offsetBlueY);
    for (int a = 0; a < raw.width(); a++){
        for (int b = 0; b < raw.height(); b++){
            output(a,b,0) = red(a,b);
            output(a,b,1) = green(a,b);
            output(a,b,2) = blue(a,b);
        }
    }
    return output;
    
}
