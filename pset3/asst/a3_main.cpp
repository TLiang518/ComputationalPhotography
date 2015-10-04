/* -----------------------------------------------------------------
 * File:    a3_main.cpp
 * Author:  Michael Gharbi <gharbi@mit.edu>
 * Created: 2015-09-30
 * -----------------------------------------------------------------
 * 
 * 
 * 
 * ---------------------------------------------------------------*/


#include "Image.h"
#include "basicImageManipulation.h"
#include "demosaic.h"
#include "align.h"
#include <iostream>
#include <sstream>
#include <iomanip>
#include <vector>


using namespace std;


// This is a way for you to test your functions. 
// We will only grade the contents of demosaic.cpp and align.cpp
int main() {
    // // Denoise ---------------------------
    // // Load sequence
    
    std::cout << "Starting" << std::endl;
    /*
    vector<Image> seq;
    int n_images = 20;
    for (int i = 1; i <= n_images; ++i) {
        ostringstream fname;
        // fname << "./Input/aligned-ISO400/1D2N-iso400-under-";
        fname << "./Input/aligned-ISO3200/1D2N-iso3200-";
        fname << i;
        fname << ".png";
        seq.push_back(Image(fname.str()));
    }
    
    /*
    
    // Denoise
    Image out = denoiseSeq(seq);
    out.write("./Output/denoised.png");

    */

    /*
    std::cout << "Starting logSNR" << std::endl;
    
    Image SNRIm = logSNR(seq,float(1/30.0));
    SNRIm.write("./Output/snr_map.png");

    */

    Image jacqui("./Input/jacqui_larg.png");
    Image jacqui_offset_14_10("./Input/jacqui_larg_14_10.png");
    Image samoyed("./Input/samoyed.png");
    Image samoyed_moved("./Input/samoyed_moved.png");
    Image doge_small("./Input/doge_small.png");
    Image doge_2("./Input/doge_2.png");
    Image doge_3("./Input/doge_3.png");
    Image doge_4("./Input/doge_4.png");



    /*
    vector<int> image_offset = align(jacqui,jacqui,20);

    for (auto & element: image_offset){
        std::cout << element << " ";
    }
    std::cout << "Should have been 0,0" << std::endl;
    */

    /*
    vector<int> image_offset_2 = align(jacqui,jacqui_offset_14_10,20);
    for (auto & element: image_offset_2){
        std::cout << element << " ";
    }
    std::cout << "Should have been 14,10" << std::endl;
    */

    /*
    vector<int> image_offset_2 = align(doge_small, doge_2,20);
    for (auto & element: image_offset_2){
        std::cout << element << " ";
    }
    std::cout << std::endl;

    vector<int> image_offset_3 = align(doge_small,doge_3,20);
    for (auto & element: image_offset_3){
        std::cout << element << " ";
    }
    std::cout << std::endl;


    vector<int> image_offset_4 = align(doge_small, doge_4,20);
    for (auto & element: image_offset_4){
        std::cout << element << " ";
    }
    std::cout << std::endl;
    */
    
    // Demosaic ---------------------------
    /*
    Image raw("./Input/raw/signs-small.png");
    Image green = basicGreen(raw, 1);
    green.write("./Output/demosaic_green.png");
    Image red = basicRorB(raw, 1, 1);
    red.write("./Output/demosaic_red.png");
    Image blue = basicRorB(raw, 0, 0);
    blue.write("./Output/demosaic_blue.png");
    Image rgb = basicDemosaic(raw, 1, 1,1,0,0);
    rgb.write("./Output/demosaiced.png");

    Image edge_green = edgeBasedGreen(raw, 1);
    edge_green.write("./Output/demosaic_edge_green.png");

    Image rgb_edge = edgeBasedGreenDemosaic(raw, 1, 1,1,0,0);
    rgb_edge.write("./Output/demosaiced_edge.png");


    Image red_greenBasedRorB = greenBasedRorB(raw, edge_green, 1, 1);
    red_greenBasedRorB.write("./Output/demosaic_greenBasedRorB_red.png");
    Image blue_greenBasedRorB = greenBasedRorB(raw, edge_green, 0, 0);
    blue_greenBasedRorB.write("./Output/demosaic_greenBasedRorB_blue.png");
    Image rgb_greenBasedRorB = improvedDemosaic(raw, 1, 1,1,0,0);
    rgb_greenBasedRorB.write("./Output/demosaiced_greenBasedRorB.png");
    */

    //
    // // Sergey ---------------------------
    Image sergeyImg("./Input/Sergey/00088v_third.png");
    Image rgb2 = split(sergeyImg);
    rgb2.write("./Output/Sergey_split.png");
    Image rgbAlign = sergeyRGB(sergeyImg,10);
    rgbAlign.write("./Output/Sergey_aligned.png");

    return 0;
}

