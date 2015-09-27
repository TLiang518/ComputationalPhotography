#include "a1.h"
#include <iostream>

using namespace std;


// This is a way for you to test your functions. 
// We will only grade the contents of a1.cpp and Image.cpp
int main() {
    cout << "nothing done in a1_main.cpp, debug me !" << endl;

    // // Example
    //Testing Image Functions
    Image black_square("./Input/square-512.png");
    //Size of this image is 3x512x512 -- should be equal to
    //786432
    long num_elems = black_square.number_of_elements();
    std::cout << "Is num_elems as expected? " << (num_elems == 786432) << std::endl;

    long x = black_square(300);
    std::cout << "What is at 300 for black square? " << x << std::endl;

    Image white_square("./Input/white_square.png");
    long x_white = white_square(300);
    std::cout << "What is at 300 for white square? " << x_white << std::endl;

    Image green_square("./Input/green_square.png");
    long x_green = green_square(300);
    std::cout << "What is at 300 for green square? " << x_green << std::endl;

    std::cout << "green square width and height: " << green_square.width() << ", " << green_square.height() << std::endl;
    std::cout << "Testing x y z" << std::endl << green_square(100,100,2) << std::endl;

    Image im("./Input/castle_small.png");

    Image brightness_tester = brightness(im,2.0);
    brightness_tester.write("./Output/brightness_castle_small.png");

    Image brighter_green_square = contrast(im,2.0,1.0);
    brighter_green_square.write("./Output/contrast_castle_small.png");

    static const float weights_array[] = {0.5,0.5,0.5};
	std::vector<float> weights (weights_array, weights_array + sizeof(weights_array) / sizeof(weights_array[0]) );
    
    Image color2graycastle = color2gray(im, weights);
    color2graycastle.write("./Output/graycastle.png");

    std::vector<Image> lumiChromiVect = lumiChromi(im);
    Image lumiImage = lumiChromiVect.at(0);
    lumiChromiVect.at(0).write("./Output/lumi.png");
    lumiChromiVect.at(1).write("./Output/chromi.png");

	Image recombined = brightnessContrastLumi(im, 1.5, 1.5, 0.5);
	recombined.write("./Output/new_recombined.png");
    // std::vector<Image> LC = lumiChromi(im);
    // LC[0].write("./Output/castle_luminance.png");
    // LC[1].write("./Output/castle_chrominance.png");
	
	Image zebra("./Input/zebra.png");
    std::vector<Image> spanishVect = spanish(zebra);
    spanishVect.at(0).write("./Output/spanish_color.png");
    spanishVect.at(1).write("./Output/spanish_gray.png");

	Image flower("./Input/flower.png");
    Image whitebalanced = grayworld(flower);
    whitebalanced.write("./Output/whitebalanced.png");
	
}
