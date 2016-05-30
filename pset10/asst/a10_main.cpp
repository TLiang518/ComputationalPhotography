#include <iostream>
#include <cmath>
#include "a10.h"

using namespace std;


//Modelled after test functions in NPR_starter 
void testBrush(){
	/*
	Tests brush function.
	Visually confirm that produced image:
		- Has randomly placed different brush "strokes" of brush texture
		- In different colors
		- In same orientations
		- and should be in different places if you run it twice
	*/
	//Seed rand
	int nStrokes = 10;
	Image im(640,480,3);
	Image texture("Input/brush.png");

	vector<float> color {1.0,1.0,1.0};
	for (int i = 0; i < nStrokes; i++){
		for (int c = 0; c < im.channels(); c++) {
			color[c] = float(rand())/RAND_MAX;
		}
		int x = rand() % im.width();
		int y = rand() % im.height();
		brush(im, x, y, color, texture);
	}
	im.write("Output/testBrush_random.png");
}

//Modelled after test functions in NPR_starter 
void testAngle(){
	/*
	Tests computeAngles function
	*/
	
	int nStrokes = 10;
	Image im("Input/boston.png");
	Image texture("Input/brush.png");
	Image out = computeAngles(im);
	Image normalized_out = out/out.max();
	normalized_out.write("Output/angle_test_boston.png");
}


//from pset 7
void testComputeTensor() {
    // load images
    Image stata1("./Input/stata-1.png");
    Image stata2("./Input/stata-2.png");

    // compute tensors
    Image tensor1 = computeTensor(stata1);
    tensor1.write("./Output/stataTensor-1.png");
    float maxi = tensor1.max();
    if(maxi != 0) {
        tensor1 = tensor1 / maxi ;
        tensor1.write("./Output/stataTensor-1normed.png");
    }

    Image tensor2 = computeTensor(stata2);
    tensor2.write("./Output/stataTensor-2.png");
    maxi = tensor2.max();
    if(maxi != 0) {
        tensor2 = tensor2 / maxi ;
        tensor2.write("./Output/stataTensor-2normed.png");
    }
}


void testBrushRotate(int nAngles) {
	/*
	
	Tests brush function.
	Visually confirm that produced image:
		- Has randomly placed different brush "strokes" of brush texture
		- In different colors
		- In different orientations
		- and should be in different places if you run it twice
	*/

	Image brush("Input/brush.png");
	vector<Image> rotated_brushes = rotateBrushes(brush, nAngles);
	string filename_base = "Input/brush";
	for (int i  = 0; i < nAngles; i++){
		rotated_brushes[i].write(filename_base + "-" + to_string(i) + ".png");
	}

}

void testSharpnessMap(){
	Image im("Input/china.png");
	float sigma = 1.0;
	Image normalized_sharpness = sharpnessMap(im, sigma);
	normalized_sharpness.write("Output/testSharpnessMask_china.png");	
}

void testAnisotopicFilter(){
	Image im("Input/china.png");
	float sigma = 1.0;
	Image normalized_sharpness = anisotropic_gaussian(im, sigma);
	normalized_sharpness.write("Output/testAnisotopicFilter_china.png");	
}

void testSingleScalePaint(){
	Image im("Input/boston.png");
	Image texture("Input/brush.png");
	Image out = Image(im.width(), im.height(), im.channels());
	
	//Make importance image all ones
	Image importance = out;
	importance = importance + 1.0;

	//out = im;
	singleScalePaint(im, out, importance, texture, 10000, 30);

	out.write("Output/singleScalePaint_boston.png");
}


void testPainterly(){
	Image im("Input/cave_story.png");
	Image texture("Input/brush.png");
	Image out =	painterly(im, texture);
	out.write("Output/testPainterly_cave_story.png");
}

void testPainterly_Sephiroth(){
	Image im("Input/KHSeph.png");
	Image texture("Input/brush.png");
	Image out =	painterly(im, texture);
	out.write("Output/testPainterly_Sephiroth.png");
}

void testSingleScaleOrientedPaint(){
	Image im("Input/boston.png");
	Image texture("Input/brush.png");
	Image out = Image(im.width(), im.height(), im.channels());
	
	//Make importance image all ones
	Image importance = out;
	importance = importance + 1.0;

	//out = im;
	cout << "Before singleScaleOrientedPaint" << endl;
	singleScaleOrientedPaint(im, out, importance, texture);

	out.write("Output/singleScaleOrientedPaint_boston.png");
}

void testOrientedPaint_EdgeAlignment(){
	Image horizontal_stripe("Input/horizontal_stripe.png");
	Image vertical_stripe("Input/vertical_stripe.png");
	Image angled_stripe("Input/colored_angled_stripes.png");

	Image texture("Input/longBrush.png");

	Image horizontal_stripe_out = orientedPaint(horizontal_stripe, texture, 5000, 20, .1);
	horizontal_stripe_out.write("Output/testOrientedPaint_EdgeAlignment_horizontal_stripe.png");
	
	Image vertical_stripe_out = orientedPaint(vertical_stripe, texture, 5000, 20, .1);
	vertical_stripe_out.write("Output/testOrientedPaint_EdgeAlignment_vertical_stripe.png");

	Image angled_stripe_out = orientedPaint(angled_stripe, texture, 5000, 20, .1);
	angled_stripe_out.write("Output/testOrientedPaint_EdgeAlignment_angled_stripe.png");
}

void testOrientedPaint_CrossStitchEdgeAlignment(){
	Image horizontal_stripe("Input/horizontal_stripe.png");
	Image vertical_stripe("Input/vertical_stripe.png");
	Image angled_stripe("Input/colored_angled_stripes.png");

	Image texture("Input/longBrush.png");

	Image horizontal_stripe_out = orientedPaint(horizontal_stripe, texture, 5000, 20, .1, false);
	horizontal_stripe_out.write("Output/testOrientedPaint_CrossStitchEdgeAlignment_horizontal_stripe.png");
	
	Image vertical_stripe_out = orientedPaint(vertical_stripe, texture, 5000, 20, .1, false);
	vertical_stripe_out.write("Output/testOrientedPaint_CrossStitchEdgeAlignment_vertical_stripe.png");

	Image angled_stripe_out = orientedPaint(angled_stripe, texture, 5000, 20, .1, false);
	angled_stripe_out.write("Output/testOrientedPaint_CrossStitchEdgeAlignment_angled_stripe.png");
}


void testOrientedPaint(){
	Image im("Input/sunset.png");
	Image texture("Input/brush.png");

	Image out = orientedPaint(im, texture, 40000, 40, .1, true);
	out.write("Output/testOrientedPaint_sunset_OriginalStitch.png");
	
	//out = orientedPaint(im, texture, 40000, 40, .1, false);
	//out.write("Output/testOrientedPaint_stanford_CrossStitch.png");
}

int main()
{
	srand (time(NULL));
    // Test your intermediate functions
    //testSomeFunction();
    /*
    testBrush();
    testComputeTensor();
    testAngle();
    testSharpnessMap();
    testSingleScalePaint();
    testBrushRotate(10);
    testPainterly();
    testPainterly_Sephiroth();


    testSingleScaleOrientedPaint();
    testAnisotopicFilter();
    testOrientedPaint_EdgeAlignment();
    testOrientedPaint_EdgeAlignment();
    testOrientedPaint_CrossStitchEdgeAlignment();
    */
    testOrientedPaint();
    
    return EXIT_SUCCESS;
}
