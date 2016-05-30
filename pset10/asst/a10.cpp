#include <iostream>
#include <random>
#include "a10.h"

using namespace std;


void brush(Image & im, int x, int y, vector<float> color, Image &texture) {
	/*
	out: the image to draw to.
	y,x: where to draw in out.
	color: the color of the stroke.
	texture: the texture of the stroke.

	Write a function brush(out, y, x, color, texture) 
	that takes as input a mutable image out and draws
	(“splats”) a single brush stroke centered at y, x. 
	
	The appearance of the brush is specified by an 
	opacity texture and a 3-array color.
	*/

	//bounds checking

	if (x > texture.width()/2 && x < (im.width() - texture.width()/2)){
		if (y > texture.height()/2 && y < (im.height() - texture.height()/2)){
			//inner texture loop
			for (int texture_x = 0; texture_x < texture.width(); texture_x++){
				for (int texture_y = 0; texture_y < texture.height(); texture_y++){
					for (int c = 0; c < color.size(); c++){
						int new_x = x + texture_x - texture.width()/2;
						int new_y = y + texture_y - texture.height()/2;
						im(new_x, new_y, c) = im(new_x, new_y, c)*(1 - texture(texture_x, texture_y, c)) + color[c]*texture(texture_x, texture_y, c);
					}
				}
			}
		}
	}
}

void singleScalePaint(Image & im, Image & out, Image & importance, Image & texture, int N, int size, float noise){
	/* initialize random seed so we can use it throughout*/
  	srand (time(NULL));

	//First, scale the texture image so that it has maximum size size. Use the provided scaleImage function or your own method.
	float scale_factor = float(size)/max(texture.height(), texture.width());
	Image scaled_texture = scaleLin(texture, scale_factor);

	/*
	For each of N random locations y,x splat a brush in out using the above function.
	Generate random locations using rnd.randrange(start, stop), as- suming you imported
	the random module using import random as rnd; or using numpy.random.randint(low, high).
	*/

	//Create color vector
	vector<float> color(im.channels(), 0.0);

	//Since we reject a number of samples, we do not splat the 
	//required N strokes. In order to fix this, multiply N by 
	//a normalization factor based on the average probability
	//of accepting samples.
	int num_iterations = N/importance.mean();
	cout << "num_iterations is: " << num_iterations << endl;

	for (int i = 0; i < num_iterations; i++){
		
		 int x = rand() % im.width();
    	 int y = rand() % im.height();

		 int r = float(rand())/RAND_MAX;

		 if (r < importance(x,y)) {
		 	for (int c = 0; c < im.channels(); c++){
		 		//noise formula (given)
		 		//read in color at y,x,c
		 		//(1 - noise/2 + 3*n*noise) equivalent to (1-noise/2+noise*numpy(random*rand*3))
		 		float n = float(rand())/RAND_MAX;
		 		color[c] = im(x,y,c)*(1 - noise/2 + n*noise);
		 	}
		 	brush(out, x, y, color, scaled_texture);
		 }
	}
}

Image sharpnessMap(Image &im, float sigma){
	Image lumi = lumiChromi(im)[0];
	Image blurred_lumi = gaussianBlur_separable(lumi, sigma);
	Image high_freq_lumi = lumi - blurred_lumi;
	Image lumi_energy = high_freq_lumi*high_freq_lumi;
	Image sharpness = gaussianBlur_separable(lumi_energy, 4.0*sigma);
	Image normalized_sharpness = sharpness/sharpness.max();
	return normalized_sharpness;
}

Image anisotropic_gaussian(Image &im, float sigma){
	Image lumi = lumiChromi(im)[0];

	float angle_to_rotate = 60 * M_PI / 180;
	cout << "angle_to_rotate is: " << angle_to_rotate << endl;
	Image blurred_lumi = gaussianBlur_horizontal(lumi, sigma);

	Image rotated_blurred_lumi = rotate(blurred_lumi, angle_to_rotate);
	Image filtered_rotated_blurred_lumi = gaussianBlur_horizontal(rotated_blurred_lumi, sigma);
	
	cout << "image size: " << im.width() << " " << im.height() << endl;
	cout << "filtered_rotated_blurred_lumi size: " << im.width() << " " << im.height() << endl;

	Image rotated_lumi = rotate(lumi, angle_to_rotate);
	Image high_freq_lumi_final = rotated_lumi - filtered_rotated_blurred_lumi;
	
	Image final = rotate(high_freq_lumi_final, -angle_to_rotate);
	Image normalized_anisotropic_lumi = final/final.max();
	return normalized_anisotropic_lumi;
}

Image painterly(Image &im, Image &texture, int N, int size, float noise){
	Image out(im.width(), im.height(), im.channels());

	//First pass: use brushes of size size, in second pass: use brushes of size size/4
	//First pass should use constant importance map

	Image first_pass_importance(im.width(), im.height(), im.channels());
	first_pass_importance = first_pass_importance + 1.0;

	//Second finer pass should only add strokes where the image has strong high frequencies
	Image second_pass_importance(im.width(), im.height(), im.channels());
	second_pass_importance = sharpnessMap(im);

	singleScalePaint(im, out, first_pass_importance, texture, N, size, noise);
	singleScalePaint(im, out, second_pass_importance, texture, N, size/4, noise);
	return out;
}

/////////////////////////////////////////////
// Orderly Painting Portion
////////////////////////////////////////////


Image computeTensor(const Image &im, float sigmaG, float factorSigma) {
    // // --------- HANDOUT  PS07 ------------------------------
    // Compute xx/xy/yy Tensor of an image. (stored in that order)
    cout << "in computeTensor" << endl;
    vector<Image> lumi_chromi = lumiChromi(im);
    Image lumi = lumi_chromi[0];
    Image chromi = lumi_chromi[1];
    //Using a Gaussian with standard deviation sigmaG, blur the luminance
    //to control the scale at which corners are extracted. A little bit
    //of blur helps smooth things out and help extract stable mid-scale corners.
    Image blurred_lumi = gaussianBlur_separable(lumi, sigmaG);
    Image gradientX_lumi = gradientX(blurred_lumi);
    Image gradientY_lumi = gradientY(blurred_lumi);

    //Structure tensor image
    //Where channel 0 is Ix2
    //and channel 1 is IxIy
    //Channel 2 is Iy2
    Image perPixelContributions(im.width(), im.height(), 3);

    for (int i = 0; i < im.width(); i++){
        for (int j = 0; j < im.height(); j++){
            perPixelContributions(i,j,0) = pow(gradientX_lumi(i,j),2);
            perPixelContributions(i,j,1) = gradientX_lumi(i,j)*gradientY_lumi(i,j);
            perPixelContributions(i,j,2) = pow(gradientY_lumi(i,j),2);
        }
    }

    Image structure_tensor = gaussianBlur_separable(perPixelContributions, sigmaG*factorSigma);

    return structure_tensor;
}

Image computeAngles(Image &im){
	/*
	Takes an image and returns a new image of the same size where each pixel
	has been replaced by the angle between its local edge orientation and the
	horizotanl line.

	Didn't we do this for the world panoramas?
	*/
	cout << "In computeAngles" << endl;
	Image out(im.width(), im.height(), 1);
	Image tensor = computeTensor(im);

	for (int y = 0; y < im.height(); y++){
		for (int x = 0; x < im.width(); x++){
			//Per piazza suggestion of just calculating eigenvalues and vectors for a 2x2 matrix instead of using eigen
			//http://www.math.harvard.edu/archive/21b_fall_04/exhibits/2dmatrices/index.html
			/*
			For A = [a b; c d]
			*/

			float a = tensor(x,y,0);
			float b = tensor(x,y,1);
			float c = b;
			float d = tensor(x,y,2);

			//Calculate Trace and Determinant
			float trace = a + d;
			float det = a*d - b*c;

			float eigenvalue1 = trace/2 + sqrt(trace*trace/4 - det);
			float eigenvalue2 = trace/2 - sqrt(trace*trace/4 - det);

			//Where these are both column vectors
			vector<float> eigenvector1{0,0};
			vector<float> eigenvector2{0,0};
			vector<float> min_eigenvector{0,0};
			//if c is not zero, eigenvectors are:
			if (fabs(eigenvalue1) < fabs(eigenvalue2)){
				min_eigenvector[0] = eigenvalue1 - d;
				min_eigenvector[1] = c;
			}
			else {
				min_eigenvector[0] = b;
				min_eigenvector[1] = eigenvalue2 - a;
			}

			//Eigenvector corresponding to minimum eigenvalue
			//vector<float> min_eigenvector = (fabs(eigenvalue1) < fabs(eigenvalue2)) ? eigenvector1 : eigenvector2;

			//atan2 is (y,x) not (x,y)
			float angle_to_horizontal = atan2(min_eigenvector[1],min_eigenvector[0]);

			//Now set the pixel in the out image to this calculated angle to horizontal
			angle_to_horizontal = fmod(angle_to_horizontal + 2*M_PI, 2*M_PI);

			out(x,y) = angle_to_horizontal;

			}
		}

	return out;

}

vector<Image> rotateBrushes(Image &texture, int n){
	vector<Image> out;
	float theta;

	for (int i = 0; i < n; i++){
		theta = -M_PI*i/n;
		out.push_back(rotate(texture, theta));
	}

	return out;
}

//Rotates brushes 90deg to the direciton of an edge, producing a cross-stitch-like effect
vector<Image> rotateBrushesForCrossStitch(Image &texture, int n){
	cout << "in rotateBrushes for cross stitch" << endl;
	vector<Image> out;
	float theta;

	for (int i = 0; i < n; i++){
		theta = -M_PI*i/n + M_PI/2;
		out.push_back(rotate(texture, theta));
	}

	return out;
}

void singleScaleOrientedPaint(Image & im, Image & out, Image & importance, Image & texture, int N, int size, float noise, int nAngles, bool useRegularStroke){
	/* initialize random seed so we can use it throughout*/
  	srand (time(NULL));

	//First, scale the texture image so that it has maximum size size. Use the provided scaleImage function or your own method.
	float scale_factor = float(size)/max(texture.height(), texture.width());
	Image scaled_texture = scaleLin(texture, scale_factor);

	/*
	For each of N random locations y,x splat a brush in out using the above function.
	Generate random locations using rnd.randrange(start, stop), as- suming you imported
	the random module using import random as rnd; or using numpy.random.randint(low, high).
	*/

	//Since we reject a number of samples, we do not splat the 
	//required N strokes. In order to fix this, multiply N by 
	//a normalization factor based on the average probability
	//of accepting samples.
	int num_iterations = N/importance.mean();
	cout << "num_iterations is: " << num_iterations << endl;


	//Get the angle images
	Image angle_image = computeAngles(im);

	//Generate the rotated texture images
	vector<Image> scaled_textures = (useRegularStroke == true) ? rotateBrushes(scaled_texture, nAngles) : rotateBrushesForCrossStitch(scaled_texture, nAngles);
	
	cout << "after rotate brushes. Used regular brush texture? " << (useRegularStroke == true) << endl;
	//Create color vector
	vector<float> color(im.channels(), 0.0);


	for (int i = 0; i < num_iterations; i++){
		
		 int x = rand() % im.width();
    	 int y = rand() % im.height();

		 int r = float(rand())/RAND_MAX;

		 if (r < importance(x,y)) {
		 	for (int c = 0; c < im.channels(); c++){
		 		//noise formula (given)
		 		//read in color at y,x,c
		 		//(1 - noise/2 + 3*n*noise) equivalent to (1-noise/2+noise*numpy(random*rand*3))
		 		float n = float(rand())/RAND_MAX;
		 		color[c] = im(x,y,c)*(1 - noise/2 + n*noise);
		 	}
		 	int texture_index = int(round(nAngles*angle_image(x,y)/(M_PI))) % nAngles;
		 	brush(out, x, y, color, scaled_textures[texture_index]);
		 }
	}
}


Image orientedPaint(Image &im, Image &texture, int N, int size, float noise, bool useRegularStroke){
	Image out(im.width(), im.height(), im.channels());

	//First pass: use brushes of size size, in second pass: use brushes of size size/4
	//First pass should use constant importance map

	Image first_pass_importance(im.width(), im.height(), 1);
	first_pass_importance = first_pass_importance + 1.0;

	//Second finer pass should only add strokes where the image has strong high frequencies
	Image second_pass_importance(im.width(), im.height(), 1);
	second_pass_importance = sharpnessMap(im);

	cout << "Before singleScaleOrientedPaint round 1" << endl;
	singleScaleOrientedPaint(im, out, first_pass_importance, texture, N, size, noise, 36, useRegularStroke);
	cout << "Before singleScaleOrientedPaint round 2" << endl;
	singleScaleOrientedPaint(im, out, second_pass_importance, texture, N, size/4, noise, 36, useRegularStroke);

	return out;
}