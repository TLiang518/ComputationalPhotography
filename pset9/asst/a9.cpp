#include <a9.h>
#include <algorithm>
#include "timing.h"

using namespace Halide;

Func HarrisCorners(const Image<float> &input,
    int schedule_index,
    vector<int> schedule_parameters,
    float sigmaG,
    float factorSigma,
    float truncate,
    float k,
    float thresh,
    int maxiDiam
    ) 
{

    // --------- HANDOUT  PS08 ------------------------------
    // Re-implement the Harris Corner detector using Halide. By fusing
    // different stages of the processing pipeline and properly scheduling your
    // opertations, you should be able to get a drastic speed-up over the
    // reference implementation given in 'reference_implementation.cpp'.
    //
    // Make sure you debug your intermediate Func (by
    // realizing them over some domain for example). You can use the
    // normalize_image helper function to remap the pixel value to the [0,1]
    // range, valid for display.
    //
    // Start by building a correct implementation of the algorithm,
    // Then walk your way up to a fast implementation by trying
    // various schedules.
    //
    // Return the output Func, properly scheduled. Chek the reference
    // implementation to see how the input parameters are used.
    //
    // You can use scheule_index to switch between different schedules and
    // compare them, this will come in handy when you implement an autotuner
    // (6.865). You can also pass schedule parameters as a vector of ints.
    

    // --------- Algorithm ----------------------------------
    Var x("x");
    Var y("y");
    Var c("c");

    std::cout << "IN HARRIS CORNERS" << endl;

    int w = input.width();
    int h = input.height();
    int z = input.channels();

    //Lumi
    Func lumi("lumi");
    lumi(x,y) = 0.3f*input(x, y, 0) + 0.6f*input(x, y, 1)+ 0.1f*input(x, y, 2);

    ///////////////////////////////////////////////////////////////////////////////////

    //Gaussian Blur
    std::cout << "BEFORE GAUSSIAN BLUR" << endl;

    Func GB("Gaussian");

    int radius = sigmaG * truncate;
    int fwidth = 2 * radius + 1;

    // compute the Gaussian kernel
    Func GKernelUnNorm("GKernelUnnorm");
    Func GKernelSum   ("GKernelSum");
    Func GKernel      ("GKernel");
    RDom gaussian_rx(0, fwidth);
    GKernelUnNorm(x) = exp(-((x-radius)*(x-radius))/(2.0f*sigmaG*sigmaG));
    GKernelSum   (x) = sum(GKernelUnNorm(gaussian_rx));
    GKernel      (x) = GKernelUnNorm(x)/GKernelSum(0);

    // Schedule the kernel to be computed root, so we do not recompute it every
    // time we use it.
    GKernelUnNorm.compute_root();
    GKernelSum   .compute_root();
    GKernel      .compute_root();
    Image<float> kernel = GKernel.realize(fwidth);

    //Load image
    //Use lumi

    // Clamp the image (boundary conditions)
    Func gaussian_clamped("gaussian_clamped");
    gaussian_clamped(x,y) = lumi(clamp(x,0,input.width() - 1), clamp(y, 0, input.height() - 1));
    
    // Blur in x
    Func gauss_x("gauss_x");
    gauss_x(x, y) = sum(GKernel(gaussian_rx.x)*gaussian_clamped(x + gaussian_rx.x - radius, y));

    // Blur in y
    Func gauss_y("gauss_y");
    gauss_y(x, y) = sum(GKernel(gaussian_rx.x)*gauss_x(x, y + gaussian_rx.x - radius));

    GB(x, y) = gauss_y(x,y);

    std::cout << "AFTER GAUSSIAN BLUR" << endl;

    ///////////////////////////////////////////////////////////////////////////////////

    // Sobel X and Y

    std::cout << "BEFORE SOBEL" << endl;

    Func clamped("clamped");
    clamped(x,y) = GB(clamp(x,0,input.width() - 1), clamp(y, 0, input.height() -1));

    Func sobelx("sobelx");
    Func sobely("sobely");

    sobelx(x,y) = -clamped(x-1,y-1) - 2.0f*clamped(x-1,y) - clamped(x-1,y+1) + clamped(x+1, y-1) + 2*clamped(x+1,y) + clamped(x+1,y+1);    
    sobely(x,y) = -clamped(x-1,y-1) - 2.0f*clamped(x,y-1) - clamped(x+1,y-1) + clamped(x-1,y+1) + 2*clamped(x,y+1) + clamped(x+1,y+1);

/*    Image<float> outputx;
    Image<float> outputy;
    outputx = sobelx.realize(w,h);
    outputy = sobely.realize(w,h);
    
    Image<float> normalized_outputx = normalize_image(outputx);
    Image<float> normalized_outputy = normalize_image(outputy);

    save_image(normalized_outputx, "Output/outputx.png");
    save_image(normalized_outputy, "Output/outputy.png");*/

    std::cout << "AFTER SOBEL" << endl;


    /////////////////////////////////////////////////// 
   
    std::cout << "BEFORE TENSOR" << endl;

    //Tensor
    Func tensor("tensor");
    tensor(x,y,c) = 0.0f;
    tensor(x,y,0) = sobelx(x,y)*sobelx(x,y);
    tensor(x,y,1) = sobelx(x,y)*sobely(x,y);
    tensor(x,y,2) = sobely(x,y)*sobely(x,y);
    
    std::cout << "AFTER TENSOR" << endl;

/*    Image<float> halide_tensor;
    halide_tensor = tensor.realize(w,h,z);
    
    Image<float> normalized_halide_tensor = normalize_image(halide_tensor);
    save_image(normalized_halide_tensor, "Output/normalized_halide_tensor.png");*/

    //Blurred Tensor

    ////////////////////////////////////////////////////

    // Gaussian Blur #2

    std::cout << "before tensor gaussian blur" << endl;

    float tensor_sigma = sigmaG*factorSigma;
    int tensor_radius = tensor_sigma*truncate;
    int tensor_fwidth = 2 * tensor_radius + 1;

    // compute the tensor_tensor_gaussian kernel
    Func tensor_GKernelUnNorm("tensor_GKernelUnnorm");
    Func tensor_GKernelSum   ("tensor_GKernelSum");
    Func tensor_GKernel      ("tensor_GKernel");
    RDom tensor_gaussian_rx(0, tensor_fwidth);
    tensor_GKernelUnNorm(x) = exp(-((x-tensor_radius)*(x-tensor_radius))/(2.0f*tensor_sigma*tensor_sigma));
    tensor_GKernelSum   (x) = sum(tensor_GKernelUnNorm(tensor_gaussian_rx));
    tensor_GKernel      (x) = tensor_GKernelUnNorm(x)/tensor_GKernelSum(0);

    // Schedule the kernel to be computed root, so we do not recompute it every
    // time we use it.
    tensor_GKernelUnNorm.compute_root();
    tensor_GKernelSum   .compute_root();
    tensor_GKernel      .compute_root();
    Image<float> tensor_kernel = tensor_GKernel.realize(tensor_fwidth);

    //Load image
    //Use lumi

    // Clamp the image (boundary conditions)
    Func tensor_gaussian_clamped("tensor_gaussian_clamped");
    tensor_gaussian_clamped(x,y,c) = tensor(clamp(x,0,input.width() - 1), clamp(y, 0, input.height() - 1), c);

    // Blur in x
    Func tensor_gauss_x("tensor_gauss_x");
    tensor_gauss_x(x, y, c) = sum(tensor_GKernel(tensor_gaussian_rx.x)*tensor_gaussian_clamped(x + tensor_gaussian_rx.x - tensor_radius, y, c));

    // Blur in y
    Func tensor_gauss_y("tensor_gauss_y");
    tensor_gauss_y(x, y, c) = sum(tensor_GKernel(tensor_gaussian_rx.x)*tensor_gauss_x(x, y + tensor_gaussian_rx.x - tensor_radius, c));

    Func blurred_tensor("blurred_tensor");
    blurred_tensor(x, y, c) = tensor_gauss_y(x,y,c);

    cout << "After tensor gaussian blur" << endl;

    /////////////////////////////////////////////////////

    //Response

    cout << "Before response" << endl;

    Func response("response");
    response(x,y) = blurred_tensor(x,y,0) * blurred_tensor(x,y,2) - blurred_tensor(x,y,1) * blurred_tensor(x,y,1) - k * ((blurred_tensor(x,y,0) + blurred_tensor(x,y,2)) * (blurred_tensor(x,y,0) + blurred_tensor(x,y,2)));

/*    Image<float> output_response;
    output_response = response.realize(w,h);
    
    Image<float> normalized_output_response = normalize_image(output_response);
    save_image(normalized_output_response, "Output/normalized_output_response.png");*/

    cout << "After response" << endl;

    ////////////////////////////////////////////////////

    //Maximum Filter
    /*
    RDom r1(0, 10), r2(0, 10), r3(0, 10);
    Expr e = minimum(r1, product(r2, sum(r3, r1 + r2 + r3)));
    */

    Func maxi_filter("maxi_filter");
    RDom maxi_rdom(-maxiDiam/2, maxiDiam, -maxiDiam/2, maxiDiam,"maxi_rdom");
    maxi_filter(x,y) = maximum(maxi_rdom, response(clamp(x + maxi_rdom.x, 0, w-1), clamp(y + maxi_rdom.y, 0, h-1)));

/*    Image<float> output_maxi_filter;
    output_maxi_filter = maxi_filter.realize(w,h);
    
    Image<float> normalized_maxi_filter = normalize_image(output_maxi_filter);
    save_image(normalized_maxi_filter, "Output/normalized_maxi_filter.png");*/

    ////////////////////////////////////////////////////

    // Non-max. suppression

    Func output("output");
    output(x,y) = select(response(x,y) == maxi_filter(x,y) && response(x,y) > thresh, 1.0f, 0.0f);
/*
    Image<float> output_final;
    output_final = output.realize(w,h);
    
    Image<float> normalized_output_final = normalize_image(output_final);
    save_image(normalized_output_final, "Output/normalized_output_final.png");*/


    // --------- Schedule -----------------------------------
    // Useful structure for you to switch between various schedule
    // pass schedule_index as argument to your HarrisCorners.
    // =========
    // IMPORTANT
    // =========
    // At the end, put your best schedule in the "schedule_index=0" case.
    
    // declare new Var that will be useful when we get to tiling for the indices of the tile and for pixel indices inside the tile
    Var tx("tx"), ty("ty"), xi("xi"), yi("yi"); 

    switch(schedule_index){
        case 0:
        {

            output.tile(x, y, tx, ty, xi, yi, 256, 32).parallel(ty).vectorize(xi, 16); 

            cout << "IN MY CASE" << endl;

            maxi_filter.compute_at(output, tx);
            maxi_filter.vectorize(x, 16);

            response.compute_at(output, tx);
            response.vectorize(x, 16);

            blurred_tensor.compute_at(output, tx);
            blurred_tensor.vectorize(x, 16);

            tensor_gauss_x.compute_at(output, tx);
            tensor_gauss_x.vectorize(x,16);

            tensor_gauss_y.compute_at(output, tx);
            tensor_gauss_y.vectorize(x,16);

            tensor_gaussian_clamped.compute_at(output, tx);
            tensor_gaussian_clamped.vectorize(x, 16);

            tensor.compute_at(output, tx);
            tensor.vectorize(x, 16);

            sobelx.compute_at(output, tx);
            sobelx.vectorize(x, 16);
            sobely.compute_at(output, tx);
            sobely.vectorize(x, 16);

            GB.compute_at(output, tx);
            GB.vectorize(x, 16);

            gauss_x.compute_at(output, tx);
            gauss_x.vectorize(x, 16);
            gauss_y.compute_at(output, tx);
            gauss_y.vectorize(x, 16);
            gaussian_clamped.compute_at(output, tx);
            gaussian_clamped.vectorize(x, 16);
            
            lumi.compute_at(output, tx);
            lumi.vectorize(x, 16);
            clamped.compute_at(output, tx);
            clamped.vectorize(x, 16);

            break;
          
        }
        
        default:
        {
            // The default, all root schedule
            // run HarrisCorners(input, -1) for example to get there.
            apply_auto_schedule(output);
            break;
        }
    }

    return output;
}

// -------------------------------------------------------------------------------------
// This applies a compute_root() schedule to all the Func's that are consumed
// by the calling Func. DO NOT EDIT IT.
void apply_auto_schedule(Func F) {
    map<string,Internal::Function> flist = Internal::find_transitive_calls(F.function());
    flist.insert(std::make_pair(F.name(), F.function()));
    map<string,Internal::Function>::iterator fit;
    for (fit=flist.begin(); fit!=flist.end(); fit++) {
        Func f(fit->second);
        f.compute_root();
        // cout << "Warning: applying default schedule to " << f.name() << endl;
    }
    cout << endl;
}


void apply_my_schedule(Func F) {
    map<string,Internal::Function> flist = Internal::find_transitive_calls(F.function());
    flist.insert(std::make_pair(F.name(), F.function()));
    map<string,Internal::Function>::iterator fit;
    for (fit=flist.begin(); fit!=flist.end(); fit++) {
        Func f(fit->second);
        f.compute_root();
        cout << "Applying my schedule to " << f.name() << endl;
    }
    cout << endl;
}

// -------------------------------------------------------------------------------------

void autoschedule_harris(const Image<float> &input) {
    // --------- HANDOUT  PS08 ------------------------------
    // An example autotuner, that does not do much yet... Build your own !
    //

    cout << "** Autotuning begins **" << endl
         << "=======================" << endl;
    int w = input.width();
    int h = input.height();

    vector<float> timings;
    vector<int> schedule_idxs;
    vector<vector<int >> schedule_params;
    
    vector<int> params;

    for (int idx = 1; idx < 20; ++idx) {
        for (int idy = 1; idy < 20; ++idy){
            params.clear();
            params.push_back(idx*16);
            params.push_back(idy*16);

            Func f = HarrisCorners(input, 0, params);

            cout << endl 
                 << "---------------------------------------" << endl;
            cout 
                << "* Autoschedule | " <<  idx << endl;

            cout << "  Params: " ;
            for (int i = 0; i < params.size(); ++i) {
                cout << params[i] << " ";
            }
            cout << endl;

            float time = profile(f, w, h);
            timings.push_back(time);
            schedule_idxs.push_back(idx);
            schedule_params.push_back(params);

            cout << endl 
                 << "---------------------------------------" << endl;
        }            
    }

    auto mini = std::min_element(timings.begin(), timings.end());
    int min_idx = mini-timings.begin();
    cout << endl
         << "==========================" << endl
         << "** Autotuning complete **" << endl
         << "   - best schedule: " << min_idx << endl
         << "   - time: " << timings[min_idx] << " ms" << endl
         << "   - params: ";
    for (int i = 0; i < schedule_params[min_idx].size(); ++i) {
        cout << schedule_params[min_idx][i] << " ";
    }
    cout << endl;
}
