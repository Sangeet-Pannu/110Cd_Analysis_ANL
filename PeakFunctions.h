#ifndef PEAKFUNCTIONS_H
#define PEAKFUNCTIONS_H

#include "TFitResult.h"
#include <TFile.h>
#include <TCanvas.h>
#include <TRandom.h>
#include <iostream>
#include <cmath>
#include <TROOT.h>
#include <TF1.h>
#include <TH1.h>

double SkewedGaussianWithBackground1(double *x, double *par) {
    double xx = x[0];
    double mean = par[0];
    double sigma = par[1];
    double alpha = par[2]; // skewness
    double A = par[3]; // amplitude
    double kL = par[4]; // left tail parameter
    double kR = par[5]; // right tail parameter
    double b0 = par[6]; // background intercept
    double b1 = par[7]; // background slope

    double delta = (xx - mean) / sigma;
    double absDelta = std::abs(delta);

    double gauss = A * exp(-0.5 * delta * delta);
    double leftTail = exp(-kL * absDelta);
    double rightTail = exp(-kR * absDelta);

    double peak;
    if (delta < 0) {
        peak = gauss * leftTail * (1 + erf(alpha * delta / sqrt(2)));
    } else {
        peak = gauss * rightTail * (1 + erf(alpha * delta / sqrt(2)));
    }

    double background = b0 + b1 * xx;

    return peak + background;
}

double SkewedGaussianWithBackground2(double *x, double *par) {
    double xx = x[0];
    // Parameters for the first skewed Gaussian
    double mean1 = par[0];
    double A1 = par[1]; // amplitude

    // Parameters for the second skewed Gaussian
    double mean2 = par[2];
    double A2 = par[3]; // amplitude

    // Shared shape parameters
    double sigma = par[4];
    double alpha = par[5]; // skewness
    double kL = par[6]; // left tail parameter
    double kR = par[7]; // right tail parameter

    // Background parameters
    double b0 = par[8]; // background intercept
    double b1 = par[9]; // background slope

    // Skewed Gaussian 1
    double delta1 = (xx - mean1) / sigma;
    double absDelta1 = std::abs(delta1);
    double gauss1 = A1 * exp(-0.5 * delta1 * delta1);
    double leftTail1 = exp(-kL * absDelta1);
    double rightTail1 = exp(-kR * absDelta1);
    double peak1 = (delta1 < 0) ? gauss1 * leftTail1 * (1 + erf(alpha * delta1 / sqrt(2)))
                                : gauss1 * rightTail1 * (1 + erf(alpha * delta1 / sqrt(2)));

    // Skewed Gaussian 2
    double delta2 = (xx - mean2) / sigma;
    double absDelta2 = std::abs(delta2);
    double gauss2 = A2 * exp(-0.5 * delta2 * delta2);
    double leftTail2 = exp(-kL * absDelta2);
    double rightTail2 = exp(-kR * absDelta2);
    double peak2 = (delta2 < 0) ? gauss2 * leftTail2 * (1 + erf(alpha * delta2 / sqrt(2)))
                                : gauss2 * rightTail2 * (1 + erf(alpha * delta2 / sqrt(2)));

    // Background
    double background = b0 + b1 * xx;
    
    return peak1 + peak2 + background;
}

double SkewedGaussianWithBackground2N2Sigma(double *x, double *par) {
    double xx = x[0];
    // Parameters for the first skewed Gaussian
    double mean1 = par[0];
    double A1 = par[1]; // amplitude

    // Parameters for the second skewed Gaussian
    double mean2 = par[2];
    double A2 = par[3]; // amplitude

    // Shared shape parameters
    double sigma = par[4];
    double sigma2 = par[5];
    double alpha = par[6]; // skewness
    double kL = par[7]; // left tail parameter
    double kR = par[8]; // right tail parameter

    // Background parameters
    double b0 = par[9]; // background intercept
    double b1 = par[10]; // background slope

    // Skewed Gaussian 1
    double delta1 = (xx - mean1) / sigma;
    double absDelta1 = std::abs(delta1);
    double gauss1 = A1 * exp(-0.5 * delta1 * delta1);
    double leftTail1 = exp(-kL * absDelta1);
    double rightTail1 = exp(-kR * absDelta1);
    double peak1 = (delta1 < 0) ? gauss1 * leftTail1 * (1 + erf(alpha * delta1 / sqrt(2)))
                                : gauss1 * rightTail1 * (1 + erf(alpha * delta1 / sqrt(2)));

    // Skewed Gaussian 2
    double delta2 = (xx - mean2) / sigma2;
    double absDelta2 = std::abs(delta2);
    double gauss2 = A2 * exp(-0.5 * delta2 * delta2);
    double leftTail2 = exp(-kL * absDelta2);
    double rightTail2 = exp(-kR * absDelta2);
    double peak2 = (delta2 < 0) ? gauss2 * leftTail2 * (1 + erf(alpha * delta2 / sqrt(2)))
                                : gauss2 * rightTail2 * (1 + erf(alpha * delta2 / sqrt(2)));

    // Background
    double background = b0 + b1 * xx;

    return peak1 + peak2 + background;
}

double SkewedGaussianWithBackground3(double *x, double *par) {
    double xx = x[0];

    // Parameters for the first skewed Gaussian
    double mean1 = par[0];
    double A1 = par[1]; // amplitude

    // Parameters for the second skewed Gaussian
    double mean2 = par[2];
    double A2 = par[3]; // amplitude

    // Parameters for the third skewed Gaussian
    double mean3 = par[4];
    double A3 = par[5]; // amplitude

    // Shared shape parameters
    double sigma = par[6];
    double alpha = par[7]; // skewness
    double kL = par[8]; // left tail parameter
    double kR = par[9]; // right tail parameter

    // Background parameters
    double b0 = par[10]; // background intercept
    double b1 = par[11]; // background slope

    // Skewed Gaussian 1
    double delta1 = (xx - mean1) / sigma;
    double absDelta1 = std::abs(delta1);
    double gauss1 = A1 * exp(-0.5 * delta1 * delta1);
    double leftTail1 = exp(-kL * absDelta1);
    double rightTail1 = exp(-kR * absDelta1);
    double peak1 = (delta1 < 0) ? gauss1 * leftTail1 * (1 + erf(alpha * delta1 / sqrt(2)))
                                : gauss1 * rightTail1 * (1 + erf(alpha * delta1 / sqrt(2)));

    // Skewed Gaussian 2
    double delta2 = (xx - mean2) / sigma;
    double absDelta2 = std::abs(delta2);
    double gauss2 = A2 * exp(-0.5 * delta2 * delta2);
    double leftTail2 = exp(-kL * absDelta2);
    double rightTail2 = exp(-kR * absDelta2);
    double peak2 = (delta2 < 0) ? gauss2 * leftTail2 * (1 + erf(alpha * delta2 / sqrt(2)))
                                : gauss2 * rightTail2 * (1 + erf(alpha * delta2 / sqrt(2)));

    // Skewed Gaussian 3
    double delta3 = (xx - mean3) / sigma;
    double absDelta3 = std::abs(delta3);
    double gauss3 = A3 * exp(-0.5 * delta3 * delta3);
    double leftTail3 = exp(-kL * absDelta3);
    double rightTail3 = exp(-kR * absDelta3);
    double peak3 = (delta3 < 0) ? gauss3 * leftTail3 * (1 + erf(alpha * delta3 / sqrt(2)))
                                : gauss3 * rightTail3 * (1 + erf(alpha * delta3 / sqrt(2)));

    // Background
    double background = b0 + b1 * xx;

    return peak1 + peak2 + peak3 + background;
}

double SkewedGaussianPeak(double *x, double *par) {
    double xx = x[0];
    double mean = par[0];
    double sigma = par[1];
    double alpha = par[2]; // skewness
    double A = par[3]; // amplitude
    double kL = par[4]; // left tail parameter
    double kR = par[5]; // right tail parameter

    double delta = (xx - mean) / sigma;
    double absDelta = std::abs(delta);

    double gauss = A * exp(-0.5 * delta * delta);
    double leftTail = exp(-kL * absDelta);
    double rightTail = exp(-kR * absDelta);

    double peak;
    if (delta < 0) {
        peak = gauss * leftTail * (1 + erf(alpha * delta / sqrt(2)));
    } else {
        peak = gauss * rightTail * (1 + erf(alpha * delta / sqrt(2)));
    }

    return peak;
}

double background(double *x,double* par){
    return (par[0]*x[0]) + par[1];
}

double backgroundWStep(double *x,double* par){

    double x0 = par[0]; //Step function off-set
    double beta = par[1]; //effective beta
    double Amp = par[2]; //Size of Step
   // double B = par[3];
    return Amp*(1.0/(1.0+exp(-beta*(x[0]-x0))));
}

double SkewedGaussianWithBackgroundStep(double *x, double *par) {
    double xx = x[0];
    // Parameters for the first skewed Gaussian
    double mean1 = par[0];
    double A1 = par[1]; // amplitude

    // Parameters for the second skewed Gaussian
    double mean2 = par[2];
    double A2 = par[3]; // amplitude

    // Shared shape parameters
    double sigma = par[4];
    double alpha = par[5]; // skewness
    double kL = par[6]; // left tail parameter
    double kR = par[7]; // right tail parameter

    // Background parameters
    double x0 = par[8]; // background intercept
    double beta = par[9]; // background slope
    double Amp = par[10];
    //double B = par[11];
    // Skewed Gaussian 1
    double delta1 = (xx - mean1) / sigma;
    double absDelta1 = std::abs(delta1);
    double gauss1 = A1 * exp(-0.5 * delta1 * delta1);
    double leftTail1 = exp(-kL * absDelta1);
    double rightTail1 = exp(-kR * absDelta1);
    double peak1 = (delta1 < 0) ? gauss1 * leftTail1 * (1 + erf(alpha * delta1 / sqrt(2)))
                                : gauss1 * rightTail1 * (1 + erf(alpha * delta1 / sqrt(2)));

    // Skewed Gaussian 2
    double delta2 = (xx - mean2) / sigma;
    double absDelta2 = std::abs(delta2);
    double gauss2 = A2 * exp(-0.5 * delta2 * delta2);
    double leftTail2 = exp(-kL * absDelta2);
    double rightTail2 = exp(-kR * absDelta2);
    double peak2 = (delta2 < 0) ? gauss2 * leftTail2 * (1 + erf(alpha * delta2 / sqrt(2)))
                                : gauss2 * rightTail2 * (1 + erf(alpha * delta2 / sqrt(2)));

    // Background
    double background = Amp*(1.0/(1+exp(-beta*(xx-x0))));;

    return peak1 + peak2 + background;
}


#endif