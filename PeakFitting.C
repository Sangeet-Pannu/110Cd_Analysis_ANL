#include <TROOT.h>
#include <TF1.h>
#include <TH1.h>
#include <TFile.h>
#include <TCanvas.h>
#include <TRandom.h>
#include <iostream>
#include <cmath>

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

double background(double *x,double* par){
    return (par[0]*x[0]) + par[1];
}

void FitSkewedGaussianPeak() {
    int numPeaks;
    TFile *file = TFile::Open("152EuCal.root");
    if (!file || file->IsZombie()) {
        std::cerr << "Error opening file." << std::endl;
        return;
    }
    TFitResultPtr fitResult;
    TH1D *hist = nullptr;
    file->GetObject("gE", hist);
    if (!hist) {
        std::cerr << "Histogram 'LabEnoTG_C' not found." << std::endl;
        return;
    }

    TCanvas *c1 = new TCanvas("c1", "Fit Skewed Gaussian with Two Tails", 800, 600);
    hist->Draw();
    c1->Modified();
    c1->Update();
  //  gPad->WaitPrimitive();
  //  gPad->WaitPrimitive();
    double fitMin, fitMax, peakCentroid;
    std::cout << "Enter fitting domain (min max): ";
    std::cin >> fitMin >> fitMax;
      std::cout << "Enter the number of peaks to fit (1 or 2): ";
        std::cin >> numPeaks;

        TF1 *fitFunc;
        if (numPeaks == 1) {
            fitFunc = new TF1("fitFunc", SkewedGaussianWithBackground1, fitMin, fitMax, 8);
            double peakCentroid;
            std::cout << "Enter peak centroid: ";
            std::cin >> peakCentroid;
            fitFunc->SetParameters(peakCentroid, 1, 0, 10000, 1, 1, hist->GetMinimum(), 0);
            fitFunc->SetParNames("Centroid","Sigma","Skew","amplitude","left tail par","right tail par","linear Intercept","Linear Slope");
            fitFunc->SetParLimits(0, peakCentroid-3, peakCentroid+3);
            fitFunc->SetParLimits(1, 1, 10);
            fitFunc->SetParLimits(3, 0, 1e+10);
        } else if (numPeaks == 2) {
            fitFunc = new TF1("fitFunc", SkewedGaussianWithBackground2, fitMin, fitMax, 10);
            double peakCentroid1, peakCentroid2;
            std::cout << "Enter first peak centroid: ";
            std::cin >> peakCentroid1;
            std::cout << "Enter second peak centroid: ";
            std::cin >> peakCentroid2;
            fitFunc->SetParameters(peakCentroid1, hist->GetMaximum(), peakCentroid2, hist->GetMaximum() / 2,
                                   1, 0, 1, 1, hist->GetMinimum(), 0);
            fitFunc->SetParNames("Centroid 1 ","amplitude 1","Centroid 2","amplitude 2","Sigma overall","skewness","left tail par","right tail par","linear Intercept","Linear Slope");
            fitFunc->SetParLimits(0, peakCentroid1-3, peakCentroid1+3); 
            fitFunc->SetParLimits(2, peakCentroid2-3, peakCentroid2+3);
            fitFunc->SetParLimits(4, 1, 7);
        }
    auto fitHistogram = [&](bool userSetLimits = false) {
        if (userSetLimits) {
            std::cout << "(Enter 1 or 2 for limiting parameters or fixing them or -1 to skip), For 1: Enter parameter limits (low high) for each parameter, For 2: Enter fixed parameter.:\n";
            for (int i = 0; i <10; ++i) {
                double low, high;
                int parW;
                std::cout << "Parameter " << i << ": ";
                std::cin >> parW >> low >> high;
                if(parW == -1) continue;
                if(parW==1)
                  {  
                    fitFunc->SetParLimits(i, low, high);        
                  }
                else{
                    fitFunc->FixParameter(i,low);
                }

            }
        }
        hist->SetAxisRange(fitMin-25,fitMax+25);
        hist->Draw();
        fitResult = hist->Fit(fitFunc, "SREM");

        TF1 *bkfunc = new TF1("background", background, fitMin, fitMax, 2);
        bkfunc->SetParameter(0,fitFunc->GetParameter(7));
        bkfunc->SetParameter(1,fitFunc->GetParameter(6));
        bkfunc->SetLineColor(kRed);
        bkfunc->Draw("same");
        c1->Modified();
        c1->Update();
    };

    fitHistogram();

    char refit;
    do {
        std::cout << "Do you want to refit with parameter limits? (y/n): ";
        std::cin >> refit;
        if (refit == 'y') {
            fitHistogram(true);
        }
    } while (refit == 'y');
    auto covMatrix = fitResult->GetCovarianceMatrix();


    if(numPeaks==1)
    {
        double peakArea = fitFunc->Integral(fitMin, fitMax) - (fitFunc->GetParameter(6) * (fitMax - fitMin) + fitFunc->GetParameter(7) * 0.5 * (fitMax * fitMax - fitMin * fitMin));
        double peakAreaError = fitFunc->IntegralError(fitMin, fitMax, fitFunc->GetParameters(), covMatrix.GetMatrixArray());
        double centroid = fitFunc->GetParameter(0);
        double centroidError = fitFunc->GetParError(0);

        std::cout << "Peak Area: " << peakArea << " ± " << peakAreaError << "\n";
        std::cout << "Centroid: " << centroid << " ± " << centroidError << "\n";


    }
 
    
    if (numPeaks == 2) {
        
        double amp_pk1 = fitFunc->GetParameter(1);
        double amp_pk2 = fitFunc->GetParameter(3);
        fitFunc->SetParameter(3,0); // we remove the second peak area.
        double peakArea1 = fitFunc->Integral(fitMin, fitMax) - (fitFunc->GetParameter(8) * (fitMax - fitMin) + fitFunc->GetParameter(9) * 0.5 * (fitMax * fitMax - fitMin * fitMin));
        double peakAreaError1 = fitFunc->IntegralError(fitMin, fitMax, fitFunc->GetParameters(), covMatrix.GetMatrixArray());
        double centroid1 = fitFunc->GetParameter(0);
        double centroidError1 = fitFunc->GetParError(0);

        std::cout << "Peak Area 1: " << peakArea1 << " ± " << peakAreaError1 << "\n";
        std::cout << "Centroid 1: " << centroid1 << " ± " << centroidError1 << "\n";

        fitFunc->SetParameter(1,0);
        fitFunc->SetParameter(3,amp_pk2);
        double peakArea2 = fitFunc->Integral(fitMin, fitMax) - (fitFunc->GetParameter(8) * (fitMax - fitMin) + fitFunc->GetParameter(9) * 0.5 * (fitMax * fitMax - fitMin * fitMin));
        double peakAreaError2 = fitFunc->IntegralError(fitMin, fitMax, fitFunc->GetParameters(), covMatrix.GetMatrixArray());
        double centroid2 = fitFunc->GetParameter(2);
        double centroidError2 = fitFunc->GetParError(2);

        std::cout << "Peak Area 2: " << peakArea2 << " ± " << peakAreaError2 << "\n";
        std::cout << "Centroid 2: " << centroid2 << " ± " << centroidError2 << "\n";
    }






    char fitAnother;
    std::cout << "Do you want to fit another peak? (y/n): ";
    std::cin >> fitAnother;

    if (fitAnother == 'y') {
        FitSkewedGaussianPeak();
    }

    delete c1;
    delete fitFunc;
    file->Close();
    delete file;
}

int main() {
    TApplication app("FitSkewedGaussianPeak", nullptr, nullptr);
    FitSkewedGaussianPeak();
    app.Run();
    return 0;
}
