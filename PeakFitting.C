#include "PeakFunctions.h"


void FitSkewedGaussianPeak() {
    int numPeaks;
    double minp, maxp;
    double fitMin, fitMax, peakCentroid;

    TFitResultPtr fitResult;
    TH1D *hist = nullptr;

    TCanvas *c1 = new TCanvas("c1", "Fit Skewed Gaussian with Two Tails", 800, 600);

    TFile *file = TFile::Open("/home/sangeetpannu/Phd_110Cd_Project/ROOT_FILES/Completed/TimeSubbedSingles.root");
    if (!file || file->IsZombie()) {
        std::cerr << "Error opening file." << std::endl;
        return;
    }
    
    file->GetObject("LabEinTG", hist);
    if (!hist) {
        std::cerr << "Histogram 'LabEinTG' not found." << std::endl;
        return;
    }

    hist->Draw();
    c1->Modified();
    c1->Update();

    

    while(1)
    {
        std::cout << "Enter domain range for viewing (min max): ";
        std::cin >> minp >> maxp;
        if(minp == maxp)break;
        hist->SetAxisRange(minp,maxp);
         c1->Modified();
         c1->Update();
        
    }
 
    
    std::cout << "Enter fitting domain (min max): ";
    std::cin >> fitMin >> fitMax;
      std::cout << "Enter the number of peaks to fit (1 or 2 or 3(2 peaks with different Sigmas)): ";
        std::cin >> numPeaks;

        TF1 *fitFunc;
        if (numPeaks == 1) {
            fitFunc = new TF1("fitFunc", SkewedGaussianWithBackground1, fitMin, fitMax, 8);
            double peakCentroid;
            std::cout << "Enter peak centroid: ";
            std::cin >> peakCentroid;
            fitFunc->SetParameters(peakCentroid, TMath::Sqrt(9.0 + 4. * peakCentroid / 1000.) / 10., 0, 10000, 1, 1, hist->GetMinimum(), 0);
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
            fitFunc->SetParameters(peakCentroid1, hist->GetMaximum(), peakCentroid2, hist->GetMaximum() / 2,1, 0, 1, 1, hist->GetMinimum(), 0);
            fitFunc->SetParNames("Centroid 1 ","amplitude 1","Centroid 2","amplitude 2","Sigma overall","skewness","left tail par","right tail par","linear Intercept","Linear Slope");
            fitFunc->SetParLimits(0, peakCentroid1-3, peakCentroid1+3); 
            fitFunc->SetParLimits(1, 0, 1e+10);
            fitFunc->SetParLimits(2, peakCentroid2-3, peakCentroid2+3);
            fitFunc->SetParLimits(3, 0, 1e+10);
            fitFunc->SetParLimits(4, 1, 20);

        }  else if (numPeaks == 3) {
            fitFunc = new TF1("fitFunc", SkewedGaussianWithBackground2N2Sigma, fitMin, fitMax, 11);
            double peakCentroid1, peakCentroid2;
            std::cout << "Enter first peak centroid: ";
            std::cin >> peakCentroid1;
            std::cout << "Enter second peak centroid: ";
            std::cin >> peakCentroid2;
            fitFunc->SetParameters(peakCentroid1, hist->GetMaximum()/2.0, peakCentroid2, hist->GetMaximum(),6, 3, 0, 1, 1, hist->GetMinimum(), -4);
            fitFunc->SetParNames("Centroid 1 ","amplitude 1","Centroid 2","amplitude 2","Sigma 1", "Sigma 2","skewness","left tail par","right tail par","linear Intercept","Linear Slope");
            fitFunc->SetParLimits(0, peakCentroid1-3, peakCentroid1+3); 
          //  fitFunc->SetParLimits(1, 0, 1e+10);
            fitFunc->SetParLimits(2, peakCentroid2-3, peakCentroid2+3);
           // fitFunc->SetParLimits(3, 0, 1e+10);
            fitFunc->SetParLimits(8, 0, 1);
            fitFunc->SetParLimits(9, 0, 1);
            fitFunc->SetParLimits(4, 1, 20);
            fitFunc->SetParLimits(5, 1, 20);
     //       fitFunc->SetParLimits(10, -10, 0);

        } else if (numPeaks == 4) {
            fitFunc = new TF1("fitFunc", SkewedGaussianWithBackgroundStep, fitMin, fitMax, 11);
            double peakCentroid1, peakCentroid2;
            std::cout << "Enter first peak centroid: ";
            std::cin >> peakCentroid1;
            std::cout << "Enter second peak centroid: ";
            std::cin >> peakCentroid2;
            fitFunc->SetParameters(peakCentroid1, hist->GetMaximum(), peakCentroid2, hist->GetMaximum() / 2,4, 0, 1, 1, peakCentroid1, -0.8,10000);
            fitFunc->SetParNames("Centroid 1 ","amplitude 1","Centroid 2","amplitude 2","Sigma overall","skewness","left tail par","right tail par","bk-step offset","beta","amplitude");
            fitFunc->SetParLimits(0, peakCentroid1-3, peakCentroid1+3); 
          //  fitFunc->SetParLimits(1, 0, 1e+10);
            fitFunc->SetParLimits(2, peakCentroid2-3, peakCentroid2+3);
           // fitFunc->SetParLimits(3, 0, 1e+10);
            //fitFunc->SetParLimits(8, 0, 1);
            fitFunc->SetParLimits(11, 0, 1e+08);
           // fitFunc->SetParLimits(9, 0, 1);
        //    fitFunc->SetParLimits(4, 1, 8);
            //fitFunc->SetParLimits(5, 1, 20);
        } else{
            std::cout << "Entered to many peaks to fit" << std::endl;
            return;
        }

    
auto fitHistogram = [&](bool userSetLimits = false) {
        if (userSetLimits) {
            std::cout << "(Enter 1 or 2 for limiting parameters or fixing them or -1 to skip), For 1: Enter parameter limits (low high) for each parameter, For 2: Enter fixed parameter.:\n";
            for (int i = 0; i <11; ++i) {
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
        fitResult = hist->Fit(fitFunc, "SRE");


        if (numPeaks == 1)
        {

            hist->SetAxisRange(fitFunc->GetParameter(0)-5,fitFunc->GetParameter(0)+5);

            int bin = fitMax-fitMin;
            double Yrange =  hist->GetMaximum();
            hist->SetAxisRange(0,Yrange+Yrange/3.0,"Y");
            hist->SetAxisRange(fitMin,fitMax);
            TF1 *Fit1 = new TF1("Fit 1",SkewedGaussianPeak,fitMin,fitMax,6);
            Fit1->SetParameters(fitFunc->GetParameter(0),fitFunc->GetParameter(1),fitFunc->GetParameter(2),fitFunc->GetParameter(3),fitFunc->GetParameter(4),fitFunc->GetParameter(5));
            Fit1->SetLineColor(kGreen);
            Fit1->Draw("same");
          
            TF1 *bkfunc = new TF1("background", background, fitMin, fitMax, 2);
            bkfunc->SetParameter(0,fitFunc->GetParameter(7));
            bkfunc->SetParameter(1,fitFunc->GetParameter(6));
            bkfunc->SetLineColor(kRed);
            bkfunc->Draw("same");

            c1->Modified();
            c1->Update();

            TCanvas *c2 = new TCanvas("c2", "Fit (low range view)", 800, 600);
            c2->cd();
            hist->Draw();
            fitResult->Draw("same");
            double ylow = bkfunc->Eval(fitMin)+bkfunc->Eval(fitMin);
            double yhigh = bkfunc->Eval(fitMax)-bkfunc->Eval(fitMax)/4;
            bkfunc->SetLineColor(kRed);
            bkfunc->Draw("same");
            hist->SetAxisRange(fitMin-10,fitMax+10);
            hist->SetAxisRange(yhigh,ylow,"Y");
            c2->Modified();
            c2->Update();
         //   c1->cd();

        }else if(numPeaks == 2){

            hist->SetAxisRange(fitFunc->GetParameter(0)-5,fitFunc->GetParameter(0)+5);
            double Yrange =  hist->GetMaximum();
            hist->SetAxisRange(0,Yrange+Yrange/3.0,"Y");
            hist->SetAxisRange(fitMin-25,fitMax+25);
            TF1 *Fit1 = new TF1("Fit 1",SkewedGaussianPeak,fitMin,fitMax,6);
            Fit1->SetParameters(fitFunc->GetParameter(0),fitFunc->GetParameter(4),fitFunc->GetParameter(5),fitFunc->GetParameter(1),fitFunc->GetParameter(6),fitFunc->GetParameter(7));
            Fit1->SetLineColor(kGreen);
            Fit1->Draw("same");

            TF1 *Fit2 = new TF1("Fit 2",SkewedGaussianPeak,fitMin,fitMax,6);
            Fit2->SetParameters(fitFunc->GetParameter(2),fitFunc->GetParameter(4),fitFunc->GetParameter(5),fitFunc->GetParameter(3),fitFunc->GetParameter(6),fitFunc->GetParameter(7));
            Fit2->SetLineColor(kMagenta);
            Fit2->Draw("same");

            TF1 *bkfunc = new TF1("background", background, fitMin, fitMax, 2);
            bkfunc->SetParameter(0,fitFunc->GetParameter(9));
            bkfunc->SetParameter(1,fitFunc->GetParameter(8));
            bkfunc->SetLineColor(kRed);
            bkfunc->Draw("same");
            c1->Modified();
            c1->Update();

        }else if(numPeaks == 3){

            hist->SetAxisRange(fitFunc->GetParameter(0)-5,fitFunc->GetParameter(0)+5);
            double Yrange =  hist->GetMaximum();
            hist->SetAxisRange(0,Yrange+Yrange/3.0,"Y");
            hist->SetAxisRange(fitMin-25,fitMax+25);

            TF1 *Fit1 = new TF1("Fit 1",SkewedGaussianPeak,fitMin,fitMax,6);
            Fit1->SetParameters(fitFunc->GetParameter(0),fitFunc->GetParameter(4),fitFunc->GetParameter(6),fitFunc->GetParameter(1),fitFunc->GetParameter(7),fitFunc->GetParameter(8));
            Fit1->SetLineColor(kGreen);
            Fit1->Draw("same");

            TF1 *Fit2 = new TF1("Fit 2",SkewedGaussianPeak,fitMin,fitMax,6);
            Fit2->SetParameters(fitFunc->GetParameter(2),fitFunc->GetParameter(5),fitFunc->GetParameter(6),fitFunc->GetParameter(3),fitFunc->GetParameter(7),fitFunc->GetParameter(8));
            Fit2->SetLineColor(kMagenta);
            Fit2->Draw("same");

            TF1 *bkfunc = new TF1("background", background, fitMin, fitMax, 2);
            bkfunc->SetParameter(0,fitFunc->GetParameter(10));
            bkfunc->SetParameter(1,fitFunc->GetParameter(9));
            bkfunc->SetLineColor(kRed);
            bkfunc->Draw("same");
            c1->Modified();
            c1->Update();
      

        } else if(numPeaks == 4){

            hist->SetAxisRange(fitFunc->GetParameter(0)-5,fitFunc->GetParameter(0)+5);
            double Yrange =  hist->GetMaximum();
            hist->SetAxisRange(0,Yrange+Yrange/3.0,"Y");
            hist->SetAxisRange(fitMin-25,fitMax+25);
            TF1 *Fit1 = new TF1("Fit 1",SkewedGaussianPeak,fitMin,fitMax,6);
            Fit1->SetParameters(fitFunc->GetParameter(0),fitFunc->GetParameter(4),fitFunc->GetParameter(5),fitFunc->GetParameter(1),fitFunc->GetParameter(6),fitFunc->GetParameter(7));
            Fit1->SetLineColor(kGreen);
            Fit1->Draw("same");

            TF1 *Fit2 = new TF1("Fit 2",SkewedGaussianPeak,fitMin,fitMax,6);
            Fit2->SetParameters(fitFunc->GetParameter(2),fitFunc->GetParameter(4),fitFunc->GetParameter(5),fitFunc->GetParameter(3),fitFunc->GetParameter(6),fitFunc->GetParameter(7));
            Fit2->SetLineColor(kMagenta);
            Fit2->Draw("same");

            TF1 *bkfunc = new TF1("background", backgroundWStep, fitMin, fitMax, 3);
            bkfunc->SetParameter(0,fitFunc->GetParameter(8));
            bkfunc->SetParameter(1,fitFunc->GetParameter(9));
            bkfunc->SetParameter(2,fitFunc->GetParameter(10));
           // bkfunc->SetParameter(3,fitFunc->GetParameter(11));
            bkfunc->SetLineColor(kRed);
            bkfunc->Draw("same");
            c1->Modified();
            c1->Update();

        }


        


        
    };

    fitHistogram();

    

     // Retrieve chi-squared and ndf
    double chi2 = fitResult->Chi2();
    int ndf = fitResult->Ndf();
    double chi2red = chi2 / ndf;

    // Print the results
    std::cout << "Chi2: " << chi2 << std::endl;
    std::cout << "NDF: " << ndf << std::endl;
    std::cout << "Reduced Chi2: " << chi2red << std::endl;

 
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
        fitFunc->SetParameter(6,0);
        fitFunc->SetParameter(7,0);
        double peakArea = fitFunc->Integral(fitMin, fitMax);
        double peakAreaError = fitFunc->IntegralError(fitMin, fitMax, fitFunc->GetParameters(), covMatrix.GetMatrixArray(),1.0E-6);
        double centroid = fitFunc->GetParameter(0);
        double centroidError = fitFunc->GetParError(0);
        std::cout << "Peak Area: " << peakArea << " ± " << peakAreaError << "\n";
        std::cout << "Centroid: " << centroid << " ± " << centroidError << "\n";


    }
 
    
    if (numPeaks == 2) {
        
        double amp_pk1 = fitFunc->GetParameter(1);
        double amp_pk2 = fitFunc->GetParameter(3);
        fitFunc->SetParameter(3,0); // we remove the second peak area.
        fitFunc->SetParameter(8,0);
        fitFunc->SetParameter(9,0);
        double peakArea1 = fitFunc->Integral(fitMin, fitMax);
        double peakAreaError1 = fitFunc->IntegralError(fitMin, fitMax, fitFunc->GetParameters(), covMatrix.GetMatrixArray(),1.0E-6);
        double centroid1 = fitFunc->GetParameter(0);
        double centroidError1 = fitFunc->GetParError(0);

        std::cout << "Peak Area 1: " << peakArea1 << " ± " << peakAreaError1 << "\n";
        std::cout << "Centroid 1: " << centroid1 << " ± " << centroidError1 << "\n";

        fitFunc->SetParameter(1,0);
        fitFunc->SetParameter(3,amp_pk2);
        double peakArea2 = fitFunc->Integral(fitMin, fitMax);
        double peakAreaError2 = fitFunc->IntegralError(fitMin, fitMax, fitFunc->GetParameters(), covMatrix.GetMatrixArray(),1.0E-6);
        double centroid2 = fitFunc->GetParameter(2);
        double centroidError2 = fitFunc->GetParError(2);

        std::cout << "Peak Area 2: " << peakArea2 << " ± " << peakAreaError2 << "\n";
        std::cout << "Centroid 2: " << centroid2 << " ± " << centroidError2 << "\n";
    }

    if (numPeaks == 3) {
        
        double amp_pk1 = fitFunc->GetParameter(1);
        double amp_pk2 = fitFunc->GetParameter(3);
        fitFunc->SetParameter(10,0);
        fitFunc->SetParameter(9,0);
        fitFunc->SetParameter(3,0); // we remove the second peak area.
        double peakArea1 = fitFunc->Integral(fitMin, fitMax);
        double peakAreaError1 = fitFunc->IntegralError(fitMin, fitMax, fitFunc->GetParameters(), covMatrix.GetMatrixArray(),1.0E-6);
        double centroid1 = fitFunc->GetParameter(0);
        double centroidError1 = fitFunc->GetParError(0);

        std::cout << "Peak Area 1: " << peakArea1 << " ± " << peakAreaError1 << "\n";
        std::cout << "Centroid 1: " << centroid1 << " ± " << centroidError1 << "\n";

        fitFunc->SetParameter(1,0);
        fitFunc->SetParameter(3,amp_pk2);
        double peakArea2 = fitFunc->Integral(fitMin, fitMax);
        double peakAreaError2 = fitFunc->IntegralError(fitMin, fitMax, fitFunc->GetParameters(), covMatrix.GetMatrixArray(),1.0E-6);
        double centroid2 = fitFunc->GetParameter(2);
        double centroidError2 = fitFunc->GetParError(2);

        std::cout << "Peak Area 2: " << peakArea2 << " ± " << peakAreaError2 << "\n";
        std::cout << "Centroid 2: " << centroid2 << " ± " << centroidError2 << "\n";
    }

    if (numPeaks == 4) {
        
        double amp_pk1 = fitFunc->GetParameter(1);
        double amp_pk2 = fitFunc->GetParameter(3);
        fitFunc->SetParameter(3,0); // we remove the second peak area.
        fitFunc->SetParameter(8,0);
        fitFunc->SetParameter(9,0);
        fitFunc->SetParameter(10,0);
        fitFunc->SetParameter(11,0);
        double peakArea1 = fitFunc->Integral(fitMin, fitMax);
        double peakAreaError1 = fitFunc->IntegralError(fitMin, fitMax, fitFunc->GetParameters(), covMatrix.GetMatrixArray());
        double centroid1 = fitFunc->GetParameter(0);
        double centroidError1 = fitFunc->GetParError(0);

        std::cout << "Peak Area 1: " << peakArea1 << " ± " << peakAreaError1 << "\n";
        std::cout << "Centroid 1: " << centroid1 << " ± " << centroidError1 << "\n";

        fitFunc->SetParameter(1,0);
        fitFunc->SetParameter(3,amp_pk2);
        double peakArea2 = fitFunc->Integral(fitMin, fitMax);
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

    if(numPeaks==1) c1->SaveAs(Form("Fit_Peak1_%f.root",fitFunc->GetParameter(0)));

    if(numPeaks==2 || numPeaks == 4) c1->SaveAs(Form("Fit_Peak1_%f_Peak2_%f.root",fitFunc->GetParameter(0),fitFunc->GetParameter(2)));

    if(numPeaks==3) c1->SaveAs(Form("Fit_Peak1_%f_Peak2_%f_DiffSigmas.root",fitFunc->GetParameter(0),fitFunc->GetParameter(2)));

    delete c1;
   // delete c2;
    delete fitFunc;
    file->Close();
    delete file;
}

