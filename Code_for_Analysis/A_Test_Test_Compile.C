#include "A_Test_Test_Test.h"
#include "TH1.h"
#include "velocity_distribution_2000_Ave.h"

void A_Test_Test_Compile()
{
    
    double Try_array[3];
    for(int KKK=1; KKK<3; KKK++)
    {
        Try_array[KKK]=2;
    }
    
    for(int KKK=0; KKK<3; KKK++)
    {
        cout << "Try_array[KKK]: " << Try_array[KKK] << endl;
    }

    TH1F   *Test_Integral = new TH1F("Test_Integral","Test_Integral",40,0,10);
    Test_Integral->Fill(5);Test_Integral->Fill(9.1);
    Test_Integral->Fill(9.5);Test_Integral->Fill(7.8);
    cout << "Integral(): " << Test_Integral->Integral(1,40) << endl;
    
    double Vecolity[2000];double Possiblity[2000];
    double sum; for(int j=0;j<2000;j++){sum = sum + velo_dist_Ave[j][3];}
    for(int j=0;j<2000;j++)
    {
        float v = 0.5*(velo_dist_Ave[j][1]+velo_dist_Ave[j][2]);
        Vecolity[j] = v;
        Possiblity[j] = velo_dist_Ave[j][3]/sum;
        //cout << "Possiblity[j]: " << Possiblity[j] << endl;
    }
    TH1F   *Flux_HIST = new TH1F("Flux_HIST","Flux_HIST",2000,0,791);
    for(int kkk=0;kkk<2000;kkk++){Flux_HIST->SetBinContent(kkk+1,Possiblity[kkk]);}

    cout << "GetBinContent(): " << Flux_HIST->GetBinCenter(Flux_HIST->GetMaximumBin()) << endl;
    Int_t binx        = Flux_HIST->GetXaxis()->FindBin(800);
    Int_t binx_Number = Flux_HIST->GetNbinsX();
    double Bin_content = Flux_HIST->GetBinContent(binx);
    
    cout << "binx: " << binx << endl;
    cout << "Bin_content: " << Bin_content  << endl;
    cout << "binx_Number: " << binx_Number << endl;
    int myarray[6]{10, 4, 14, 84, 1, 3};

    std::cout << std::find(std::begin(myarray), std::end(myarray), 1) << endl;
    std::cout << std::find(std::begin(myarray), std::end(myarray), 3) << endl;
    std::cout << std::end(myarray) << endl;
    std::cout << std::find(std::begin(myarray), std::end(myarray), 5) << endl;
    std::cout << std::find(std::begin(myarray), std::end(myarray), 7) << endl;

    if (std::find(std::begin(myarray), std::end(myarray), 5) == std::end(myarray))
    {
        std::cout << std::find(std::begin(myarray), std::end(myarray), 5) << endl;
        std::cout << "It doesn't exists5";
    }
    if (std::find(std::begin(myarray), std::end(myarray), 7) != std::end(myarray))
        std::cout << "It exists7";

    
    double Position1[2]={1,2};
    vector<double> Position2={10,5,8};
    cout << "Test3: " << Test3(1) << endl;
    cout << "The_smallest_in_the_vector(Position2): " << The_smallest_in_the_vector(Position2) << endl;
    
    gRandom = new TRandom3(0);
    gRandom->SetSeed(0);

    TF1 *Two_Sides    = new TF1("Two_Sides","x",0,2);
    double Side_confirmed = Two_Sides->GetRandom();//Give out 0,1;
    
    cout << "Side_confirmed: " << Side_confirmed << endl;
    
    int A=0;
    int B=0;
    if(A==0) B=1;
    if(B==1) cout << "ok: " << endl;

    double Try[3]={1,2,3};
    double DR[3]={1,2,3};
    double *Try_Point=Try;
    
    cout << "Try_Point[0]: " << Try_Point[0] << endl;
    cout << "Try_Point[1]: " << Try_Point[1] << endl;
    cout << "Try_Point[2]: " << Try_Point[2] << endl;

    Try_Point = PAP(1,Try_Point,DR);
    
    cout << "Try_Point[0]: " << Try_Point[0] << endl;
    cout << "Try_Point[1]: " << Try_Point[1] << endl;
    cout << "Try_Point[2]: " << Try_Point[2] << endl;

    double BBB[2]={1,2};
    TF1 *fa = new TF1("fa","GGG(x)",0,100000);
    cout << "fa->GetX(1,0.1,20):" << fa->GetX(770,0,100000) << endl;
    
    int C=0;int D=0;
    if(C<2) D=2;
    if(C>-1) D=3;
    
    cout << "D: " << D << endl;
    return 0;

}

