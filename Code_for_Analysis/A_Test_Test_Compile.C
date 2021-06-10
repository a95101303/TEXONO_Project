#include "A_Test_Test_Test.h"

void A_Test_Test_Compile()
{
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

    return 0;

}
