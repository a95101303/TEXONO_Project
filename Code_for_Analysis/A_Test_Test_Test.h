double Test3(int A)
{
    double B;
    if(A==1)B=2;
    return B;
}

double Test2(double *A)
{
    cout << "TypeA: " << typeid(*A).name() << endl;
    cout << "A: " << *(A+1) << endl;
    cout << "Position0: " << A[0] << endl;
    cout << "Position1: " << A[1] << endl;
    return 0;
}

double Test(double *Position1)
{
    Test2(Position1);
    cout << "Type: " << typeid(Position1).name() << endl;
    return 0;
}

double Test1(vector<double> Position1)
{
    cout << "Position1: " << Position1[0] << endl;
    cout << "Type: " << typeid(Position1).name() << endl;
    return Position1[0];
}

double The_smallest_in_the_vector(vector<double> Vector)
{
    double Return_the_smallest_one=0;
    if(Vector.size()>0)
    {
        std::sort(Vector.begin(), Vector.end());
        Return_the_smallest_one = Vector[0];
    }
    else{Return_the_smallest_one=0;}
    
    return Return_the_smallest_one;
}
