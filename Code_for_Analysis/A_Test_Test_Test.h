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

