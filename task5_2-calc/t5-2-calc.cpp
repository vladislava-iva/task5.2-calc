#include <iostream>
#include <cmath>
#include <string>
#include <iomanip>
using namespace std;

const int n = 6;
double xi[] = {5, 15, 25, 35, 45, 55};
double yi[] = {2.2, 2.4, 2.6, 2.7, 2.8, 2.9};

void kram2(double a00, double a01, double a10, double a11,double b0, double b1, double& x, double& y)
{
    double D = a00*a11 - a01*a10;
    double D1 = b0*a11 - a01*b1;
    double D2 = a00*b1 - b0 *a10;
    cout << "de="<< D <<"\n";
    cout << "de1="<<D1<<"\n";
    cout << "de2=" <<D2<<"\n";
    x = D1 / D;
    y = D2 / D;
}

void gauss3(double A[3][4],double& ra,double& rb, double& rc)
{
    cout << "начальная расширенная матрица:\n";
    for(int r=0;r<3;r++){
        for(int j=0;j<=3;j++) cout <<A[r][j]<<" ";
        cout << endl;
    }
    for(int col = 0; col < 3; col++){
        int piv = col;

        for(int row =col+1; row < 3; row++){
            if(fabs(A[row][col]) > fabs(A[piv][col])) piv = row;
        }

        for(int j = 0; j <= 3; j++) swap(A[col][j], A[piv][j]);
        double d = A[col][col];
        for(int j = col; j <= 3; j++) A[col][j] /= d;

        for(int row = 0; row < 3; row++){
            if(row == col) continue;
            double cf = A[row][col];
            for(int j = col; j <= 3; j++) A[row][j] -= cf * A[col][j];
        }
    }

    cout << "матрица после диагонализации:\n";
    for(int r=0;r<3;r++){
        for(int j=0;j<=3;j++) cout << A[r][j]<<" ";
        cout << endl;
    }
    ra = A[0][3];
    rb = A[1][3];
    rc = A[2][3];
}

double yApprox(const string& tip, double a, double b, double c, double x)
{
    if(tip == "lin") return a*x + b;
    if(tip == "step") return b * pow(x, a);
    if(tip == "pok") return b * exp(a * x);
    return a*x*x + b*x + c;
}

void vivTabl(const string& tip, double a, double b, double c = 0)
{
    cout << "\nтаблица значений:\n";
    cout << "xi yi y_approx (dy)^2"<< "\n";
    for(int i = 0; i < n; i++){
        double yf = yApprox(tip, a, b, c, xi[i]);
        double dy = yf - yi[i];
        cout << xi[i] << " "<<yi[i]<<" " << fixed << setprecision(4) << yf <<" "<< fixed << setprecision(4) << dy*dy << "\n";
    }
}

double pogr(const string& tip, double a, double b, double c = 0)
{
    double s = 0;
    for(int i = 0; i < n; i++){
        double dy = yApprox(tip, a, b, c, xi[i]) - yi[i];
        s += dy*dy;
    }
    return s;
}

int main()
{
    cout << fixed << setprecision(4);

    double sx=0, sy=0, sx2=0, sxy=0, sx3=0, sx4=0, sx2y=0;
    for(int i = 0; i < n; i++){
        sx += xi[i];
        sy += yi[i];
        sx2 += xi[i]*xi[i];
        sxy += xi[i]*yi[i];
        sx3 += xi[i]*xi[i]*xi[i];
        sx4 += xi[i]*xi[i]*xi[i]*xi[i];
        sx2y += xi[i]*xi[i]*yi[i];
    }

    cout << "\nвар 14: xi={5,15,25,35,45,55}, yi={2.2,2.4,2.6,2.7,2.8,2.9}\n";
    cout << "\nпредварительные суммы:\n";
    cout << "сумма xi=" << sx << "\n";
    cout << "сумма yi=" << sy << "\n";
    cout << "сумма xi^2=" << sx2 << "\n";
    cout << "сумма xi*yi=" << sxy << "\n";
    cout << "сумма xi^3=" << sx3 << "\n";
    cout << "сумма xi^4=" << sx4 << "\n";
    cout << "сумма xi^2*yi=" << sx2y << "\n";

    double pS[4];

    cout << "\n\n1.1. y=a*x+b\n";
    cout << "\nсистема нормальных уравнений МНК:\n";

    cout  << sx2 << "*a+" << sx << "*b=" << sxy << "\n";
    cout  << sx << "*a+" << n << "*b=" << sy << "\n";

    cout << "\nрешение по формулам Крамера:\n";
    double a1, b1;
    kram2(sx2, sx, sx, n, sxy, sy, a1, b1);
    cout << "\na=de1/de=" << a1 << " (округл. " << round(a1*100)/100 << ")\n";
    cout << "b=de2/de=" << b1 << " (округл. " << round(b1*100)/100 << ")\n";
    cout << "\nрез-тат: y=" << a1 << "*x+(" << b1 << ")\n";
    cout << "округлённо: y=" << round(a1*100)/100 << "*x+" << round(b1*100)/100 << "\n";

    vivTabl("lin", a1, b1);
    pS[0] = pogr("lin", a1, b1);
    cout << "\nсуммарная погрешность S1=" << pS[0] << "\n";

    cout << "\n\n1.2. y=be*x^a\n";

    cout << "\nлогарифмируем: ln(y)=a*ln(x)+ln(be)\n";
    cout << "замена: Y=ln(y), X=ln(x), b=ln(be)\n";
    double sX=0, sY=0, sXX=0, sXY=0;
    cout << "\nпреобразованная таблица:\n";
    cout << "X=ln(xi)" << " " << "Y=ln(yi)" << " " << "X^2" << " " << "X*Y" << "\n";

    for(int i = 0; i < n; i++){
        double X = log(xi[i]), Y = log(yi[i]);
        cout << X << " " << Y << " " << X*X << " " << X*Y << "\n";
        sX += X; sY += Y; sXX += X*X; sXY += X*Y;
    }

    cout << "сумма X=" << sX << "сумма Y=" << sY << "cумма X^2=" << sXX << "сумма XY=" << sXY << "\n";
         
    cout << "\ncистема:\n";
    cout << sXX << "*a+" << sX << "*b=" << sXY << "\n";
    cout << sX << "*a+" << n << "*b=" << sY << "\n";
    cout << "\nрешение по формулам Крамера:\n";

    double a2, b2;
    kram2(sXX, sX, sX, n, sXY, sY, a2, b2);
    double be2 = exp(b2);
    cout << "\na=" << a2 << " (округл. " << round(a2*100)/100 << ")\n";
    cout << "b=" << b2 << " след. be=e^b=" << be2 << " (округл. " << round(be2*100)/100 << ")\n";
    cout << "\nрез-тат: y=" << be2 << "*x^" << a2 << "\n";
    cout << "округлённо: y=" << round(be2*100)/100 << "*x^" << round(a2*100)/100 << "\n";

    vivTabl("step", a2, be2);
    pS[1] = pogr("step", a2, be2);
    cout << "\nсуммарная погрешность S2=" << pS[1] << "\n";

    cout << "\n\n1.3. y=be*e^(a*x)\n";
    cout << "\nлогарифмируем: ln(y)=a*x+ln(be)\n";
    cout << "замена: Y=ln(y), b=ln(be)\n";

    double sY3=0, sXY3=0;
    cout << "\nпреобразованная таблица:\n";
    cout << "xi Y=ln(yi) xi*Y" << "\n";
    for(int i = 0; i < n; i++){
        double Y = log(yi[i]);
        cout << xi[i]<< " " << Y  << " " << xi[i]*Y << "\n";

        sY3 += Y;
        sXY3 += xi[i]*Y;
    }
    
    cout << "сумма Y=" << sY3 << " сумма x*Y=" << sXY3 << "\n";
    cout << "(сумма xi=" << sx << ", сумма xi^2=" << sx2 << ")\n";
    cout << "\ncистема:\n";
    cout << sx2 << "*a+" << sx << "*b=" << sXY3 << "\n";
    cout  << sx << "*a+" << n << "*b=" << sY3 << "\n";

    cout << "\nрешение по формулам Крамера:\n";
    double a3, b3;
    kram2(sx2, sx, sx, n, sXY3, sY3, a3, b3);
    double be3 = exp(b3);
    cout << "\na=" << a3 << " (округл. " << round(a3*10000)/10000 << ")\n";
    cout << "b=" << b3 << " след. be=e^b=" << be3 << " (округл. " << round(be3*100)/100 << ")\n";
    cout << "\nрез-тат: y=" << be3 << " * e^(" << a3 << "*x)\n";
    cout << "округлённо: y=" << round(be3*100)/100<< " * e^(" << round(a3*10000)/10000 << "*x)\n";

    vivTabl("pok", a3, be3);
    pS[2] = pogr("pok", a3, be3);
    cout << "\nсуммарная погрешность S3=" << pS[2] << "\n";

    cout << "\n\n1.4. y=a*x^2+b*x+c\n";
    cout << "\nсистема нормальных уравнений:\n";
    cout << sx4 << "*a+" << sx3 << "*b+" << sx2 << "*c=" << sx2y << "\n";
    cout << sx3 << "*a+" << sx2 << "*b+" << sx << "*c=" << sxy << "\n";
    cout << sx2 << "*a+" << sx << "*b+" << n << "*c=" << sy << "\n";

    cout << "\nрешение методом Гаусса-Жордана:\n";
    double M[3][4] = {
        {sx4, sx3, sx2, sx2y},
        {sx3, sx2, sx, sxy },
        {sx2, sx, (double)n, sy}
    };
    double a4, b4, c4;
    gauss3(M, a4, b4, c4);
    cout << "\na=" << a4 << " (округл. " << round(a4*100)/100 << ")\n";
    cout << "b=" << b4 << " (округл. " << round(b4*100)/100 << ")\n";
    cout << "c=" << c4 << " (округл. " << round(c4*100)/100 << ")\n";
    cout << "\nрезультат: y=" << a4 << "*x^2+("<< b4 << ")*x + (" << c4 << ")\n";
    cout << "округлённо: y=" << round(a4*100)/100<< "*x^2+" << round(b4*100)/100<< "*x+" << round(c4*100)/100 << "\n";
    vivTabl("kv", a4, b4, c4);
    pS[3] = pogr("kv", a4, b4, c4);
    cout << "\nсуммарная погрешность S4 = " << pS[3] << "\n";

    cout << "\n\n3. сравнение погрешностей\n\n";
    string names[4] = {
        "линейная y=a*x+b",
        "степенная y=be*x^a",
        "показательная y=be*e^(ax)",
        "квадратичная. y=a*x^2+b*x+c"
    };
    int best = 0;
    for(int i = 1; i < 4; i++) if(pS[i] < pS[best]) best = i;
    for(int i = 0; i < 4; i++){
        cout << "S" << i+1 << "=" << " " << pS[i] << " " << names[i];
        if(i == best) cout << " min";
        cout << "\n";
    }
    cout << "\nвывод: наилучшая аппроксимация - " << names[best] << "\n\n";
    return 0;
}