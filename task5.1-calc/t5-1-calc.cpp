#include <iostream>
#include <cmath>
#include <vector>
#include <string>

using namespace std;

const int n = 6;
double xi[] = {5, 15, 25, 35, 45, 55};
double yi[] = {2.2, 2.4, 2.6, 2.7, 2.8, 2.9};

void kram2(double a00, double a01, double a10, double a11, double b0, double b1, double& x, double& y) {
    double det = a00*a11-a01*a10;
    double det1 = b0 * a11 -a01*b1;
    double det2 = a00*b1 -b0 * a10;
    cout << "de = " <<det<<"\n";
    cout << "de1 = " <<det1<<"\n";
    cout << "de2 = " <<det2<<"\n";
    x=det1/det;
    y=det2/det;
}

void gauss3(double A[3][4], double& ra, double& rb, double& rc) {
    for (int col = 0; col < 3; col++) {
        int piv = col;
        for (int row = col+1; row <3; row++){
            if (fabs(A[row][col]) >fabs(A[piv][col])) piv = row;
        }
        for (int j = 0; j <= 3; j++) swap(A[col][j], A[piv][j]);
        double di = A[col][col];
        for (int j = col; j <= 3; j++) A[col][j] /= di;
        for (int row = 0; row < 3; row++) {
            if (row == col) continue;
            double cf = A[row][col];
            for (int j = col; j <= 3; j++) A[row][j] -= cf * A[col][j];
        }
    }
    ra = A[0][3]; rb = A[1][3]; rc = A[2][3];
}

void vivTabl(const string& tip, double a, double b, double c = 0) {
    cout << "\nтаблица значений:\n";
    cout << "xi yi y_approx (dy)^2\n";
    
    for (int i = 0; i < n; i++) {
        double yf;
        if (tip == "lin") yf = a*xi[i] + b;
        else if (tip == "step") yf = b * pow(xi[i], a);
        else if (tip == "pok") yf = b * exp(a * xi[i]);
        else yf = a*xi[i]*xi[i] + b*xi[i] + c;
        
        double dy = yf - yi[i];
        
        cout  << xi[i] << " " << yi[i] << " " << fixed  << yf << " " << dy*dy << "\n";
    }
}

double pogr(const string& tip, double a, double b, double c = 0) {
    double s = 0;
    for (int i = 0; i < n; i++) {
        double yf;
        if (tip == "lin") yf = a*xi[i] + b;
        else if (tip == "step") yf = b * pow(xi[i], a);
        else if (tip == "pok") yf = b * exp(a * xi[i]);
        else yf = a*xi[i]*xi[i] + b*xi[i] + c;
        
        double dy = yf - yi[i];
        s += dy * dy;
    }
    return s;
}

int main() {
    cout << fixed;
    cout << "\n";
  

    double s_x=0, s_y=0, s_xx=0, s_xy=0, s_x3=0, s_x4=0, s_x2y=0;
    for (int i = 0; i < n; i++) {
        s_x += xi[i];
        s_y += yi[i];
        s_xx += xi[i]*xi[i];
        s_xy += xi[i]*yi[i];
        s_x3 += xi[i]*xi[i]*xi[i];
        s_x4 += xi[i]*xi[i]*xi[i]*xi[i];
        s_x2y += xi[i]*xi[i]*yi[i];
    }

    double pS[4];

    cout<< "1.1 y = a*x + b"<<endl;
    cout << "\nиз условия dS/da=0, dS/db=0 получаем систему:\n";
    cout << "сумма xi^2*a + сумма xi * b = сумма xi*yi\n";
    cout << "сумма xi*a+ n * b = сумма yi\n\n";
    cout  << s_xx << "*a+" << s_x << "*b= " << s_xy << "\n";
    cout  << s_x << "*a+" << n << "*b= " << s_y << "\n\n";
    cout << "решение по формулам Крамера:\n";

    double a1, b1;
    kram2(s_xx, s_x, s_x, n, s_xy, s_y, a1, b1);
    cout << "\na=de1/de=" << a1 << "\n";
    cout << "b=de2/de=" << b1 << "\n";
    cout << "\nрезультат: y=" << a1 << "*x + (" << b1 << ")\n";
    cout << "округлённо: y=" << round(a1*100)/100 << "*x+(" << round(b1*100)/100 << ")\n";
    
    vivTabl("lin", a1, b1);
    pS[0] = pogr("lin", a1, b1);
    cout << "суммарная погрешность S1=" << pS[0] << "\n";

    cout<<"\n1.2 y = be * x^a"<<endl;
    cout << "\nлогарифмируем обе части: ln(y)=a*ln(x)+ln(be)\n";
    cout << "обозначим: Y=ln(y), X=ln(x), b=ln(be)\n";
    cout << "тогда задача сводится к линейной: Y=a*X+b\n\n";

    double sX=0, sY=0, sXX=0, sXY=0;
    cout << "преобразованная таблица X=ln(xi), Y=ln(yi):\n";
    cout << "X=ln(xi) Y=ln(yi)\n";
    for (int i = 0; i < n; i++) {
        double X = log(xi[i]), Y = log(yi[i]);
        cout << X << " " << Y << "\n";
        sX += X; sY += Y; sXX += X*X; sXY += X*Y;
    }
    cout << "\nсумма X= " << sX << "\n";
    cout << "сумма Y= " << sY << "\n";
    cout << "сумма X^2= " << sXX << "\n";
    cout << "сумма XY= " << sXY << "\n\n";

    cout << "система:\n";
    cout  << sXX << "*a +" << sX << "*b=" << sXY << "\n";
    cout  << sX << "*a+" << n << "*b = " << sY << "\n\n";
    cout << "решение по формулам Крамера:\n";
    double a2, b2;
    kram2(sXX, sX, sX, n, sXY, sY, a2, b2);
    double be2 = exp(b2);
    cout << "\na = " << a2 << "\n";
    cout << "b = " << b2 << " след-но be = e^b= " << be2 << "\n";
    cout << "\nрезультат: y=" << be2 << "*x^" << a2 << "\n";
    cout << "округлённо: y=" << round(be2*100)/100 << "*x^" << round(a2*100)/100 << "\n";
    
    vivTabl("step", a2, be2);
    pS[1] = pogr("step", a2, be2);
    cout << "суммарная погрешность S2 = " << pS[1] << "\n";

    cout<<"\n1.3 y = be*e^(a*x)"<<endl;
    cout << "\nлогарифмируем обе части: ln(y) = a*x+ln(be)\n";
    cout << "обозначим: Y=ln(y), b=ln(be)\n";
    cout << "тогда задача сводится к линейной: Y=a*x+b\n\n";

    double sX3=0, sY3=0, sXX3=0, sXY3=0;
    cout << "преобразованная таблица x, Y=ln(yi):\n";
    cout << "xi Y=ln(yi)\n";
    for (int i = 0; i < n; i++) {
        double Y = log(yi[i]);
        cout << xi[i] << " " << Y << "\n";
        sX3 += xi[i]; sY3 += Y;
        sXX3 += xi[i]*xi[i]; sXY3 += xi[i]*Y;
    }
    cout << "\nсумма x=" << sX3 << "\n";
    cout << "сумма Y=" << sY3 << "\n";
    cout << "сумма x^2=" << sXX3 << "\n";
    cout << "сумма xY=" << sXY3 << "\n\n";

    cout << "система:\n";
    cout << sXX3 << "*a+" << sX3 << "*b=" << sXY3 << "\n";
    cout << sX3 << "*a+" << n << "*b=" << sY3 << "\n\n";
    cout << "решение по формулам Крамера:\n";
    double a3, b3;

    kram2(sXX3, sX3, sX3, n, sXY3, sY3, a3, b3);
    double be3 = exp(b3);
    cout << "\na=" << a3 << "\n";
    cout << "b=" << b3 << " след-но be=e^b=" << be3 << "\n";
    cout << "\nрезультат: y=" << be3 << " *e^(" << a3 << "*x)\n";
    cout << "округлённо: y=" << round(be3*100)/100 << "*e^(" << round(a3*100)/100 << "*x)\n";
    
    vivTabl("pok", a3, be3);
    pS[2] = pogr("pok", a3, be3);
    cout << "суммарная погрешность S3 = " << pS[2] << "\n";

    cout<<"\n1.4 y=a*x^2 + b*x + c"<<endl;
    cout << "\nиз условия dS/da=0, dS/db=0, dS/dc=0 получаем систему 3x3:\n";
    cout << "сумма x^4*a + сумма x^3*b+сумма x^2*c =сумма x^2*yi\n";
    cout << "сумма x^3*a + сумма x^2*b+сумма x*c =сумма x*yi\n";
    cout << "сумма x^2*a + сумма x *b +n*c=сумма yi\n\n";
    cout << s_x4 << "*a + " << s_x3 << "*b + " << s_xx << "*c = " << s_x2y << "\n";
    cout << s_x3 << "*a + " << s_xx << "*b + " << s_x << "*c = " << s_xy << "\n";
    cout << s_xx << "*a + " << s_x << "*b + " << n << "*c = " << s_y << "\n\n";
    cout << "решение методом Гаусса-Жордана:\n";

    double M[3][4] = {
        {s_x4, s_x3, s_xx, s_x2y},
        {s_x3, s_xx, s_x, s_xy },
        {s_xx, s_x, (double)n, s_y }
    };
    double a4, b4, c4;
    gauss3(M, a4, b4, c4);
    cout << "\na = " << a4 << "\n";
    cout << "b = " << b4 << "\n";
    cout << "c = " << c4 << "\n";
    cout << "\nрезультат: y=" << a4 << "*x^2 + (" << b4 << ")*x+(" << c4 << ")\n";
    cout << "округлённо: y=" << round(a4*100)/100
         << "*x^2 + " << round(b4*100)/100
         << "*x + " << round(c4*100)/100 << "\n";
    
    vivTabl("kv", a4, b4, c4);
    pS[3] = pogr("kv", a4, b4, c4);
    cout << "суммарная погрешность S4=" << pS[3] << "\n";

    cout<<"\n3. сравнение погрешностей"<<endl;
    string names[4] = {"линейная y=a*x+b","степенная y=be*x^a","показательная y = be*e^(ax)","квадратичная y = a*x^2+b*x+c"};
    int best = 0;
    for (int i = 1; i < 4; i++) if (pS[i] < pS[best]) best = i;

    for (int i = 0; i < 4; i++) {
        cout << "S" << i+1 << " = " << pS[i] << " " << names[i] << (i == best ? " - мин" : "") << "\n";
    }
    cout << "\nвывод: наилучшая аппроксимация " << names[best] << "\n\n";

    return 0;
}
