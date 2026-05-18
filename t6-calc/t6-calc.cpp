#define _USE_MATH_DEFINES
#include <iostream>
#include <fstream>
#include <cmath>
#include <cstdlib>
#include <ctime>
#include <vector>
#include <string>
#include <iomanip>
using namespace std;

const int n=14; //вар.

double randAB(double a,double b) {
    return a + (b-a)*(double)rand()/RAND_MAX;
}

void runGnuplot(const string& scr) {
    system(("gnuplot " + scr).c_str());
}

void zad1() {//площадь треугольника
    cout << "\nзадание 1\n";

    double ax,bx,ay,by,S_exact;

    auto f_kus=[](double x) -> double {
        if (x<=n) return 10.0*x/n;
        else return 10.0*(x-20.0)/(n-20.0);
    };

    auto f1_plot=[](double x){ //верхняя
        return 10.0*x/n; 
    };

    auto f2_plot=[](double x){ //нижняя граница
        return 10.0*(x-20.0)/(n-20.0); 
    };

    ax=0; bx=20; ay=0; by=10;
    double S_rect=(bx-ax)*(by-ay);
    S_exact=0.5*bx*by;

    cout << fixed << setprecision(4);
    cout << "S=" << S_exact << "\n\n";
    cout << "N M S_MC |err| otn%"<< "\n";

    vector<int> Ns={100,1000,10000,100000};

    for (int N : Ns) { //лежит ли y между f1(x) и f2(x)
        int M=0;//кол-во попаданий точек
        vector<double> xi_in,yi_in,xi_out,yi_out;

        for (int i=0; i < N; i++) {
            double x=randAB(ax,bx);
            double y=randAB(ay,by);
            double f1 = f1_plot(x);
            double f2 = f2_plot(x);

            if ( (y >= f1 && y <= f2) || (y >= f2 && y <= f1) ) {
                M++;
                xi_in.push_back(x);  
                yi_in.push_back(y);
            } else {
                xi_out.push_back(x); 
                yi_out.push_back(y);
            }
        }

        double S_mc=(double)M/N*S_rect;//тк M/N примерно равно S/Sпрямоуг
        double err_abs=fabs(S_mc-S_exact); //абсолютная погрешность = |S − Sточная|
        double err_otn=err_abs/S_exact*100.0;//относительная = абсолютная / Sточная *100%

        cout << N<< " " << M<< " " << S_mc<< " " << err_abs<< " "  << err_otn << "%\n";

        if (N == 10000) {
            ofstream fin("z1_in.dat"),fout("z1_out.dat"),fcrv("z1_crv.dat");
            for (size_t k=0; k < xi_in.size();  k++) fin  << xi_in[k]  << " " << yi_in[k]  << "\n";
            for (size_t k=0; k < xi_out.size(); k++) fout << xi_out[k] << " " << yi_out[k] << "\n";

            for (int k=0; k<=200; k++) {
                double xv=ax + (bx-ax)*k/200.0;
                fcrv << xv << " " << f1_plot(xv) << " " << f2_plot(xv) << "\n";
            }

            ofstream fvx("z1_vx.dat");
            fvx << n<< " " << 10.0 << "\n";
            fvx << 0.0 << " " << 0.0  << "\n";
            fvx << 20.0<< " "<< 0.0  << "\n";

            ofstream gp("z1.gnu");
            gp << "set terminal pngcairo size 900,600 enhanced font 'Arial,11'\n";
            gp << "set output 'z1.png'\n";
            gp << "set title 'Задание 1: Монте-Карло,треугольник n=" << n
               << "\\nS_MC=" << fixed << setprecision(2) << S_mc
               << " (точно=100),N=10000' font 'Arial,12'\n";
            gp << "set xlabel 'x'\nset ylabel 'y'\n";
            gp << "set grid lw 1 lc rgb '#cccccc'\n";
            gp << "set xrange [-0.5:21]\nset yrange [-0.5:11.5]\n";
            gp << "set key top right\n";
            
            gp << "plot 'z1_out.dat' using 1:2 with points pt 1 ps 0.3"<< " lc rgb '#BBBBBB' title 'вне треугольника',\\\n";
            gp << "     'z1_in.dat'  using 1:2 with points pt 1 ps 0.3"<< " lc rgb '#CC2200' title 'внутри',\\\n";
            gp << "     'z1_crv.dat' using 1:2 with lines lw 2.5 lc rgb '#2266CC'"<< " title 'f1(x)=10x/" << n << "',\\\n";
            gp << "     'z1_crv.dat' using 1:3 with lines lw 2.5 lc rgb '#228833'"<< " title 'f2(x)=10(x-20)/(" << n << "-20)',\\\n";
            gp << "     'z1_vx.dat' using 1:2 with points pt 7 ps 2"<< " lc rgb '#000000' title 'вершины'\n";
            gp.close();
            runGnuplot("z1.gnu");
        }
    }
}
void zad2() {//интеграл
    cout << "\nзадание 2\n";

    auto g=[](double x){ 
        return sqrt(29.0-14.0*cos(x)*cos(x)); //при n=14
    };

    double a=0,b=5;
    double gmax=sqrt(29.0);

    double exact=0;//интеграл g(x) от а до b
    int K=1000000;
    for (int i=0;i<K;i++) {
        double x=a + (b-a)*(i+0.5)/K;
        exact += g(x);
    }
    exact *= (b-a)/K; //hi
    

    cout << "точное значение (прямоугольники): " << exact << "\n";

    vector<int> Ns={100,1000,10000,100000};

    for (int N : Ns) {
        double sum1=0;
        int M=0;
        double S_rect=(b-a)*gmax;
        vector<double> xi_in,yi_in,xi_out,yi_out;
        for (int i=0;i<N;i++) {
            double x=randAB(a,b);
            double y=randAB(0,gmax);\

            if (y<=g(x)) {
                M++;
                xi_in.push_back(x); 
                yi_in.push_back(y);
            } else {
                xi_out.push_back(x); 
                yi_out.push_back(y);
            }

            sum1 += g(randAB(a,b));
        }

        double I_geom =(double)M/N*S_rect;
        double I_stat =(b-a)*sum1/N;
        cout << "N=" << N
             << "геом=" << fixed << setprecision(5) << I_geom
             << "стат=" << I_stat
             << "|ошибка геом|=" << fabs(I_geom-exact)
             << "|ошибка стат|=" << fabs(I_stat-exact) << "\n";

        if (N == 10000) {
            ofstream fin("z2_in.dat"),fout("z2_out.dat"),fcrv("z2_crv.dat");
            for (size_t k=0;k<xi_in.size();k++)  fin  << xi_in[k]  <<"\t"<< yi_in[k]  <<"\n";
            for (size_t k=0;k<xi_out.size();k++) fout << xi_out[k] <<"\t"<< yi_out[k] <<"\n";

            for (int k=0;k<=300;k++) {
                double xv=a+(b-a)*k/300.0;
                fcrv << xv <<"\t"<< g(xv) <<"\n";
            }

            ofstream gp("z2.gnu");
            gp << "set terminal pngcairo size 900,600 enhanced font 'Arial,11'\n";
            gp << "set output 'z2.png'\n";
            gp << "set title 'Задание 2: Монте-Карло,интеграл,N=10000' font 'Arial,12'\n";
            gp << "set xlabel 'x'\nset ylabel 'y'\n";
            gp << "set grid lw 1 lc rgb '#cccccc'\n";
            gp << "set key top right\n";
            gp << "plot 'z2_out.dat' using 1:2 with points pt 1 ps 0.3 lc rgb '#AAAAAA' title 'вне',\\\n";
            gp << "     'z2_in.dat'  using 1:2 with points pt 1 ps 0.3 lc rgb '#CC2200' title 'внутри',\\\n";
            gp << "     'z2_crv.dat' using 1:2 with lines lw 2 lc rgb '#2266CC' title 'sqrt(29-14cos²x)'\n";
            gp.close();
            runGnuplot("z2.gnu");
        }
    }
}

void zad3() {//вычисление пи
    cout << "\nзадание 3\n";
    double R=n;
    cout << "R=" << R <<endl;

    vector<int> Ns={100,1000,10000,100000};
    for (int N :Ns) {
        int M=0;

        vector<double> xi_in,yi_in,xi_out,yi_out;

        for (int i=0;i<N;i++) {
            double x=randAB(-R,R);
            double y=randAB(-R,R);
            if ((x)*(x) + (y)*(y) < R*R) {//точка попала в круг
                M++;
                xi_in.push_back(x); yi_in.push_back(y);
            } else {
                xi_out.push_back(x); yi_out.push_back(y);
            }
        }
        double pi_mc=4.0*M/N;//Sкруга = pi*R^2, Sкв=4R^2, Sкр/Sкв=pi/4
        double err_abs=fabs(pi_mc-M_PI);
        double err_otn=err_abs/M_PI*100.0;
        cout << "N=" << N << " M=" << M
             << "pi=" << fixed << setprecision(6) << pi_mc
             << "|err|=" << err_abs
             << "otn=" << err_otn << "%\n";

        if (N == 10000) {
            ofstream fin("z3_in.dat"),fout("z3_out.dat"),fcirc("z3_circ.dat");
            for (size_t k=0;k<xi_in.size();k++)  fin  << xi_in[k]  <<"\t"<< yi_in[k]  <<"\n";
            for (size_t k=0;k<xi_out.size();k++) fout << xi_out[k] <<"\t"<< yi_out[k] <<"\n";
            for (int k=0;k<=360;k++) {
                double phi=2.0*M_PI*k/360.0;
                fcirc << R*cos(phi) <<"\t"<< R*sin(phi) <<"\n";
            }
            ofstream gp("z3.gnu");
            gp << "set terminal pngcairo size 700,700 enhanced font 'Arial,11'\n";
            gp << "set output 'z3.png'\n";
            gp << "set title 'Задание 3: π ≈ " << fixed << setprecision(4) << 4.0*(double)M/N
               << ",N=10000' font 'Arial,12'\n";
            gp << "set xlabel 'x'\nset ylabel 'y'\n";
            gp << "set grid lw 1 lc rgb '#cccccc'\n";
            gp << "set size ratio 1\n";
            gp << "set xrange [" << -R-1 << ":" << R+1 << "]\n";
            gp << "set yrange [" << -R-1 << ":" << R+1 << "]\n";
            gp << "set key top right\n";
            gp << "plot 'z3_out.dat' using 1:2 with points pt 1 ps 0.3 lc rgb '#AAAAAA' title 'вне круга',\\\n";
            gp << "     'z3_in.dat'  using 1:2 with points pt 1 ps 0.3 lc rgb '#CC2200' title 'в круге',\\\n";
            gp << "     'z3_circ.dat' using 1:2 with lines lw 2 lc rgb '#2266CC' title 'R=" << R << "'\n";
            gp.close();
            runGnuplot("z3.gnu");
        }
    }
}

void zad4() {//площадь фигуры в полярных координатах
    cout << "\nзадание 4\n";//A*cos^2 al + B*sin^2 al = p^2
    double A=n;
    double B=n-10.0;
    cout << "A=" << A << " B=" << B << "\n";

    auto rho=[&](double phi){ 
        return sqrt(A*cos(phi)*cos(phi) + B*sin(phi)*sin(phi)); 
    };

    double rmax=sqrt(max(A,B));//граница фигуры

    double ax=-rmax,bx=rmax;
    double ay=-rmax,by=rmax;
    double S_rect=(bx-ax)*(by-ay);
    double S_exact=M_PI*sqrt(A)*sqrt(B);

    vector<int> Ns={100,1000,10000,100000};
    for (int N : Ns) {
        int M=0;
        vector<double> xi_in,yi_in,xi_out,yi_out;
        for (int i=0;i<N;i++) {
            double x=randAB(ax,bx);
            double y=randAB(ay,by);
            double ri =sqrt(x*x + y*y);
            double phi=atan2(y,x);
            double rho_phi=rho(phi);

            if (ri < rho_phi) {
                M++;
                xi_in.push_back(x); yi_in.push_back(y);
            } else {
                xi_out.push_back(x); yi_out.push_back(y);
            }

        }
        double S_mc=(double)M/N*S_rect;
        double err_abs=fabs(S_mc-S_exact);
        double err_otn=err_abs/S_exact*100.0;
        cout << "N=" << N << " M=" << M
             << "S=" << fixed << setprecision(4) << S_mc
             << "|err|=" << err_abs
             << "otn=" << err_otn << "%\n";

        if (N == 10000) {
            ofstream fin("z4_in.dat"),fout("z4_out.dat"),fcrv("z4_crv.dat");
            for (size_t k=0;k<xi_in.size();k++)fin<< xi_in[k]  <<"\t"<< yi_in[k]  <<"\n";
            for (size_t k=0;k<xi_out.size();k++) fout << xi_out[k] <<"\t"<< yi_out[k] <<"\n";

            for (int k=0;k<=500;k++) {
                double phi=2.0*M_PI*k/500.0;
                double r=rho(phi);
                fcrv << r*cos(phi) <<"\t"<< r*sin(phi) <<"\n";
            }

            ofstream gp("z4.gnu");
            gp << "set terminal pngcairo size 700,700 enhanced font 'Arial,11'\n";
            gp << "set output 'z4.png'\n";
            gp << "set title 'Задание 4: полярная фигура,N=10000' font 'Arial,12'\n";
            gp << "set xlabel 'x'\nset ylabel 'y'\n";
            gp << "set grid lw 1 lc rgb '#cccccc'\n";
            gp << "set size ratio 1\n";
            gp << "set key top right\n";
            gp << "plot 'z4_out.dat' using 1:2 with points pt 1 ps 0.3 lc rgb '#AAAAAA' title 'вне',\\\n";
            gp << "     'z4_in.dat'  using 1:2 with points pt 1 ps 0.3 lc rgb '#CC2200' title 'внутри',\\\n";
            gp << "     'z4_crv.dat' using 1:2 with lines lw 2 lc rgb '#2266CC' title 'ρ(φ)'\n";
            gp.close();
            runGnuplot("z4.gnu");
        }
    }
}

int main() {
    srand((unsigned)time(0));
    cout << fixed << setprecision(4);
    cout << "вариант n=" << n << "\n";
    zad1();
    zad2();
    zad3();
    zad4();
    return 0;
}