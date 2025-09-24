#include <iostream>
#include <cmath>
#include<vector>
#include<algorithm>
#include<iterator>

#include "sample_src/Shekel/ShekelProblem.hpp"

using namespace std;

class MethodPiyavsky {
    double L;
    double epsilon; 
    double a, b;
    vector<double> x_values;
    vector<double> z_values;
    TShekelProblem *problem;
public:
    /*double objective_function(double x)
    {
        //double res = x*x;
        //double res = 2 / (x * x) - 1 / x;
        //double res = x * x - cos(18 * x);
    }
    double derivative(double x)
    {
        //return 2 * x;
        //return 2 * x + sin(18 * x) * 18;
    }*/

    double objective_function(double x) 
    {
        return problem->ComputeFunction({ x });
    }

    double derivative(double x) 
    {
        return problem->ComputeFunctionDerivatives({ x })[0];
    }

    double lipschitz_constant(double a, double b)
    {
        //double val1 = abs(derivative(a));
        //double val2 = abs(derivative(b));
        //return max(val1, val2);
        return problem->GetLipschitzConstant();
    }

    MethodPiyavsky(double epsilon_, double a_=0, double b_=0) : epsilon(epsilon_), a(a_), b(b_)
    {
        problem = new TShekelProblem();
        if (a == 0 && b == 0) 
        {
            vector<double> lb;
            vector<double> ub;
            problem->GetBounds(lb, ub);
            a = lb[0];
            b = ub[0];
        }
        L = lipschitz_constant(a, b);
        x_values.push_back(a);
        x_values.push_back(b);

        z_values.push_back(objective_function(a));
        z_values.push_back(objective_function(b));
    }

    void algorithm()
    {
        while (true)
        {
            //Для (xi-1, xi) вычиcл R(i)
            vector<double> R_values;
            for (size_t i = 1; i < x_values.size(); ++i) 
            {
                double delta = x_values[i] - x_values[i - 1];
                double R = 0.5 * L * delta - 0.5 * (z_values[i] + z_values[i - 1]);
                R_values.push_back(R);
            }     

            double max_val = R_values[0];
            size_t t = 0;
            for (size_t i = 0;i < R_values.size(); i++) 
            {
                if (R_values[i] > max_val)
                {
                    max_val = R_values[i];
                    t = i;
                }
            }

            double x_k1 = (x_values[t] + x_values[t + 1]) / 2.0 - (z_values[t + 1] - z_values[t]) / (2.0 * L); 
            //double x_k1 = (x_values[t] + x_values[t - 1]) / 2.0 - (z_values[t] - z_values[t-1]) / (2.0 * L);
            double z_k1 = objective_function(x_k1);    // Зн-е функции в новой точке

            // Вставляем новую точку и значение функции в вектора
            x_values.insert(x_values.begin() + t + 1, x_k1);
            sort(x_values.begin(), x_values.end());

            z_values.insert(z_values.begin() + t + 1, z_k1);

            //условие остановки
            if (x_values[t + 1] - x_values[t] <= epsilon) 
            {
                size_t min_index = 0;
                double min_value = z_values[0];
                for (size_t i = 1; i < z_values.size(); i++) 
                {
                    if (z_values[i] < min_value) {
                        min_value = z_values[i];
                        min_index = i;
                    }
                }
                double min_x = x_values[min_index];
                cout << "Точка минимума x = " << min_x << endl;
                cout <<"Кол-во точек:" << x_values.size() << endl;
                double min_z = z_values[min_index];
                cout << "Значение функции в точке: " << min_z << endl;
                break;
            }
        }
    }
};



int main() 
{
    setlocale(LC_ALL, "Russian");
    TShekelProblem problem;
    double epsilon = 0.01;
    double a, b;
    cout << "Введите левую границу интервала (a): ";
    cin >> a;
    cout << "Введите правую границу интервала (b): ";
    cin >> b;
    MethodPiyavsky method(epsilon, a, b);
    method.algorithm();

    return 0;

}




/*int main()
{
    setlocale(LC_ALL, "Russian");
    double epsilon = 0.001;
    double a, b;

    //cout << "Введите точность: ";
    //cin >> epsilon;
    cout << "Введите левую границу интервала (a): ";
    cin >> a;
    cout << "Введите правую границу интервала (b): ";
    cin >> b;
    MethodPiyavsky method(epsilon, a, b);
    method.algorithm();

    return 0;
}*/
/*int main()
{
    setlocale(LC_ALL, "Russian");
    TShekelProblem problem;
    vector<double> point = { 2.5 };
    //double value = problem.Compute(0, point);
    double value = problem.ComputeFunction(point);
    cout << "Значение функции Шекеля в точке x = " << point[0] << " : " << value << "\n";

    return 0;
}*/


