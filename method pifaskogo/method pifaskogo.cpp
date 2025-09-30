﻿#include <iostream>
#include <cmath>
#include<vector>
#include<algorithm>
#include<iterator>
#include <fstream>
#include "sample_src/Shekel/ShekelProblem.hpp"
using namespace std;

class MethodPiyavsky {
    double L;
    double epsilon;
    double a, b;
    vector<double> x_values;
    vector<double> z_values;
    TShekelProblem* problem;
    double min_x;
    double min_value;
    size_t point_count;
    double true_opt_point;
    double true_opt_value;
public:
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
        return problem->GetLipschitzConstant();
    }

    MethodPiyavsky(double epsilon_, double a_ = 0, double b_ = 0, int function_number = 0) : epsilon(epsilon_), a(a_), b(b_)
    {
        problem = new TShekelProblem(function_number);
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

        vector<double> true_point = problem->GetOptimumPoint();
        true_opt_point = true_point[0];
        true_opt_value = problem->GetOptimumValue();
    }

    ~MethodPiyavsky()
    {
        delete problem;
    }

    void getResults(double& x_min, double& f_min, size_t& count, double& true_x, double& true_f)
    {
        x_min = min_x;
        f_min = min_value;
        count = point_count;
        true_x = true_opt_point;
        true_f = true_opt_value;
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
                min_value = z_values[0];
                for (size_t i = 1; i < z_values.size(); i++)
                {
                    if (z_values[i] < min_value) {
                        min_value = z_values[i];
                        min_index = i;
                    }
                }
                min_x = x_values[min_index]; 
                point_count = x_values.size();

                /*double min_x = x_values[min_index];
                cout << "Точка минимума x = " << min_x << endl;
                cout << "Кол-во точек: " << x_values.size() << endl;
                cout << "Значение функции в точке: " << min_value << endl;

                vector<double> true_opt_point = problem->GetOptimumPoint();
                double true_opt_value = problem->GetOptimumValue();
                cout << "Истинная точка минимума: " << true_opt_point[0] << endl;
                cout << "Истинное значение функции: " << true_opt_value << endl;
                cout << "Погрешность по x: " << abs(min_x - true_opt_point[0]) << endl;
                cout << "Погрешность по f: " << abs(min_value - true_opt_value) << endl;
                */
                break;
            }
        }
    }
};

int main()
{
    setlocale(LC_ALL, "Russian");
    TShekelProblem problem;
    int function_number;
    double epsilon = 0.01;
    double a, b;
    double min_x, min_value;
    size_t point_count;
    double true_x;
    double true_f;
    //cout << "Введите номер функции Шекеля (0-999): ";
    //cin >> function_number;
    cout << "Введите левую границу интервала (a): ";
    cin >> a;
    cout << "Введите правую границу интервала (b): ";
    cin >> b;
    //MethodPiyavsky method(epsilon, a, b, function_number);
    //method.algorithm();

    ofstream file("results.csv");
    if (!file.is_open())
    {
        cout << "Ошибка с файлом";
        return 1;
    }
    file << "Номер ф-ции; Точка минимума; Знач ф-ции; Кол-во точек; Ист.знач.минимума; Ист.знач.ф-ции; Погрешность по x;Погрешность по f " << endl;
    for (function_number = 0;function_number < 999;function_number++)
    {
        MethodPiyavsky method(epsilon, a, b, function_number);
        method.algorithm();
        method.getResults(min_x, min_value, point_count, true_x, true_f);
        double error_x = abs(min_x - true_x);
        double error_f = abs(min_value - true_f);
        file << function_number << ";" << min_x << ";" << min_value << ";" << point_count << ";"<< true_x <<";"<< true_f<<";"<<error_x<<";"<<error_f<< endl;
    }

    file.close();
    cout << "Файл успешно создан и заполнен";
    return 0;
}
