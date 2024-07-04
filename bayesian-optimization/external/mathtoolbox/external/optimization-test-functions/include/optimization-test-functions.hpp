/*
 MIT License

 Copyright (c) 2019 Yuki Koyama

 Permission is hereby granted, free of charge, to any person obtaining a copy
 of this software and associated documentation files (the "Software"), to deal
 in the Software without restriction, including without limitation the rights
 to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
 copies of the Software, and to permit persons to whom the Software is
 furnished to do so, subject to the following conditions:

 The above copyright notice and this permission notice shall be included in all
 copies or substantial portions of the Software.

 THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
 AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
 OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
 SOFTWARE.
 */

#ifndef OPTIMIZATION_TEST_FUNCTIONS_HPP_
#define OPTIMIZATION_TEST_FUNCTIONS_HPP_

#include <cassert>
#include <iostream>
#if defined(OTF_WITH_EIGEN)
#include <Eigen/Core>
#endif

namespace otf
{
    enum class FunctionType
    {
        Beale,
        Rosenbrock,
        Sphere,
    };

    inline double GetValue(const double x[], const int n, const FunctionType type)
    {
        switch (type)
        {
            case FunctionType::Beale:
            {
                assert(n == 2);
                const double a = 1.500 - x[0] + x[0] * x[1];
                const double b = 2.250 - x[0] + x[0] * x[1] * x[1];
                const double c = 2.625 - x[0] + x[0] * x[1] * x[1] * x[1];
                return a * a + b * b + c * c;
            }
            case FunctionType::Rosenbrock:
            {
                assert(n >= 2);
                double value = 0.0;
                for (int i = 0; i < n - 1; ++ i)
                {
                    value += 100.0 * (x[i + 1] - x[i] * x[i]) * (x[i + 1] - x[i] * x[i]) + (1.0 - x[i]) * (1.0 - x[i]);
                }
                return value;
            }
            case FunctionType::Sphere:
            {
                double value = 0.0;
                for (int i = 0; i < n; ++ i)
                {
                    value += x[i] * x[i];
                }
                return value;
            }
        }
    }

    inline void GetGrad(const double x[],
                        const int n,
                        const FunctionType type,
                        double grad[])
    {
        switch (type)
        {
            case FunctionType::Beale:
            {
                assert(n == 2);
                const double a = 1.500 - x[0] + x[0] * x[1];
                const double b = 2.250 - x[0] + x[0] * x[1] * x[1];
                const double c = 2.625 - x[0] + x[0] * x[1] * x[1] * x[1];
                grad[0] = 2.0 * a * (- 1.0 + x[1]) + 2.0 * b * (- 1.0 + x[1] * x[1]) + 2.0 * c * (- 1.0 + x[1] * x[1] * x[1]);
                grad[1] = 2.0 * a * x[0] + 2.0 * b * 2.0 * x[0] * x[1] + 2.0 * c * 3.0 * x[0] * x[1] * x[1];
                return;
            }
            case FunctionType::Rosenbrock:
            {
                assert(n >= 2);
                grad[0] = - 400.0 * x[0] * (x[1] - x[0] * x[0]) + 2.0 * (x[0] - 1.0);
                for (int i = 1; i < n - 1; ++ i)
                {
                    grad[i] = 200.0 * (x[i] - x[i - 1] * x[i - 1]) - 400.0 * x[i] * (x[i + 1] - x[i] * x[i]) + 2.0 * (x[i] - 1.0);
                }
                grad[n - 1] = 200.0 * (x[n - 1] - x[n - 2] * x[n - 2]);
                return;
            }
            case FunctionType::Sphere:
            {
                for (int i = 0; i < n; ++ i)
                {
                    grad[i] = 2.0 * x[i];
                }
                return;
            }
        }
    }

    inline void GetSolution(const int n,
                            const FunctionType type,
                            double x[])
    {
        switch (type)
        {
            case FunctionType::Beale:
            {
                assert(n == 2);
                x[0] = 3.0;
                x[1] = 0.5;
                return;
            }
            case FunctionType::Rosenbrock:
            {
                assert(n >= 2);
                for (int i = 0; i < n; ++ i)
                {
                    x[i] = 1.0;
                }
                return;
            }
            case FunctionType::Sphere:
            {
                for (int i = 0; i < n; ++ i)
                {
                    x[i] = 0.0;
                }
                return;
            }
        }
    }

#if defined(OTF_WITH_EIGEN)
    inline double GetValue(const Eigen::VectorXd& x, const FunctionType type)
    {
        return GetValue(x.data(), x.size(), type);
    }

    inline Eigen::VectorXd GetGrad(const Eigen::VectorXd& x, const FunctionType type)
    {
        Eigen::VectorXd grad(x.size());
        GetGrad(x.data(), x.size(), type, grad.data());
        return grad;
    }

    inline Eigen::VectorXd GetSolution(const int n, const FunctionType type)
    {
        Eigen::VectorXd x(n);
        GetSolution(n, type, x.data());
        return x;
    }
#endif
}

#endif
