// Simple overloads for complex number operations

#include <cmath>

class Complex
{
public:
    double real;
    double imaginary;

    Complex(double real, double imaginary)
    {
        this->real = real;
        this->imaginary = imaginary;
    }

    Complex& operator+=(const Complex& rhs)
    {
        real += rhs.real;
        imaginary += rhs.imaginary;

        return *this;
    }

    Complex& operator-=(const Complex& rhs)
    {
        real -= rhs.real;
        imaginary -= rhs.imaginary;

        return *this;
    }

    Complex& operator*=(const Complex& rhs)
    {
        double newReal = (real * rhs.real) - (imaginary * rhs.imaginary);
        double newImaginary = (real * rhs.imaginary) + (imaginary * rhs.real);

        real = newReal;
        imaginary = newImaginary;

        return *this;
    }

    Complex& operator/=(const Complex& rhs)
    {
        double divisor = (rhs.real * rhs.real) + (rhs.imaginary * rhs.imaginary);
        double newReal = ((real * rhs.real) + (imaginary * rhs.imaginary)) / divisor;
        double newImaginary = ((rhs.real * imaginary) - (real * rhs.imaginary)) / divisor;
        
        real = newReal;
        imaginary = newImaginary;

        return *this;
    }

    Complex conjugate()
    {
        return Complex(real, -imaginary);
    }

    double magnitude()
    {
        return std::sqrt(std::pow(real, 2.0) + std::pow(imaginary, 2.0));
    }
};