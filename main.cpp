#include <array>
#include <cmath>
#include <cstdint>
#include <fstream>
#include <future>
#include <iostream>
#include <thread>
#include <vector>

#include "complex.hpp"

void printUsageMessage()
{
    std::cout << "Julia set visualizer\n"
              << "Usage: ./julia <image width> <julia constant real> <julia constant imaginary> <threads>\n"
              << "Parameters: <image size>: width AND height of the output image. max: 65535\n"
              << "            <julia constant real>: real component of the julia constant\n"
              << "            <julia constant imaginary>: imaginary component of the julia constant\n"
              << "            <threads>: number of parallel computations to split the calculations into.\n"
              << "                       useful for red-ucing program runtime"
              << std::endl;
}

// Header for a TGA file. Multi-byte values are little endian.
// This must be modifed after the program starts to inject the image size.
// Note that the output image will be uncompressed.
std::array<uint8_t, 18> TGA_HEADER = {
    0x00, // ID field length
    0x00, // no color map
    0x02, // true-color image, uncompressed
    0x00, 0x00, 0x00, 0x00, 0x00, // no color map
    0x00, 0x00, // x-origin
    0x00, 0x00, // y-origin
    0x00, 0x00, // image width (uninitialized)
    0x00, 0x00, // image height (uninitialized)
    0x18, // pixel depth (RGB, 24bpp)
    0x00  // image descriptor
};

// Footer for a TGA file.
const std::array<uint8_t, 26> TGA_FOOTER = {
    0x00, 0x00, 0x00, 0x00, // extentsion area offset, not present
    0x00, 0x00, 0x00, 0x00, // developer directory offset, not present
    'T', 'R', 'U', 'E', 'V', 'I', 'S', 'I', 'O', 'N', '-', 'X', 'F', 'I', 'L', 'E', // TGA version 2 specifier
    '.', '\0' // required
};

// Wrapper for pixel coordinates.
struct Pixel
{
    int x;
    int y;
};

// An RGB color value, but in reverse order (blue, green, red) due to little endian byte ordering.
using ColorBGR = std::array<uint8_t, 3>;

const double findMinEscapeRadius(const Complex& c);
void juliaFunction(Complex& z, const Complex& c);
Complex mapCoord(const int x, const int y, const int imageSize, const double escapeRadius);
const ColorBGR mapColor(const int iterations, const int maxIterations);
const Pixel mapPixel(const int index, const int imageSize);
void workerThread(
    const int lowIndex, const int highIndex, const int imageSize,
    const Complex juliaConstant, const double escapeRadius, const double escapeRadiusSquared,
    std::promise<std::vector<uint8_t>>* returnValue, const int threadId
);

// parameters: image size, julia constant, number of threads to compute with
int main(int argc, char** argv)
{
    if (argc != 5)
    {
        printUsageMessage();
        return -1;
    }

    const int imageSize = std::atoi(argv[1]);

    TGA_HEADER[12] = static_cast<uint8_t>(imageSize % 256);
    TGA_HEADER[13] = static_cast<uint8_t>(imageSize / 256);
    TGA_HEADER[14] = TGA_HEADER[12];
    TGA_HEADER[15] = TGA_HEADER[13];

    // imageData is the sequence of bytes representing the output file.
    std::vector<uint8_t> imageData;
    imageData.insert(
        imageData.end(),
        std::make_move_iterator(TGA_HEADER.begin()),
        std::make_move_iterator(TGA_HEADER.end())
    );

    const Complex juliaConstant(
        static_cast<double>(std::atof(argv[2])),
        static_cast<double>(std::atof(argv[3]))
    );
    const double escapeRadius = findMinEscapeRadius(juliaConstant);
    const double escapeRadiusSquared = escapeRadius * escapeRadius;

    const int threads = std::atoi(argv[4]);
    const int pixelsPerThread = (imageSize * imageSize) / threads;
    std::vector<std::vector<uint8_t>> imageFragments(threads);
    std::vector<std::thread> threadPool;
    std::vector<std::promise<std::vector<uint8_t>>> promisePool(threads);
    std::vector<std::future<std::vector<uint8_t>>> futurePool;

    // splits the workload of the calculations between multiple threads
    for (int i = 0; i < threads; i++)
    {
        // this is just more concise than calculating upperBound out of the loop
        int upperBound = (i + 1) * pixelsPerThread;
        if (i == threads - 1)
            upperBound = imageSize * imageSize; // helps prevent an integer rounding error when threads != 2^n

        futurePool.push_back(promisePool[i].get_future());
        threadPool.push_back(std::thread(
            workerThread,
            i * pixelsPerThread, upperBound, imageSize,
            juliaConstant, escapeRadius, escapeRadiusSquared,
            &promisePool[i], i
        ));
    }

    // waits for the threads to return and aggregates the calculated color values into the final image
    // threads are joined in the order that they were created, meaning the image data is also in order
    for (int i = 0; i < threads; i++)
    {
        threadPool[i].join();
        auto returnValue = futurePool[i].get();
        imageData.insert(
            imageData.end(),
            std::make_move_iterator(returnValue.begin()),
            std::make_move_iterator(returnValue.end())
        );
    }

    imageData.insert(
        imageData.end(),
        std::make_move_iterator(TGA_FOOTER.begin()),
        std::make_move_iterator(TGA_FOOTER.end())
    );

    std::ofstream file("out.tga");
    for (uint8_t ch : imageData)
    {
        file << ch;
    }
    file.close();

    return 0;
}

// This function finds the minimum escape radius for a complex constant c. Starts at 5.
// Mathematically, this function approximates the minimum value of R such that R^2 - R >= |c|
const double findMinEscapeRadius(const Complex& c)
{
    const double magnitude = c.magnitude();
    const double step = 0.01;
    double R = 5.0;

    while ((R * R) - R >= magnitude)
        R -= step;

    return R;
}

// NOTE: The reference to z is modifed. (Modifying the reference to z appears to improve performance as opposed to returning a new value.)
// This function performs a single iteration of the recursive function f(z) = z^2 + c,
// where z and c are complex numbers, and c is a constant.
// Mathematically speaking, this function is used to define a Julia set.
// If the function converges for a value of z, z is in the Julia set. Otherwise, it is not in the Julia set.
inline void juliaFunction(Complex& z, const Complex& c)
{
    z *= z;
    z += c;
}

// This function maps a pixel coordinate to a coordinate on the complex plane.
// Pixel coordinates start from the bottom left of the image. x increases to the right, and y increases upward.
// The complex plane appears how one would expect: the center of the image represents 0 + 0i. The axes are as one would expect.
// The mapping is: [(0, 0), (imageSize, imageSize)] -> [(-escapeRadius, -escapeRadius), (escapeRadius, escapeRadius)]
Complex mapCoord(const int x, const int y, const int imageSize, const double escapeRadius)
{
    Complex result(0.0, 0.0);
    int half = imageSize / 2;

    double ratio = static_cast<double>(x) / half;
    if (ratio < 1.0)
        result.real = -1.0 * escapeRadius * (1.0 - ratio);
    else
        result.real = escapeRadius * (ratio - 1.0);

    ratio = static_cast<double>(y) / half;
    if (ratio < 1.0)
        result.imaginary = -1.0 * escapeRadius * (1.0 - ratio);
    else
        result.imaginary = escapeRadius * (ratio - 1.0);

    return result;
}

// This function maps a number of iterations to a color using the cubehelix coloring scheme.
// The ratio of iterations to the maximum iterations is used as input to the coloring function.
// For more information, see "A colour scheme for the display of astronomical intensity images" by D. A. Green (2011).
const ColorBGR mapColor(const int iterations, const int maxIterations)
{
    ColorBGR result = {0x00, 0x00, 0x00};

    // function parameters
    const double lambda = std::pow(static_cast<double>(iterations) / maxIterations, 0.4);
    const double startColor = 3.0;
    const double rotations = 0.5;
    const double hue = 1.0;

    // intermediary values
    const double a = hue * lambda * (1 - lambda) / 2;
    const double phi = 2.0 * M_PI * ((startColor / 3.0) + (rotations * lambda));

    // calculation of color values
    const double red = lambda + a * ((-0.14861 * std::cos(phi)) + (1.78277 * std::sin(phi)));
    const double green = lambda + a * ((-0.29227 * std::cos(phi)) + (-0.90649 * std::sin(phi)));
    const double blue = lambda + a * (1.97294 * std::cos(phi));

    // convert [0, 1] -> [0, 255]
    result[0] = static_cast<uint8_t>(blue * 255);
    result[1] = static_cast<uint8_t>(green * 255);
    result[2] = static_cast<uint8_t>(red * 255);

    return result;
}

// This function maps an index into an array of size imageSize^2 to a pixel coordinate.
// This is done since iterating over x and y values in nested loops when the threads parameter
// is not a power of 2 becomes too complex.
// The mapping is [0, (imageSize^2) - 1] -> [(0, 0), (imageSize - 1, imageSize - 1)]
const Pixel mapPixel(const int index, const int imageSize)
{
    Pixel pixel{index % imageSize, index / imageSize};
    return pixel;
}

// This function is the worker thread that calculations are divided into.
// The calculations are as follows:
//     Map an index from [lowIndex, highIndex] to a pixel coordinate
//     Map the pixel coordinate to its corresponding coordinate on the complex plane
//     Perform f(z) = z^2 + c until it converges or diverges
//     Color the corresponding pixel according to how much it diverges
void workerThread(
    const int lowIndex, const int highIndex, const int imageSize,
    const Complex juliaConstant, const double escapeRadius, const double escapeRadiusSquared,
    std::promise<std::vector<uint8_t>>* returnValue, const int threadId
) {
    std::vector<uint8_t> imageFragment;

    for (int i = lowIndex; i < highIndex; i++)
    {
        auto pixel = mapPixel(i, imageSize);
        auto coord = mapCoord(pixel.x, pixel.y, imageSize, escapeRadius);
        int iteration = 0;
        int maxIteration = 1000; 

        // if the point converges, color a black pixel.
        // else, return a color value proportional to the number of iterations it took to diverge
        while (iteration < maxIteration && coord.magnitudeSquared() <= escapeRadiusSquared)
        {
            juliaFunction(coord, juliaConstant);
            iteration++;
        }

        if (iteration == maxIteration)
        {
            for (int i = 0; i < 3; i++)
                imageFragment.push_back(0x00);
        }
        else
        {
            auto color = mapColor(iteration, maxIteration);
            for (size_t i = 0; i < color.size(); i++)
                imageFragment.push_back(color[i]);
        }
    }

    returnValue->set_value(imageFragment);
}