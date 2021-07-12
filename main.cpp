#include <array>
#include <cmath>
#include <cstdint>
#include <cstdlib>
#include <fstream>
#include <future>
#include <iterator>
#include <iostream>
#include <string>
#include <thread>
#include <vector>

#include "complex.hpp"

const std::string USAGE_MESSAGE = "Usage: ./julia <image width> <julia constant real> <julia constant imaginary> <threads>";

// little endian
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

const std::array<uint8_t, 26> TGA_FOOTER = {
    0x00, 0x00, 0x00, 0x00, // extentsion area offset, not present
    0x00, 0x00, 0x00, 0x00, // developer directory offset, not present
    'T', 'R', 'U', 'E', 'V', 'I', 'S', 'I', 'O', 'N', '-', 'X', 'F', 'I', 'L', 'E', // TGA version 2 specifier
    '.', '\0' // required
};

// [(0, 0), (imageSize, imageSize)] -> [(-escapeRadius, -escapeRadius), (escapeRadius, escapeRadius)]
Complex mapCoord(const int x, const int y, const int imageSize, const double escapeRadius)
{
    Complex result(0.0, 0.0);
    int half = imageSize / 2;

    double ratio = static_cast<double>(x) / half;
    if (ratio < 1.0)
    {
        result.real = -1.0 * escapeRadius * (1.0 - ratio);
    }
    else
    {
        result.real = escapeRadius * (ratio - 1.0);
    }

    ratio = static_cast<double>(y) / half;
    if (ratio < 1.0)
    {
        result.imaginary = -1.0 * escapeRadius * (1.0 - ratio);
    }
    else
    {
        result.imaginary = escapeRadius * (ratio - 1.0);
    }

    return result;
}

const double findMinEscapeRadius(const Complex& c)
{
    double magnitude = c.magnitude();
    double R = 5.0;
    double step = 0.01;

    while ((R * R) - R >= magnitude)
    {
        R -= step;
    }

    return R;
}

// the reference to z is modified. this speeds up things immensely
inline void juliaFunction(Complex& z, const Complex& c)
{
    z *= z;
    z += c;
}

// low -> blue; high -> red
const std::array<uint8_t, 3> mapColor(const int iterations, const int maxIterations)
{
    std::array<uint8_t, 3> result = {0x00, 0x00, 0x00};

    // cubehelix coloring ("A colour scheme for the display of astronomical intensity images" by D. A. Green, 2011)
    const double lambda = std::pow(static_cast<double>(iterations) / maxIterations, 0.4);
    const double startColor = 3.0;
    const double rotations = 0.5;
    const double hue = 1.0;

    // intermediary values
    const double a = hue * lambda * (1 - lambda) / 2;
    const double phi = 2.0 * M_PI * ((startColor / 3.0) + (rotations * lambda));

    const double red = lambda + a * ((-0.14861 * std::cos(phi)) + (1.78277 * std::sin(phi)));
    const double green = lambda + a * ((-0.29227 * std::cos(phi)) + (-0.90649 * std::sin(phi)));
    const double blue = lambda + a * (1.97294 * std::cos(phi));

    // convert [0, 1] -> [0, 255]
    result[0] = static_cast<uint8_t>(blue * 255);
    result[1] = static_cast<uint8_t>(green * 255);
    result[2] = static_cast<uint8_t>(red * 255);

    return result;
}

// [0, (imageSize^2) - 1] -> [(0, 0), (imageSize - 1, imageSize - 1)]
const std::pair<int, int> mapPixel(const int index, const int imageSize)
{
    std::pair<int, int> pixel;
    pixel.first = index % imageSize;
    pixel.second = index / imageSize;
    return pixel;
}

void workerThread(
    const int lowIndex, const int highIndex, const int imageSize,
    const Complex juliaConstant, const double escapeRadius, const double escapeRadiusSquared,
    std::promise<std::vector<uint8_t>>* returnValue, const int threadId
) {
    std::vector<uint8_t> imageFragment;

    // std::cout << threadId << ": lowIndex: " << lowIndex << std::endl;
    // std::cout << threadId << ": highIndex: " << highIndex << std::endl;
    // std::cout << threadId << ": imageSize: " << imageSize << std::endl;
    // std::cout << threadId << ": juliaConstant.real: " << juliaConstant.real << std::endl;
    // std::cout << threadId << ": juliaConstant.imaginary: " << juliaConstant.imaginary << std::endl;
    // std::cout << threadId << ": escapeRadius: " << escapeRadius << std::endl;
    // std::cout << threadId << ": escapeRadiusSquared: " << escapeRadiusSquared << std::endl;

    for (int i = lowIndex; i < highIndex; i++)
    {
        auto pixel = mapPixel(i, imageSize);
        auto coord = mapCoord(pixel.first, pixel.second, imageSize, escapeRadius);
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

    //std::cout << "imageFragment.size(): " << imageFragment.size() << std::endl;
    returnValue->set_value(imageFragment);
}

// parameters: image size, julia constant, number of threads to compute with
int main(int argc, char** argv)
{
    if (argc != 5)
    {
        std::cout << USAGE_MESSAGE << std::endl;
        return -1;
    }

    const int imageSize = std::atoi(argv[1]);

    TGA_HEADER[12] = static_cast<uint8_t>(imageSize % 256);
    TGA_HEADER[13] = static_cast<uint8_t>(imageSize / 256);
    TGA_HEADER[14] = TGA_HEADER[12];
    TGA_HEADER[15] = TGA_HEADER[13];
    std::vector<uint8_t> imageData;
    imageData.insert(
        imageData.end(),
        std::make_move_iterator(TGA_HEADER.begin()),
        std::make_move_iterator(TGA_HEADER.end())
    );

    Complex juliaConstant(
        static_cast<double>(std::atof(argv[2])),
        static_cast<double>(std::atof(argv[3]))
    );
    const double escapeRadius = findMinEscapeRadius(juliaConstant);
    const double escapeRadiusSquared = escapeRadius * escapeRadius;

    const int threads = std::atoi(argv[4]);
    std::vector<std::vector<uint8_t>> imageFragments(threads);
    std::vector<std::thread> threadPool;
    std::vector<std::promise<std::vector<uint8_t>>> promisePool(threads);
    std::vector<std::future<std::vector<uint8_t>>> futurePool;
    const int pixelsPerThread = (imageSize * imageSize) / threads;

    for (int i = 0; i < threads - 1; i++)
    {
        futurePool.push_back(promisePool[i].get_future());
        threadPool.push_back(std::thread(
            workerThread,
            i * pixelsPerThread,
            ((i + 1) * pixelsPerThread),
            imageSize,
            juliaConstant,
            escapeRadius,
            escapeRadiusSquared,
            &promisePool[i],
            i
        ));
    }

    futurePool.push_back(promisePool[promisePool.size() - 1].get_future());
    threadPool.push_back(std::thread(
        workerThread,
        (threads - 1) * pixelsPerThread,
        imageSize * imageSize,
        imageSize,
        juliaConstant,
        escapeRadius,
        escapeRadiusSquared,
        &promisePool[promisePool.size() - 1],
        promisePool.size() - 1
    ));

    // for (int y = 0; y < imageSize; y++)
    // {
    //     for (int x = 0; x < imageSize; x++)
    //     {
    //         auto coord = mapCoord(x, y, imageSize, escapeRadius);
    //         int iteration = 0;
    //         int maxIteration = 1000; 

    //         // if the point converges, color a black pixel.
    //         // else, return a color value proportional to the number of iterations it took to diverge
    //         while (iteration < maxIteration && coord.magnitudeSquared() <= escapeRadiusSquared)
    //         {
    //             juliaFunction(coord, juliaConstant);
    //             iteration++;
    //         }

    //         if (iteration == maxIteration)
    //         {
    //             for (int i = 0; i < 3; i++)
    //                 imageData.push_back(0x00);
    //         }
    //         else
    //         {
    //             auto color = mapColor(iteration, maxIteration);
    //             for (size_t i = 0; i < color.size(); i++)
    //                 imageData.push_back(color[i]);
    //         }
    //     }
    // }

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