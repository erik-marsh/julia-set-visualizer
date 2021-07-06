#include <array>
#include <cmath>
#include <cstdint>
#include <cstdlib>
#include <fstream>
#include <iterator>
#include <iostream>
#include <vector>

#include "complex.hpp"

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
Complex mapCoord(int x, int y, int imageSize, double escapeRadius)
{
    Complex result(0.0, 0.0);
    int half = imageSize / 2;

    double ratio = (double) x / half;
    if (ratio < 1.0)
    {
        result.real = -1.0 * escapeRadius * (1.0 - ratio);
    }
    else
    {
        result.real = escapeRadius * (ratio - 1.0);
    }

    ratio = (double) y / half;
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

double findMinEscapeRadius(Complex& c)
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
void juliaFunction(Complex& z, Complex& c)
{
    z *= z;
    z += c;
}

std::array<uint8_t, 3> mapColor(int iterations)
{
    std::array<uint8_t, 3> result = {0x00, 0x00, 0x00};
    result[0] = static_cast<uint8_t>(iterations);
    result[1] = static_cast<uint8_t>(iterations >> 8);
    result[2] = static_cast<uint8_t>(iterations >> 16);
    return result;
}

// parameters: image size, julia set function
int main(int argc, char** argv)
{
    if (argc != 2)
    {
        return -1;
    }

    int imageSize = std::atoi(argv[1]);

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

    Complex juliaConstant(-0.8, 0.156);
    double escapeRadius = findMinEscapeRadius(juliaConstant);
    double escapeRadiusSquared = escapeRadius * escapeRadius;

    for (int y = 0; y < imageSize; y++)
    {
        for (int x = 0; x < imageSize; x++)
        {
            auto coord = mapCoord(x, y, imageSize, escapeRadius);
            int iteration = 0;
            int maxIteration = 1000; 

            // if the point converges, color a black pixel.
            // else, return a color value proportional to the number of iterations it took to diverge
            while (iteration <= maxIteration && coord.magnitudeSquared() <= escapeRadiusSquared)
            {
                juliaFunction(coord, juliaConstant);
                iteration++;
            }

            if (iteration == maxIteration)
            {
                for (int i = 0; i < 3; i++)
                    imageData.push_back(0x80);
            }
            else
            {
                auto color = mapColor(iteration);
                for (size_t i = 0; i < color.size(); i++)
                    imageData.push_back(color[i]);
            }
        }
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