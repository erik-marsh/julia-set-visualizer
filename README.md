# Julia set visualizer

C++ program that outputs a colored Julia fractal. The program evaluates the fractal by iterating the recursive function f(z) = z^2 + c, where z and c are complex numbers. In particular, z is a point on the complex plane. For a more detailed and technial explanation, see [this Wikipedia article](https://en.wikipedia.org/wiki/Julia_set).

The colors represent the number of iterations the point on the complex plane takes to diverge from a certain threshold. Darker colors represent less iterations, and brighter colors represent more iterations. Black indicates a point that converges (theoretically, that feature is not yet implemented correctly). The coloring scheme is known as cubehelix coloring. You can read [a paper on cubehelix coloring](https://astron-soc.in/bulletin/11June/289392011.pdf), or visit a [related page](http://www.mrao.cam.ac.uk/~dag/CUBEHELIX/).

## Usage
`make`

`./julia <image size> <julia constant real part> <julia constant imaginary part> <number of threads>`

For example, `./julia 2048 -0.8 0.156 4` will output a file named `out.tga` that has sixe 2048x2048. This image shows the Julia fractal for f(z) = z^2 - 0.8 + 156i