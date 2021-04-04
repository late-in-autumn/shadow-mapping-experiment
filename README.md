# CSCI 580 HW 5

## Author Information
- Sijie Bu (3102360188)
- sbu@usc.edu

## Project Information
- Built and Tested Visual Studio 2019.
  - The compiler option for debug information had to be changed to `/Zi` as the default option originated from migrating the VC6 project (`/ZI`) was incompatible with current compilers.
- Methods throw a `std::bad_alloc` exception if there is insufficient memory for some of the intermediate variables.
  - Though in practice this is highly unlikely and this this behavior is there to be a safety net.
- Includes a more generic version of the provided `ARRAY()` helper function, `MATRIX_OFFSET()`, and it is also inlined.
- Includes a template helper function for converting degrees to radiants, `TO_RADIANT()`.
  - To help with implementing this helper, the `PI` macro definition was moved from `rend.cpp` to `rend.h`, after consulting with TAs over the Slack channel.
- Contains template methods for floating point comparison.
  - These are re-implementations of floating point comparison algorithms that were originally described in the book The Art of Computer Programming: Seminumerical algorithms by Donald Knuth. See comments in the assignment code for the exact algorithm implementation.
- Contains a template function for clamping values because the `std::clamp()` function is not available with the Microsoft toolchain.
- Contains private subroutines for shading calculation, sorting vertices, calculating edge and plane equasions, computing cross products, and matrix multiplication.
- Contains additional constants for identity matrix, the eye vector, full-intensity color coefficients (for modified Gouruand shading), and the number of quadrants in a square (4, used by the texture functions).

## Notes about the procedural texture function
- The procedural texture function divides the square UV space into four quadrants:
  - 0 <= u < 0.5 && 0 <= v < 0.5
  - 0.5 <= u <= 1 && 0 <= v < 0.5
  - 0 <= u < 0.5 && 0.5 <= v <= 1
  - 0.5 <= u <= 1 && 0.5 <= v <= 1
- Similar to the image-based texture function, the procedural texture function will initialize itself when first invoked:
  - During the process, the function will use the `rand()` function to pick four random colors, one for each quadrant.
  - And subsequent invocations will use the four selected colors until the program exits or being forced to re-initialize by setting the `reset` flag.
- However, during testing, I observed that the random colors that were picked by the function seem to lack some variety.
  - I believe that this is because that the `rand()` function is not a cryptographic-grade PRNG, and therefore, generates pseudorandom numbers with relatively high level of predictability.