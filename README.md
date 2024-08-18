# Triangle Boolean Subtraction Algorithm

This repository contains the implementation of a Triangle Boolean Subtraction algorithm developed to create 
holes in navigation mesh polygons. The algorithm is designed to handle cases where single precision floating
point is required for performance reasons, making it suitable for use in game development.

## Overview
The need for this algorithm arose from the necessity to subtract one triangle from another to create obstacles 
in navigation meshes. This process involves complex geometric calculations, and after extensive research and 
testing, this solution was developed.

## Features

- **Triangle Subtraction**: The algorithm performs boolean subtraction on triangles, simplifying the problem to 20 possible cases.
- **Optimized for Performance**: The implementation uses single precision floating point for faster performance, crucial in real-time applications like games.
- **Reliable and Tested**: The algorithm has been rigorously tested with billions of cases, ensuring stability and accuracy.

## Usage
The algorithm is implemented in C/C++, and can be integrated into projects requiring navigation mesh manipulation or other geometric processing tasks.

## Installation
The best way to use it is to integrate it within your application.

## Implementation Details
For a detailed explanation of the algorithm, including step-by-step implementation, please refer to 
[this blog post](___).

## Short Overview
The algorithm starts by checking for overlap between the triangles using SAT, reorients them if necessary, and 
proceeds with collinearity checks. Once the overlap is confirmed, an overlap polygon is created, cleaned up, and 
inserted into the subject triangle. Additional checks ensure that the points are correctly processed and handled 
to avoid stacking or precision errors.

The algorithm constructs tables to keep track of vertex information, including whether a vertex is outside, if it's a 
crossing point, and its relationship to other vertices. These tables are then used to generate the final resulting polygons.

## License
This project is licensed under the MIT License - see the [LICENSE](LICENSE) file for details.
