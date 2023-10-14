## Introduction
**LNLib is a C++ NURBS Algorithms Library.** <br/>These algorithms are primary referenced from [The NURBS Book 2nd Edition](https://link.springer.com/book/10.1007/978-3-642-97385-7). <br/>The APIs are re-designed to make it more friendly to users.

## Run LNLib
Please run build.bat first and it will construct solution by Cmake.

## Features
Basic Elements:
- UV
- XYZ
- XYZW
- Matrix4d

NURBS Algorithms:
- Basis Function Computation
- Bezier/B-Spline/NURBS Curve and Surface
- Curve and Surface Decomposition
- Knot Insertion/Refinement/Removal
- Degree Elevation and Reduction
- Curve/Surface Point Inversion
- Surface Tangent Vector Inversion
- Curve Reparameterization
- Curve Transform and Reverse
- Create Arc/Conic Curve
- Create Bilinear/Cylindrical/Ruled/Revolved Surface
- Global/Local Curve/Surface Interpolation and Approximation
- Curve Modification in Control Point Locations or Weight Values
- Clamp/UnClamp Nurbs Curve

This library is **in progress.**

## Roadmap
- Algorithms of Chapter 10 & 11 in The NURBS Book 2nd Edition
- Test cases for all algorithms

## Contributing
Welcome join this project including discussions in **Issues** and make **Pull requests**.

## Author
LNLib is created by Yuqing Liang (BIMCoder Liang).

- bim.frankliang@foxmail.com
- 微信公众号：**BIMCoder梁老师**

## License
The source code is published under [GNU General Public License v3.0](https://www.gnu.org/licenses/), the license is available [here](LICENSE).

## Primary Reference
[The NURBS Book 2nd Edition](https://link.springer.com/book/10.1007/978-3-642-97385-7) by **Les Piegl & Wayne Tiller**
