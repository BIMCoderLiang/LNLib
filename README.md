## Introduction
**LNLib is a C++ NURBS Algorithms Library.** <br/>These algorithms are primary referenced from [The NURBS Book 2nd Edition](https://link.springer.com/book/10.1007/978-3-642-97385-7). <br/>The APIs are re-designed to make it more friendly to users.

<img src="assets/LNLib.png" width=400 height=200>

## Run LNLib
Please run build.bat first to construct C++ solution by CMake.

## Features
Basic Elements:
- UV
- XYZ
- XYZW
- Matrix4d
- LNObject

Algorithms in ***The Nurbs Book***:
|Chapter|Content|
|--|--|
|***Chapter 1***  | Basis Function Computation |
|***Chapter 1 to 4***  | Bezier/B-Spline/NURBS Curve and Surface |
|***Chapter 5***  | Curve and Surface Decomposition</br>Knot Insertion/Refinement/Removal</br>Degree Elevation and Reduction |
|***Chapter 6***  | Curve/Surface Point Inversion</br>Surface Tangent Vector Inversion</br>Curve/Surface Reparameterization</br>Curve Transform and Reverse</br> Surface Swap and Reverse|
|***Chapter 7***  | Create Arc/Conic Curve |
|***Chapter 8***  | Create Bilinear/Cylindrical/Ruled/Revolved/CornerFillet Surface |
|***Chapter 9***  | Global/Local Curve/Surface Interpolation and Approximation |
|***Chapter 10***  | Create Swung/Loft/Sweep/Gordon/Coons Surface |
|***Chapter 11***  | Curve Modification in Control Point Locations or Weight Values |
|***Chapter 12***  | Curve Clamp/UnClamp/IsClamp </br> KnotVector IsUniform </br> Curve IsClosed/IsPeriodic|

Additional Algorithms:
|Description|Content|
|--|--|  
|***Basic Properties***  | Curve/Surface Curvature and Normal</br>Curve Split/Segment/Merge/Offset</br>Curve IsLinear/IsArc</br>Curve Approximate Length</br>Surface Approximate Area |
|***Curve Creation***  | Create Line/Cubic Hermite |
|***Tessellation***  | Curve Tessellation </br> Surface Triangulation|

## Visualization
[LNLibViewer](https://github.com/BIMCoderLiang/LNLibViewer) based on [VTK](https://vtk.org/)

<img src="assets/curve.png" width=400 height=400><img src="assets/surface.png" width=400 height=400>

## NURBS Fitting by Neural Network
[ND-LNLib](https://github.com/BIMCoderLiang/NURBS-Diff-with-LNLib) based on [LibTorch](https://pytorch.org/cppdocs/installing.html) (PyTorch C++ version)

<img src="assets/aicurve.jpg" width=800 height=400><img src="assets/aisurface.png" width=800 height=400>

## Understand LNLib by LLM
Welcome to use https://deepwiki.com/BIMCoderLiang/LNLib powered by Devin.

## Contributing
Welcome join this project including discussions in **Issues** and make **Pull requests**.

<a href="https://contributors-img.web.app/image?repo=BIMCoderLiang/LNLib">
  <img src="https://contributors-img.web.app/image?repo=BIMCoderLiang/LNLib" />
</a>

## Author
LNLib is created by Yuqing Liang (BIMCoder Liang).

- bim.frankliang@foxmail.com
- 微信公众号：**BIMCoder**

## License
The source code is published under [LGPL 2.1](https://www.gnu.org/licenses/), the license is available [here](LICENSE).

## Primary Reference
[The NURBS Book 2nd Edition](https://link.springer.com/book/10.1007/978-3-642-97385-7) by **Les Piegl & Wayne Tiller**
