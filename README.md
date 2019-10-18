# On the Formfactor between Two Polygons

<p align="center"><img src="caltechFF.png" width="600" /></p>

## Abstract
Form factors are used in radiosity to describe the fraction of diffusely reflected light leaving one surface and arriving at another. They are a fundamental geometric property used for computation. Many special configurations admit closed form solutions. However, the important case of the form factor between two polygons in three space has had no known closed form solution. We give such a solution for the case of general (planar, convex or concave, possibly containing holes) polygons.

## Authors
* Peter Schröder
* Pat Hanrahan

## Changes
* CMake project & build support
* MSVC compatibility fixes
* data export for "geometry for two rectangles sharing a common edge" figure

## Figure generator
```
ffData = Import["caltechFF.txt", "Table"];

ListLinePlot[Transpose[ffData], PlotRange -> {{0, 528}, {0, 1.025}},
 ImageSize -> {840 - 24, 512 - 24}, AxesLabel -> {Style["\!\(\*
StyleBox[\"\[Theta]\", \"InlineCode\"]\)", Black,
    FontFamily -> "Cambria", FontSize -> 24],
   Style["F12", Black, FontFamily -> "Cambria", FontSize -> 24]},
 PlotStyle -> Directive[Black, Thick],
 AxesStyle -> Directive[Black, Thick, 16],
 Ticks -> {Table[{i, 180 i/512 \[Degree]}, {i, 512/18, 512, 512/18}],
   Table[N[i], {i, 1/5, 1, 1/5}]}, ImageMargins -> 16]
```

## References
 * [Caltech Multi-Res Modeling Group - Publications](http://www.multires.caltech.edu/pubs/pubs.htm)
