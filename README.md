<!-- README.md is generated from README.Rmd. Please edit that file -->

# rgl.cry

<!-- badges: start -->
<!-- badges: end -->

This package provides tools for selected area electron diffraction
(SAED) pattern and crystal structure visualization by using `rgl` and
`cry` packages. In particular `cry_demo()` and `dp_demo()` read a file
in CIF formats (Crystallographic Information Framework) and display SAED
pattern and crystal structure. The function `dp_demo()` also performs
simple simulation of powder X-ray diffraction (PXRD) pattern and it can
be saved the results to a file in the working directory.

The package has been tested on several platforms, including Linux on
Crostini with a Core™ m3-8100Y Chromebook, I found that even on this
low-powered platform, the performance was acceptable.

## Installation

``` r
install.packages("rgl.cry")
```

## Example

A CIF file is read. A reciprocal lattice map with a cell widget drawn.
In this example, a file is not specified. Therefore, a system default is
used.

``` r
dp_demo()
```

A CIF file is read. A crystal structure with a axis widget drawn.

``` r
cry_demo()
```

### Utility functions

The crystal and diffraction pattern are aligned and displayed.

``` r
align("a")
align("c*")
align("30 30") # x, y (deg)
```

Select one or more atoms or reciprocal lattice points in the window. The
labels and Miller indices of the selected atoms or lattice points will
be displayed.

``` r
> dp_demo()
[1] 1
> select()
To select points, use dragging the left mouse button.
To finish, press ESC.
.
1 0 -2
1 0 0
1 0 2

> cry_demo()
[1] 2
> select()
To select points, use dragging the left mouse button.
To finish, press ESC.
...
O3 
O3 
```

### Extras

`dp_demo()` can perform PXRD pattern simulation. The result is putted as
a file on the current directory by specifying options like this:

``` r
## Output the simulation results of the PXDR pattern as a file.
dp_demo(xrd = TRUE)
```

The file looks like this:

``` zsh
% sort -n +7 rgl.cry.dp.demo.2024-02-26_000000.dat
     h  k  l        d       absF         lp twotheta
137  0  0  0      Inf 308.117522        Inf  0.00000
136 -1  0  0 5.583222  17.389409 101.955123 15.87322
138  1  0  0 5.583222  17.389409 101.955123 15.87322
128  0 -2  0 4.453238  40.628676  63.820476 19.93782
146  0  2  0 4.453238  40.628676  63.820476 19.93782
129  1 -2  0 4.142661  17.444203  54.849532 21.44963
145 -1  2  0 4.142661  17.444203  54.849532 21.44963
```

## References

<div id="refs" class="references csl-bib-body hanging-indent">

<div id="ref-sc_table_kek" class="csl-entry">

Hanashima, T. “Contents of Atomic Scattering Factors’ Table in S. Sasaki
(1987), KEK Report 87-3.” Constructed in web page in 2001.
<https://www2.kek.jp/imss/pf/tools/sasaki/sinram/sinram.html>.

</div>

<div id="ref-mol_color" class="csl-entry">

Helmenstine, Todd. 2019. “Molecule Atom Colors – CPK Colors.”
<https://sciencenotes.org/molecule-atom-colors-cpk-colors/>.

</div>

<div id="ref-enwiki:1179864711" class="csl-entry">

Wikipedia contributors. 2023. “Atomic Radius — Wikipedia, the Free
Encyclopedia.”
<https://en.wikipedia.org/w/index.php?title=Atomic_radius&oldid=1179864711>.

</div>

</div>

## Acknowledgment

I would like to express my gratitude to the rgl, cry packages, and R for
making this work possible. The data for the scattering factors is based
on the data from the KEK Report (Hanashima). For the coloring of the
atoms, I used (Helmenstine 2019), and for the atomic radii, I used
(Wikipedia contributors 2023).
