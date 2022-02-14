<p align="center">
  <img src="/img/sigma-l.PNG">
</p>

# schwarzschild-double
Partial absorption cross section of Schwarzschild BH calculation using [GSL](https://www.gnu.org/software/gsl/) and [glib](https://docs.gtk.org/glib/).

# GUI
The GUI was made with [Tcl/Tk](https://www.tcl.tk/) and [gnuplot](http://www.gnuplot.info/)

![gui](/img/gui-1.PNG)

![gui](/img/gui-2.PNG)

![tortoise coordinates](/img/tortoise.PNG)

![wave and potential](/img/wave-potential.PNG)

![transmission and reflection coefficients](/img/coefficients.PNG)

![partial absortion cross section](/img/sigma-l.PNG)

To use the GUI you need schwarzschild.exe and gnuplot.exe in the same directory.

# Compilation
If you want to compile the code you'll need the [GSL](https://www.gnu.org/software/gsl/) and [glib](https://docs.gtk.org/glib/) libraries. I've done all the development in Windows so I've used the [MSYS2](https://www.msys2.org/) system. You can use [Chocolatey](https://chocolatey.org/) or install MSYS2 directly.

To install the packages I've used [pacman](https://archlinux.org/pacman/pacman.8.html). To compile you just need a c99 compiler. I've used [CodeBlocks](https://www.codeblocks.org/) as IDE.

# Documentation
For  more information you can consult the pdf's in the *doc* folder (English version is just an automatic translation)
