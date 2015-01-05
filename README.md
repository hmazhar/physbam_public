Physbam Public
==============

A modified version of the public physbam repo to support using [Cmake](http://www.cmake.org/) and newer GCC compilers
This is using the public version of [physbam](http://physbam.stanford.edu/)

Currently tested on the following:

```
Linux:
GCC 4.4.6
GCC 4.6.3
GCC 4.7.2
GCC 4.8.1  -- (Does not compile several examples in project, Linking issue)

OSX:
GCC 4.7.3
GCC 4.9.0
Apple LLVM version 6.0 (clang-600.0.56) (based on LLVM 3.5svn)

Extra packages installed using homebrew:
brew install ffmpeg openexr fftw zlib libpng libjpeg

```


Currently does not compile on:

```

Linux:
OSX:
Windows:

VS 2013
MinGW (GCC 3.4.5)
```
Build testing using GCC 4.6.3 is availible [here](https://drone.io/github.com/hmazhar/physbam_public)

Physbam is Copyright 1999-2010:

Andrew Selle, Andy Lutimirski, Avi Robinson-Mosher, Bridget Vuong, Christopher Allocco, Craig Schroeder, Don Hatch, Douglas Enright, Duc Nguyen, Eftychios Sifakis, Eilene Hao, Elliot English, Eran Guendelman, Fen Zhao, Frank Losasso, Frederic Gibou, Geoffrey Irving, Huamin Wang, Igor Neverov, Jared Go, Jeffrey Hellrung, Jeong-Mo Hong, Jerry Talton, Jiayi Chong, Jonathan Su, Jon Gretarsson, Joseph Teran, Joyce Pan, Justin Solomon, Kevin Der, Mark A. Wicks, Michael Lentine, Michael Turitzin, Mike Rodgers, Neil Molino, Nick Rasmussen, Nipun Kwatra, Paul, James White, Rachel Weinstein, Ranjitha Kumar, Robert Bridson, Robert Travis, Ron Fedkiw, Ryan Kautzman, Sergey Koltakov, Sergey Levine, Silvia Salinas-Blemker, Tamar Shinar, Unnur Gretarsdottir, Wen Zheng, Zhaosheng Bao. All rights reserved. 



