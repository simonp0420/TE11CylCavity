# TE11CylCavity - Find ϵᵣ and Q for disk-shaped dielectric sample inserted in a TE₁₁ₚ cavity

## How to install Julia for using this package
[Julia](julialang.org) is a modern, powerful computer language for technical computing.  It has all the nice features of 
Matlab's linear algebra syntax, along with modern programming features found in languages like Python.  Best of all, 
correctly written Julia is as fast as compiled Fortran, C, or C++.

To install Julia, follow the instructions [here](https://julialang.org/install/).  On Windows, this is as simple as installing
Julia from the Windows Store (it will actually install both Juliaup and Julia).

## Running Julia

The best experience of using Julia comes from use of VS Code with the 
[Julia Extension](https://code.visualstudio.com/docs/languages/julia) as an IDE (integrated development environment).
However, if you don't want to spend the time learning this IDE, you can start Julia from an open command-line window, preferably
the Windows Terminal app.  Note that unlike Matlab, Julia requires an explicit command-line argument to enable multithreading.
So you should start your Julia session like this:

```
PS C:\Users\peter\Documents\julia\packages\TE11CylCavity> julia -t auto
               _
   _       _ _(_)_     |  Documentation: https://docs.julialang.org
  (_)     | (_) (_)    |
   _ _   _| |_  __ _   |  Type "?" for help, "]?" for Pkg help.
  | | | | | | |/ _` |  |
  | | |_| | | | (_| |  |  Version 1.11.5 (2025-04-14)
 _/ |\__'_|_|_|\__'_|  |  Official https://julialang.org/ release
|__/                   |

julia>
```

The `-t auto` argument tells Julia that you want it to determine the number of threads to use automatically, based on your CPU 
architecture.

## Initial configuration of Julia 

This only needs to be done once (unless changes are made to the `TE11CylCavity` package, more on that later).
At the Julia prompt, enter the lines shown below:

```
using Pkg: Pkg
Pkg.activate("Cavity", shared=true)
Pkg.add(url="https://github.com/simonp0420/TE11CylCavity")
Pkg.add("MKL") # This line only needed if you are running on an Intel CPU
Pkg.precompile()
```

There will be quite a bit of output after the third line.  You can ignore
any warnings regarding packages that could be updated.  From now on in any Julia script where you want to 
use the functions in `TE11CylCavity`,  the first lines of the script should be

```Julia
using Pkg: Pkg
Pkg.activate("Cavity", shared=true)
using MKL # Only if running on an Intel CPU
using TE11CylCavity
```

**Note that the first two lines above are not currently present in the scripts located in the `scripts` directory of the 
repository.  You will have to edit the scripts for your use and add them yourself!**
Once the scripts are properly edited, you can run them via a command at the Julia prompt like this:

```Julia
include("test_findcase6fresQ.jl")
```

When you are ready to exit your Julia session, you can either type Ctrl-D or `exit()` at the Julia prompt, to take you
back to the OS prompt.

### More details on updating `TE11CylCavity`

Recall that when we added `MKL` to the shared environment named "Cavity", we did not have to specify the full Github URL of
the MKL (Intel Math Kernel Library) package.  This is because `MKL` is an officially registered Julia package.  Since 
`TE11CylCavity` is not (yet?) registered, you have to install it via its Github repository URL.  
This also means that updates to the 
`TE11CylCavity` package will not be picked up automatically via the convenient `Pkg.update()` function of Julia's package
manager.  Instead, to pull in changes that were made to the repo, you have to remove the package from the shared 
environment and then re-add it as follows:

```Julia
using Pkg: Pkg
Pkg.activate("Cavity", shared=true)
Pkg.rm("TE11CylCavity")
Pkg.add(url="https://github.com/simonp0420/TE11CylCavity")
```

Remember that this only needs to be done if there are updates to the repository that you want to have in your
local Julia environment.
