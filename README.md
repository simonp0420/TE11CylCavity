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

```Julia
using Pkg: Pkg
Pkg.add(url="https://github.com/simonp0420/TE11CylCavity")
Pkg.add("MKL") # This line only needed if you are running on an Intel CPU
Pkg.precompile()
```

These lines install the packages `TE11CylCavity` and `MKL` into the Julia default package environment. There will be
quite a bit of output after the second line.  You can ignore any warnings regarding packages that could be updated.

## Running Scripts
Once you create a script named, say "myscript.jl" containing Julia code to execute, you can run it via a 
function call at the Julia prompt like this:

```Julia
include("myscript")
```

The `include` function reads the file and executes the lines as if they were typed into the Julia REPL 
(Read-Evaluate-Print-Loop, i.e. the Julia interactive command line environment, similar to Matlab's or Python's.)
The above assumes that you started Julia in the same directory as your script, or that you used the Julia `cwd` function
to move Julia into this directory.  If your script is in a different directory you can use the Julia `joinpath`
function to specify the full path to the file.  For example, to run a file "myscript.jl" located in
"C:\home\simonp\julia\scripts" you could use the following:

```Julia
include(joinpath("C:\\", "home", "simonp", "julia", "scripts", "myscript.jl"))
```

The double backslash above is needed because in strings Julia treats a single backslash character as an escape character.

To run one of the scripts included in the "scripts" directory of your installed `TE11CylCavity` package, use, e.g.

```Julia
include(pkgdir(TE11CylCavity, "scripts", "test_findcase6fresQ.jl"))
```

Prior to doing this you must have entered a `using TE11CylCavity` or `using TE11CylCavity: TE11CylCavity` command
in order for `TE11CylCavity` to be defined. The `pkgdir` function used above locates the file "test_findcase4fresQ.jl"
in the subdirectory "scripts" of the installation directory of the `TE11CylCavity` package.  


When you are ready to exit your Julia session, you can either type Ctrl-D or `exit()` at the Julia prompt, to take you
back to the OS prompt.

### More details on updating `TE11CylCavity`

Recall that when we added `MKL` to the default package environment, we did not have to specify the full Github URL of
the MKL (Intel Math Kernel Library) package.  This is because `MKL` is an officially registered Julia package.  Since 
`TE11CylCavity` is not (yet?) registered, you have to install it via its Github repository URL.  
This also means that updates to the `TE11CylCavity` package will not be picked up automatically via the convenient
`Pkg.update()` function of Julia's package manager.  Instead, to pull in changes that were made to the repo, you have
to remove the package from the environment and then re-add it as follows:

```Julia
using Pkg: Pkg
Pkg.rm("TE11CylCavity")
Pkg.add(url="https://github.com/simonp0420/TE11CylCavity")
```

Remember that this only needs to be done if there are updates to the repository that you want to have in your
local Julia environment.
