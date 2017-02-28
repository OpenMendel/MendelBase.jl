# MendelBase

MendelBase includes all the base functions of [OpenMendel](https://openmendel.github.io). It includes useful utilities such as functions to process keywords that specify the data files to use and the analysis options to perform, and functions to read data from external files.

[![](https://img.shields.io/badge/docs-current-blue.svg)](https://OpenMendel.github.io/MendelBase.jl)

## Installation

*Note: Three OpenMendel packages - (1) [SnpArrays](https://openmendel.github.io/SnpArrays.jl/latest/), (2) [Search](https://openmendel.github.io/Search.jl), and (3) [MendelBase](https://openmendel.github.io/MendelBase.jl) must be installed before any Mendel analysis packages will run. It is easiest to install if these three packages are installed in the above order and before any other OpenMendel packages.*


Within Julia, use the package manager to install MendelBase:

    Pkg.clone("https://github.com/OpenMendel/MendelBase.jl.git")

This package supports Julia v0.4 and v0.5.

## Citation

If you use [OpenMendel](https://openmendel.github.io) analysis packages in your research, please cite the following reference in the resulting publications:

*Lange K, Papp JC, Sinsheimer JS, Sripracha R, Zhou H, Sobel EM (2013) Mendel: The Swiss army knife of genetic analysis programs. Bioinformatics 29:1568-1570.*

<!--- ## Contributing
We welcome contributions to this Open Source project. To contribute, follow this procedure ... --->

## Acknowledgments

This project is supported by the National Institutes of Health under NIGMS awards R01GM053275 and R25GM103774 and NHGRI award R01HG006139.
