ma-shm
==============================

Parallel time domain spectral element method utilizing flat shell elements for wave propagation modelling
Dependencies: Matlab Parallel Computing Toolbox, gmsh

Notes: 
- gmsh binaries should be copied to \bin\external\gmsh\
- config_matlab.m should be run first

Project Organization
--------------------

    .
    ├── AUTHORS.md
    ├── LICENSE
    ├── README.md
    ├── bin
    │   └── external
    │       └── gmsh
    ├── config
    ├── data
    │   ├── external
    │   │   ├── exp
    │   │   └── num
    │   ├── interim
    │   │   ├── exp
    │   │   └── num
    │   ├── processed
    │   │   ├── exp
    │   │   └── num
    │   └── raw
    │       ├── exp
    │       └── num
    ├── docs
    │   ├── manuals
    │   ├── proposal
    │   └── published_papers
    ├── notebooks
    ├── reports
    │   ├── conference_papers
    │   ├── figures
    │   │   ├── corel_draw
    │   │   └── grapher
    │   ├── journal_papers
    │   ├── presentations
    │   └── project_reports
    └── src
        ├── data_processing
        ├── external
        ├── models
        │   └── common
        ├── tools
        └── visualization
