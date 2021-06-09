# Tensegrity Finite Element Method (TsgFEM)

### **Welcome to **TsgFEM** software!**

#### General Information
 
Our group focuses the research of integrating structure and control design. Based on the tensegrity paradigm, we design tensegrity structures to meet the specified objectives. These objectives can vary from minimizing the mass of the structure to controlling the structure to meet certain performance. This software is intended to study the statics and dynamics of tensegrity systems based on FEM. The authors would like to make this as open-source software to help other researchers who are also interested in this field.

The purpose of this software is to facilitate the statics and dynamics of Tensegrity systems based on FEM.
The software allows to perform modeling, structural design, nonlinear static and dynamic FEM simulation of any tensegrity and truss systems. Contribution of this software is mainly in these three aspects:
#### Modeling: 
1. Model any tensegrity structures by nodal coordinates and the nodes’ connectivity information.
2. Specify the constraints of nodal coordinates (restricted motions of the X-, X-, or/and Z-directions of some certain nodes). 
3. Normally, structure members in symmetric positions have the same force densities, they are allowed to be grouped in the software. 
#### Statics: 
1. Prestress modes and mechanism modes analysis by singular value decomposition of equilibrium matrix.
2. Prestress design and minimal mass calculation considering yielding and buckling constraints.
3. Stiffness and stability analysis.
4. Solve the equilibrium equations with any given external forces considering nonlinear of geometry and material.
5. Simulate the forced motion by given sequencies of some nodes or changing rest length of some members.
#### Dynamics: 
1. Linearized dynamics.
2. Modal analysis, calculate natural frequencies and mode shapes.
3. Nonlinear dynamics simulation of any tensegrity structures or truss systems with linear elastic，multilinear elastic or plastic materials.


The name of this software _TsgFEM_ is suggested to be pronounced as _Tenseg FEM_. TsgFEM allows users to perform accurate dynamics and statics analysis for any tensegrity structures, including:
1. Rigid body dynamics with acceptable errors. This is achieved by setting relatively high stiffness for bars in the dynamics simulation. 
2. Actuate FEM dynamics simulation with elastic or plastic deformations in the presence of various kinds of boundary conditions, such as fix any nodes in any direction, apply static or dynamic external forces (i.e., gravitational force, some specified forces, or arbitrary seismic vibrations). 
3. Accurate modal analysis, including natural frequency and corresponding modes. 
4. An interface towards structural control by a compact state-space form of an accurate linear model. 
5. Statics analysis, which allows one to do form-finding, static load analysis as well as critical buckling analysis. 
 
The software is being used to develop various structures not limited to tensegrity structures, such as torus, cable dome, tower, lander, D-Bar Structures. The software is also applicable to analyze tensile membranes (pure string to string network) and truss structures (pure bar to bar network).


Undergraduate linear algebra, material mechanics/continuum mechanics, finite element method, and some basic knowledge of MATLAB is required to understand the codes well. This software is developed based on:
- 64-bit Windows
- MATLAB 
- MATLAB Optimization Toolbox

Note: Win7/Win10/Mac OS/Linux/Win XP/Win Vista, the software is compatible with a MATLAB version later than 2009a. However, we encourage the user to run the software with the latest MATLAB release if possible. (More information about MATLAB versions can be found here: https://en.wikipedia.org/wiki/MATLAB). Since ‘linprog’ and ‘fmincon’ functions from MATLAB Optimization Toolbox are used in statics calculation, so this toolbox should also be installed (more information can be found here: https://www.mathworks.com/products/optimization.html). 
Since commercial software embedded a lot more material database, for research purposes, one may also use 
- ANSYS

TsgFEM also provide an interface to ANSYS. The automatically generated interface file allows users to perform the FEM simulation in ANSYS automatically and it is compatible with a ANSYS version later than 16.0. However, we also encourage the user to run the software with the latest ANSYS release if possible.


#### LICENSE

    /* This Source Code Form is subject to the terms of the Mozilla Public
     * License, v. 2.0. If a copy of the MPL was not distributed with this
     * file, You can obtain one at http://mozilla.org/MPL/2.0/.
 
---

***<font size=4>Tensegrity Finite Element Method (TsgFEM) folder contains the following parts:</font>***

---

#### setup.m 
If one wants to start using the MOTES software, run the 'setup.m' first.
Open MATLAB and run 'setup.m' file, it will:

- Add all library functions to MATLAB path, 
- Add the Software Verification and Examples,
- Add User Guide,
- Add Videos folder,
- Add JOSS paper.

Note: setup.m must be run every time MATLAB opens, before any other file.

#### JossPaper

This folder contains the journal paper corresponding to the software, including source file documents and reference. The journal paper will provide the background introduction, a summary of our work, applications, references, etc. 

#### Function Library

This folder contains:

All the library functions for tensegrity statics and dynamics analysis are organized in this folder. By following instructions of statics and dynamics analysis from ''User_Guide.pdf", one can perform the analysis.

#### Software Verification and Examples

This folder contains:

1. Statics Examples

Here, we give examples to verify and demonstrate the statics of this software.

2. Dynamics Examples

Here, we give examples to verify and demonstrate the dynamics of this software.

#### User Guide

This folder contains "User_Guide.pdf". The file provides a detailed description, how to use the software, various applications, and how to become a developer as well.

#### Videos
Some interesting tensegrity animation examples are shown in this folder.

---

### Help Desk:

We are open and willing to answer any question. Please state your problem clearly and use the following emails to contact: Muhao Chen: <muhaochen@tamu.edu>, Shuo Ma: <mashuo@zjut.edu.cn>. Thank you!

----

### Join TsgFEM Community and Contribute

#### How to contribute

Feedbacks and contributions are appreciated. Please use the same nomenclature so that everybody could be on the same page.

1. Fork it
2. Submit a pull request OR send emails to help desk.

We will reply to you ASAP.

#### Coding standards

* MATLAB (>= 2009a)
* Function input and output comments
* Use the same Nomenclature as follows

#### Nomenclature

##### Geometry: 
    N: initial node positions
    n: nodal coordinate vector
    C_b: bar connectivity
    C_s: string connectivity
    C: connectivity matrix of the structure
    ne: amount of members
    nn: amount of nodes
    a: a vector containing the index of free nodal coordinate
    b: a vector containing the index of fixed nodal coordinate
    Ia: a matrix locating the free nodal coordinate
    Ib: a matrix locating the fixed nodal coordinate
    Gp: group matrix
    l: length of members
    l0: rest length of members    
##### Statics
    A_1a: equilibrium matrix considering boundary constraints
    V2: null space of equilibrium matrix
    w0: external force
    q_gp: force density vector in group
    q: force density vector
    t_gp: force vector in group 
    t: force vector
    index_b: number of bars
    index_s: number of strings
    A_b: cross sectional area of bars
    A_s: cross sectional area of strings
    A_gp: cross sectional area of all members in group
    A: cross sectional area of all members
    E: Young's modulus of members
    rho: material density of members
    mass: mass vector of members
    substep: substeps in static analysis
    w_t: external force in substeps
    dnb_t: forced motion of fixed nodal coordinate in substeps
    l0_t: rest length of members in substeps
    t_t: force of members in substeps
    n_t: nodal coordinate vector in substeps
     
##### Dynamics
    tf: simulation time duration
    dt: simulation time step
    tspan: time steps in solving dynamics equation
    out_tspan: time steps for output
    lumped: use lumped or consistent mass matrix
    M: mass matrix
    D: damping matrix
    omega: natural frequency
    V_mode: free vibration mode
    dnb_d_t: velocity of forced motion 
    dnb_dd_t: acceleration of forced motion 
    l_t: members length in time steps
    nd_t: velocity of nodal coordinate in time steps

