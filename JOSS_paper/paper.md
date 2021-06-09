---

title: 'TsgFEM: Tensegrity Finite Element Method'


tags:
  - Tensegrity systems
  - Multibody dynamics
  - Flexible structures
  - Prestressable structures
  - Finite element method
  - Linearized tensegrity dynamics
  - Elastic and plastic deformation

authors:
  - name: Shuo Ma
    orcid: 0000-0003-3789-2893
    affiliation: 1

  - name: Muhao Chen^[corresponding author]
    orcid: 0000-0003-1812-6835
    affiliation: 2

  - name: Robert E. Skelton
    orcid: 0000-0001-6503-9115
    affiliation: 2
    
affiliations:
 - name: College of Civil Engineering, Zhejiang University of Technology, Hanzhou, Zhejiang, China
   index: 1

 - name: Department of Aerospace Engineering, Texas A&M University, College Station, Texas, USA
   index: 2

date: 4 June 2021
bibliography: paper.bib
---

# Summary


Tensegrity is a coined word by Buckminister Fuller [@fuller1982synergetics] for the art form created by Ioganson (1921) and Snelson (1948) [@lalvani1996origins] to represent a stable network of compressive members (bars/struts) and tensile members (strings/cables) [@chen2020general]. From the definition, it is straightforward to see that the fundamental property of the tensegrity is that all the bars and strings are axially loaded. Since the bars and strings are best in taking compression and tension and there is no material bending, the structure mass can be greatly reduced. Indeed, biological structures also indicate that tensegrity concepts yield the most efficient structures. For example, the bones, muscles, and elbows of animals and humans are tensegrity models. Wang et al. found that microtubules and microfilaments in the living cells work as compressive and tensile members to change the traction of the cell surfaces [@wang2001mechanical]. Simmons et al. showed that the DNA bundles are consistent with tensegrity prism [@liedl2010self]. After decades of study, many lightweight structures have been redesigned by the tensegrity paradigm. For example, Skelton and de Oliveira proved that T-Bar and D-Bar structures require less mass than a signal rod in taking compressive buckling load [@skelton2009tensegrity]. Ma et al. showed a mass efficient tensegrity cantilever structure subject to yield and buckling constraints [@ma2020design]. Barbarigos et al. designed and analyzed a lightweight pedestrian bridge [@rhode2010designing]. And many space applications are employing tensegrity solutions, for example, lightweight space habitat [@chen2021review], deployable lunar tower for space mining [@chen2020deployable], and planetary landers [@luo2017analysis;@sabelhaus2015system]. Moreover, the many advantages of tensegrity have also attracted researchers to find new ways to design soft robotics, i.e., six-bar tensegrity robot [@booth2020surface;@wang2020first], robotic spine [@sabelhaus2020model], morphing wings [@chen2020design], robotic fish [@chen2019swimming], debris capturing robot [@feng2021design].

The statics and dynamics analysis of tensegrity structures is essential to get insight into the structure properties. Based on the following assumptions, the FEM equations for the statics and dynamics of any tensegrity structures are formulated: 1). All the structural members are axially loaded and are connected by frictionless ball joints. 2). The structural members are allowed to have elastic or plastic deformation. 3). The rotation of the structure member along its longitudinal axis is neglected. 4). Each structural member is homogeneous along its length and of an equal cross-section. Thus, the mass of each structural member is distributed uniformly along its length. 5). A string can never push along its length; tension in the string is substituted to zero. Based on the FEM and Lagrangian methods with nodal vectors as the generalized coordinates given in [@chen2020general;@ma2021tensegrity], we developed TsgFEM (Tensegrity Finite Element Method). 

The name of this software, TsgFEM, is suggested to be pronounced as Tenseg FEM. TsgFEM provides two categories for the analysis of any tensegrity structure: statics and dynamics. For statics analysis, TsgFEM gives the minimum mass of the tensegrity structures by optimizing cross-sectional area in the strings and bars in the absence of external forces (self-equilibrium state) with any given external forces. The software allows solving for the minimum mass subject to buckling and yielding failures of structure members. The other important function of static analysis is that it can calculate the equilibrium configuration of tensegrity structures in any given external force, the forced motion of nodes, and the rest length change of strings. In the static analysis, the nonlinear properties of material and geometry can be considered, and the modified Newton method is used to guarantee the result converges to a stable equilibrium configuration. For the analysis of the dynamics, TsgFEM uses a second-order vector form differential equation to simulate the dynamics of any complexity of the tensegrity structure, including 1). Rigid body dynamics with acceptable errors. This is achieved by setting relatively high stiffness for bars in the dynamics simulation. 2). Accurate FEM dynamics simulation with elastic or plastic deformations in the presence of various kinds of boundary conditions, such as fix any nodes in any direction, apply static or dynamic external forces (i.e., gravitational force, some specified forces, or arbitrary seismic vibrations). 3). Accurate modal analysis, including natural frequency and corresponding modes. 4). An interface towards structural control by a compact state-space form of an accurate linear model.

# Statement of need

For the FEM analysis for tensegrity structures, little research has been conducted. For example, Zogoul et al. conducted a static study for a tensegrity bridge by FEM [@zgoul2012static]. Jensen et al. showed the finite element analysis of tensegrity structures in offshore aquaculture installations by ABAQUS [@jensen2007finite]. Kan et al. formulated the dynamics of clustered tensegrity structures by FEM [@kan2018nonlinear]. However, most of these analyses are either using commercial software or assuming the structure members are elastic ones. Few software packages have been developed for the analysis of tensegrity statics and dynamics. For example, STEDY [@tadiparthi2019stedy] is a package for conducting tensegrity dynamics with rigid bars based on the Lagrangian method. MOTES [@goyal2019motes] is a software for analysis of both statics and dynamics with rigid bars and linear elastic strings. Both of them are rigid body dynamics and developed in non-minimum Cartesian coordinates. However, for many applications, the elastic or plastic deformation of structure members cannot be neglected. It is possible to do the FEM analysis by commercial software, i.e., ANSYS and ABAQUS, but the licenses are expansive, and the modeling and setups require much experience. To this end, we derive a closed-form dynamics equation based on Lagrangian's method with a nodal vector as the generalized coordinate, allowing the user to perform simulations for bars and strings with any given strain-stress properties. The proposed dynamics equation is in a compact form and accurate.


The developed software is capable of dealing with large deformation for structures with elastic or plastic materials. The results of the software are compared with analytical solutions as well as commercial software ANSYS. The accuracy of the software is proved to be very accurate, i.e., in the modal analysis of our examples, results show that the natural frequency errors are about $10^{-14}$. In fact, this software also provides an interface to ANSYS. That is, one can modify the settings in the example source codes based on the simulation needs and run the codes. This software can generate a file that allows users to run the simulation in ANSYS automatically. The software has been used for the design and analysis of various tensegrity structures, i.e., tensegrity lunar tower, tensegrity morphing wings, tensegrity cable domes, tensegrity lander, tensegrity dolphin, where we integrate structure, control design, and signal processing to get the required performance and control law.



# References

