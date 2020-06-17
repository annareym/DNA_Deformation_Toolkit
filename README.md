
# DNA Mapipulations Toolkit
DNA twisting and stretching manipulation tools that can be used as [PLUMED Collective Variables](https://www.plumed.org/doc-v2.6/user-doc/html/colvarintro.html). The tools are implemented within [PLUMED](https://www.plumed.org/) free energy library environment and can be used to perform all-atom [steered MD](https://en.wikipedia.org/wiki/Molecular_dynamics#Steered_molecular_dynamics_(SMD)) simulations of DNA twisting and stretching transitions.

This version allows applying twist and stretch restraints to DNA fragments of arbitrary length and curvature, as opposed to [the previous version of DNA twisting tool](https://github.com/annareym/PLUMED_DNA-Twist).

## Authors

Alexey Voronov & Anna Reymer, [Reymer Lab](https://cmb.gu.se/english/about_us/staff?userId=xreyan), Department of Chemistry & Molecular Biology, University of Gothenburg, Sweden.


## Code dependencies

* [Eigen 3.3.7](http://eigen.tuxfamily.org/)
* [Autodiff 0.5.10](https://github.com/autodiff/autodiff)
* [Catch2](https://github.com/catchorg/Catch2) (needed only for unit tests)
* [Plumed 2.6.0](https://www.plumed.org/)


## How to use
DNA twisting and stretching tools work analogously to any collective variable (colvar) implemented in PLUMED. TWIST2 and STRETCH colvars monitor or control the value of total twist and total rise, correspondingly, between any chosen base pair levels _i_ and _j_ in a DNA fragment. As an imput to the colvars, namely plumed.dat file, provide 3\*2\*N atom numbers from DNA bases that will be restrained: (N restrained base levels \* 2 DNA strands \* 3 atoms per base: C1',N1/N9 and C6/C8 depending whether purines or pyrimidines, correspondignly), force constant (_KAPPA_), and desired value (_AT_) of total twist (in degrees) or stretch (in nm). To push a system into the desired conformation, an energy penalty will be added to the potential energy functional: E<sub>Deform</sub>=0.5\*k\*(x<sub>0</sub>-x)<sup>2</sup>. If you want to just monitor the total twist or total stretch values while running MD, provide only 3\*2\*N atom numbers.

#### Example of plumed.dat input file for TWIST2:
```
tw: TWIST2 ATOMS=72,74,87,991,993,1005,104,106,116,958,960,974,134,136,150,928,930,940,167,169,181,896,898,911,199,201
,214,864,866,878,231,233,243,831,833,847,261,263,277,801,803,813,294,296,308,769,771,784,326,328,341,737,739,751,358,3
60,370,704,706,720,388,390,404,674,676,686,421,423,435,642,644,657,453,455,468,610,612,624
tw_r: RESTRAINT ARG=tw KAPPA=0.15 AT=418.8
PRINT STRIDE=500 ARG=tw,tw_r.bias FILE=twist2_out.txt
```

Provided atom numbers should respect the following order: first 3 atom numbers represent base level i on Watson DNA strand (5'->3' DNA direction), second 3 atom numbers represent the complementary base level _i_ on Crick DNA strand (3'->5' DNA direction), the for base level _i+1_, _i+2_,..,_j-1_,_j_. 

Plumed.dat file will look analogously for STRETCH colvar.

## Acknowledgements
This work was supported by Swedish Foundation for Strategic Research SSF Grant [ITM17-0431](https://strategiska.se/en/research/ongoing-research/instrument-technique-and-method-development-2017/project/9774/) to Dr. Anna Reymer.

## License

GNU LGPLv3.
