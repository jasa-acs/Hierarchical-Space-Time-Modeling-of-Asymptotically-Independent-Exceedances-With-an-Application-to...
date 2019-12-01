# Hierarchical Space-Time Modeling of Asymptotically Independent Exceedances With an Application to Precipitation Data

# Author Contribution Checklist

## Data 

### Abstract

The data set consists of hourly observations at 50 rainfall stations for the years 1993 to 2014 in the mountainous areas of the Cevennes and the Alps (France). Extreme rainfall events usually occur during fall season; therefore, our application considers only observations from the September to November months. Elevations at gauge sites range from 140 to 1000 meters. Data have been extracted from the database of the French Met Office, Météo France.

### Availability

This data from Météo France is not freely available to the general public, but access can be requested for academic purposes. Data can be downloaded from Météo France's web interface: https://publitheque.meteo.fr/okapi/accueil/okapiWebPubli/index.jsp)

Restrictions: Along with the submitted manuscript, we provide the data used in the paper for checking the results during the review process, but it will not be made publicly available and may not be redistributed by the editors and reviewers. Météo France grants the data for research purposes upon request, and specifically within the framework of a cooperation agreement with HydroSciences Montpellier / Montpellier University, France.

### Description

Permissions: Data have been provided to us by Météo France within the framework of a cooperation agreement with HydroSciences Montpellier / Montpellier University, France. We are grateful to Julie Carreau (HydroSciences Montpellier) for help with extracting the dataset from the Météo France database.

Licensing information: Details on the terms and conditions for use of the data base can be found on the data provider web portal. We point out that data from the above web interface are freely available for non-commercial research purposes upon request, see the explanations here (in French):
https://donneespubliques.meteofrance.fr/?fond=faq&id_dossier=5#contenu_38

Link to data: Data can be downloaded from the following web interface, although with the above explained access restrictions:
https://publitheque.meteo.fr/okapi/accueil/okapiWebPubli/index.jsp)

To obtain the stations used in the application study in our manuscript, the following identifiers of measurement sites must be requested:

 7011004    7025001    7032003    7117001    7131001  7159001    7172002    7228001    7286004    7330004    7334003    26002003   26047001   26050001   26115001 26124001   26176001   26177001   26202001   26292002   26298001   26313001   26361001   30068001   30087002 30132003   30273001   30299001   34008001   34028003   34030002   34186001   34209002   34217001   34317001 34344001   43096001   43150001   43234005   43268005   48027003   48069001   48095004   48171001   84036001 84082001   84085003   84094001   84107002   84150001 


## Code

### Abstract

Provided code allows to estimate the hierarchical space-time models for exceedances based on gamma process convolutions and Gaussian copula models as proposed in our manuscript.

To efficiently calculate the ellipse intersection area, we follow

Hughes, G. B., and Chraibi, M. (2012), *Calculating ellipse overlap
areas*, Computing and Visualization in Science, 15, 291--301.

and we use selected parts of the code in http://github.com/chraibi/EEOver

The part of the code used for computing the bivariate normal distribution function is taken from

http://www.math.wsu.edu/faculty/genz/software/fort77/mvtdstpack.f

(its author is Alan Genz, Department of Mathematics, Washington State University)

### Description 

How delivered: collection of R and C scripts

Licensing information: GNU Lesser General Public License (LGPL) version 3

(part of the external libraries are published under this licence, such that we do not use the less
restrictive GPL licence)

Link to code/repository: github.com/cgaetan

Version information 1.0

### Optional Information

Hardware requirements: workstation running a Linux OS

Supporting software requirements : The GNU Scientific Library (GSL), (https://www.gnu.org/software/gsl/)


## Instructions for use 

Code has been tested on a workstation running a Linux OS

1\) Install the [GSL](https://www.gnu.org/software/gsl/)  library. 

For instance under a Linux Debian
distribution, open a terminal and type

sudo apt-get install libgsl0ldbl libgsl-dev

For Mac users, please follow the instructions [here](http://macappstore.org/gsl/) 

2\) Compile the source files and then link all specified object files into a shared object

For instance under  Linux OS open a terminal and type

R CMD SHLIB spt-gamma.c zsolve_quartic.c solvers-2.c Roots3And4.c
-L/usr/lib/x86\_64-linux-gnu -lgsl -lgslcblas

3\) file simulation-gamma.R contains the code for a simulated example. To run this 

R CMD BATCH simulation-gamma.R
