# Matlab module of exploration seismology

## SeismicLab

- SeismicLab is a MATLAB seismic data processing package. It was developed by the [***Signal Analysis and Imaging Group***](https://saig.physics.ualberta.ca/) of the University of Alberta, Canada.
- The package can be used to process small seismic data sets and, it is mainly intended for research and teaching purposes. Scripts to read and write SU data (the SEGY flavor used by [Seismic Un*x](www.cwp.mines.edu/cwpcodes/index.html)) are provided. A particular feature of this package is that the SEGY headers are loaded into a structure that can be easily assessed by MATLAB.
  - [GNU General Public License](http://seismic-lab.physics.ualberta.ca/gpl.html)
  - Download SeismicLab [ [SeismicLab.tar.gz](http://seismic-lab.physics.ualberta.ca/SeismicLab.tar.gz) | [SeismicLab.tar.Z](http://seismic-lab.physics.ualberta.ca/SeismicLab.tar.Z) | [SeiemicLab.tar](http://seismic-lab.physics.ualberta.ca/SeismicLab.tar)]
  - Add the script `setpath.m` in your working directory
    - download [`setpath.m`](http://seismic-lab.physics.ualberta.ca/setpath.m)
    - see [`setpath.m`](http://seismic-lab.physics.ualberta.ca/setpath.html)
    - try this one for windows [`setpath_Windows.m`](http://seismic-lab.physics.ualberta.ca/setpath_Windows.m)
  - Edit `setpath.m` to reflect the path to SeismicLab in your system
  - Start Matlab, run setpath
  - Try the [demos](http://seismic-lab.physics.ualberta.ca/help.html#A11)
- Users can download the relevant information through the website http://seismic-lab.physics.ualberta.ca/.

## CREWES

### Introduction

- **CREWES MATLABT Software Library** (CMSL) is a software package developed by Gary F. Margrave for teaching and research of exploration seismology.
- An older, less complete, free version is provided [here](https://www.crewes.org/ResearchLinks/FreeSoftware/).
- The software package accompanies the textbook [***Numerical Methods of Exploration Seismology: With Algorithms in MATLAB***](https://www.cambridge.org/core/books/numerical-methods-of-exploration-seismology/53A21CAD45D4047D117191E6BF4408E2) (NMES) by Gary F. Margrave and Michael P. Lamoureux (Cambridge University Press, 2019).

### CREWES Toolbox Version: 2104

- [CREWES Matlab Toolbox (ZIP)](https://www.crewes.org/ResearchLinks/FreeSoftware/crewes_educational.zip) 156.6 MB
- [Sample data to accompany Methods of Seismic Data Processing (ZIP)](https://www.crewes.org/ResearchLinks/FreeSoftware/NMESdata.zip) 13.65 MB
- [Guide to the CREWES Matlab toolbox (PDF)](https://www.crewes.org/ResearchLinks/FreeSoftware/NumMeth.pdf)  6.98 MB
- [Introductory seismic data processing course (PDF)](https://www.crewes.org/ResearchLinks/FreeSoftware/Methods_of_Seismic_Data_Processing.pdf)  88.09 MB

### Installation instructions

- Extract the contents of crewes_educational.zip to:
  - **Microsoft Windows:** ```%USERPROFILE%\Documents\MATLAB\crewes```
  - **Linux:** ```$HOME/Documents/MATLAB/crewes```
- Then start Matlab and look for ```Set Path``` on the Home Ribbon.
- In the window that appears, click on the button ```Add with subfolders...``` .  Locate the folder called ```crewes``` -- Click on it once to highlight it, then click ```Select Folder```.  Save the new path by clicking ```Save```. 

## SegyMAT
- SegyMAT is a set of m-files for reading and writing SEG-Y files from [MATLAB](http://mathworks.com) and [Octave](https://www.gnu.org/software/octave/), that aims to
  - completely support SEG-Y revision 0 and 1;
  - be easy to use in other projects;
  - be a Swiss Army knife dealing with the SEGY-Y format in MATLAB/Octave.

- SegyMAT is not lightning fast. SegyMAT makes heavy use of ‘structures’. Unfortunately structures are not very effective in terms of speed in MATLAB. (Or they have not been implemented very effectively in SegyMAT). However structures make the implementation and maintenance easier, and the code (hopefully) easy to read. That said, some effort has been made to optimize SegyMAT for speed.
- The latest **stable** version of SegyMAT is available from [Sourceforge](https://sourceforge.net/projects/segymat/).
- The current **development** version of SegyMAT is available from [Github](https://github.com/cultpenguin/segymat).
- The SegyMAT's documentation is available from https://segymat.readthedocs.io/en/latest/.
