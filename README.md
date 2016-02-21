# SmartSVA
SmartSVA introduces an improved Surrogate Variable Analysis algorithm that
automatically captures salient features from data in the presence of confounding factors.
Comparing to the popular SVA algorithm, SmartSVA works 10 times faster.

Installing Dependencies
=======================

Smart SVA requires "sva" and "isva" packages. Start R and install these packages 
as follows:
* Installing "sva" (try http:// if https:// URLs are not supported)

  > source("https://bioconductor.org/biocLite.R")
  
  > biocLite("sva")
  
* Installing "isva"

  > install.packages("isva")
  
You may need to install / update other required R packages for "sva" and 
"isva" to function properly.

You may obtain the lastest revision of SmartSVA by cloning this repository:

    > git clone --recursive https://github.com/ehsanbehnam/SmartSVA/

Usage
=====

For a simple usage example, refer to the sample example given with the
source code.

Contact authors
===============

Jun Chen
Chen.Jun2@mayo.edu

Ehsan Behnam
behnamgh@usc.edu

Copyright and License Information
=================================
  
This program is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.
  
This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.
  
You should have received a copy of the GNU General Public License
along with this program.  If not, see <http://www.gnu.org/licenses/>.
