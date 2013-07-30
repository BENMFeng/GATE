-![alt text](https://raw.github.com/BENMFeng/GATE/master/GATE.gif "auto-linked README.md URL")

#GATE
====
Bioinformatics Stargate --
**G**enomics integrated **A**pplications, also for **T**ranscriptomics, **E**pigenetics, Evolution experiential analysis pipeline.
Mainly suitable for generating batch jobs of applications for illumina sequencing platform.

  Copyright (c) 2012 - 2013 - BENM(Binxiao) Feng                        
  All Rights Reserved                                                   
  Send all comments to BENM - BinxiaoFeng\@gmail.com                     
                                                                        
  This program is free software: you can redistribute it and/or modify, it under the terms of the GNU General Public License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version. This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more details. You should have received a copy of the GNU General Public License along with this program.  If not, see <http://www.gnu.org/licenses/>
  
GATE is an integrative PERL pacakge to conduct application for massive parallel bioinformatics analysis. It includes the most popular software which used for different purpose. It is consist of multiple processes, including qc, denovo, predict, aln, exp, diff, as, var, peak, edit, fusion, ncrna, anno, phylogen, stat and plot. Some pipeline with star symbol in head of name is not finished yet, which under development.

User need to setup the software you need, and create an config file for running GATE. GATE would generate a pipeline bash shell script for you to run on Linux server or submit to PBS(qsub/bsub).

Donwload: git clone git://github.com/BENMFeng/GATE.git

Exmaple: bin/GATE --config=config.txt --prefix=test qc aln var

