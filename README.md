#EAGLE#

<br>
<!--<p align="center"><img src="/logo.png" 
alt="EAGLE" width="350" height="100" border="0" /></p>-->
<br>

A tool for efficient compression of FASTQ files.

## INSTALLATION ##

Simply run the following instructions at a Linux terminal:

<pre>
wget https://github.com/pratas/hawk/archive/master.zip
unzip master.zip
cd hawk-master
make cleanall 
make
</pre>

## EXECUTION

### Run Hawk

Run Hawk using:

<pre>
./Hawk -l 7 file.fastq
</pre>

## PARAMETERS

To see the possible options type
<pre>
./Hawk
</pre>
or
<pre>
./Hawk -h
</pre>
<!--
These will print the following options:
<pre>
<p>
Usage: Eagle &#60OPTIONS&#62 ... -r [FILE]  [FILE]:&#60...&#62

  -v                       verbose mode             
  -c  &#60ctx&#62                context size model       
  -i                       use inversions           
  -ea &#60pts&#62                enlarge absent           
  -en &#60pts&#62                enlarge N's              
  -s  &#60sub&#62                sub-sample               
  -o  &#60oFile&#62              output map file          
                                                    
  -r  [rFile]              reference file (database)
                                                    
  [tFile1]:&#60tFile2&#62:&#60...&#62  target file(s)</p>         
</pre>
-->

## CITATION ##

On using this software/method please cite:

Paper currently submitted.
DOI: doi-to-appear

## ISSUES ##

For any issue let us know at [issues link](https://github.com/pratas/hawk/issues).

## LICENSE ##

GPL v2.

For more information:
<pre>http://www.gnu.org/licenses/gpl-2.0.html</pre>


