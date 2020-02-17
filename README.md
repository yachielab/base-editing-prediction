# Base editing prediction model used in Sakata, Ishiguro, Mori et al. *bioRxiv* (2019)

## Installation and User Manual
This is a Python script used in Sakata, Ishiguro, Mori et al (2019) to predict frequencies of base editing patterns for a given input sequence using a model trained using amplicon sequencing data obtained for a specific base editing method. Let <img src=images/S_i.png width=14 height=14> be the nucleotide base transition status at *i* bp position relative to the PAM at the target site and <img src=images/P_si.png width=36 height=16>be the probability of  <img src=images/S_i.png width=14 height=14>. A base editor model is prepared as a profile of <img src=images/P_si.png width=36 height=16> and <img src=images/P_sjsi.png width=55 height=18> that can be prepared from amplicon sequence data of different target sites treated with the corresponding base editor (sample codes to generate a base editing model from amplicon data can be found in sample_training_codes/).


In this script, a predicted frequency of a given editing pattern for an input target sequence is calculated by the following formula:

 

<img src=images/equation.png width=250>

 

where
<img src=images/S_mn.png width=18> is a base editing pattern in a window spanning from *m* bp *n* bp relative to the PAM, which can be alternatively represented by a string of transition statuses,<img src=images/S_list.png width=140>,

<img src=images/difinition.png width=250>,



This script also enables to predict frequencies of the all possible base editing patterns for an input target sequence and generate an expected editing spectrum with total base editing frequencies across different positions relative to the PAM.



## Software Dependency

Python 3.7.0 or later



## Installation

1. Donwload the software by

   ``git clone https://github.com/yachielab/base-editing-prediction``

2. Install the necessary Python packages

   ``pip install matplotlib``
   ``pip install numpy``

   

## SYNOPSIS

```
SYNOPSIS
	base-editing-prediction.py [-i] [-o] [-s] [-e] [-m] [-f][-h] [--help]

Options:
-i <String>
	Input target sequence
-o <String>
	Editing outcome sequence
-s <Integer>
	Start position of the input target sequene relative to the PAM¥
-e <Integer>
	End position of the input target sequene relative to the PAM
-m <String>
	File path of base editing model
-f <String>
	Output file path and name. Ex. [file_path]/[output]
-h, --help
	Print the usage of this script
```

## Usage
### Sample codes
#### Example 1: Predicting a frequency of a given editing outcome pattern

````
python base-editing-prediction.py -i ACACACACACACTCTGATCATACGA -o ACACACATGTGCTCTGATCATACGA -e '-1' -m sample_models/TargetACE.csv
````

**Output :** Standard output as follows

````
Model name        : Target-ACE
Target sequence   : ACACACACACACTCTGATCATACGA
Editing outcome   : ACACACATGTGCTCTGATCATACGA
Editing frequency : 0.005603696026973482
````


#### Example 2: Predicting frequencies of all the possible base editing patterns for a given target sequence

````
python base-editing-prediction.py -i ACACACACACACTCTGATCATACGA -s '-25' -m sample_models/TargetACEmax.csv -f [file_path]/[output]
````

**Output :** Generates the following three files in \[file\_path\]

1\. **\[output\]\_allpatterns.csv**

A CSV file showing all of the possible base editing patterns for a given target sequence and their outcome frequencies. The editing patterns are
sorted by their frequencies.


````
#Model name      : Target-ACE
#Target sequence : ACACACACTCTGATCATACGA
#Start position  : -21
#End poition     : -1
Editing outcome,Editing frequency
ACATACGCTCTGATCATACGA,0.017277907935363164
ACATACACTCTGATCATACGA,0.015349088637839415
ACATATGCTCTGATCATACGA,0.014826841469432842
ATATATGCTCTGATCATACGA,0.012313483323590231
ACATATACTCTGATCATACGA,0.009293002375035244
ATATACGCTCTGATCATACGA,0.008889011068741216
ACACACGCTCTGATCATACGA,0.007599118491103579
ATATATACTCTGATCATACGA,0.007163219683861949
ACATGCGCTCTGATCATACGA,0.005779899155419436
ACATGTGCTCTGATCATACGA,0.005660730160931537
...
````

2\. **\[output\]\_spectrum.csv**

A CSV file showing total frequencies of the three possible base
transition patterns in every position across the target sequence

````
#Model name      : Target-ACE
#Target sequence : ACACACACTCTGATCATACGAGGG
#Start position  : -21
#End poition     : 2
Position to the PAM,Target nucleotide,Frequency of A, Frequency of T, Frequency of G, Frequency of C
-21,A,N.A,0,0,0
-20,C,0.0017650751761785623,0.07406627909928758,0.0006266740443025214,N.A.
-19,A,N.A,6.201996639981142e-05,0.0029083418558160917,3.7436383472460367e-06
-18,C,0.004170815198189767,0.15798486306217932,0.007738787007075881,N.A.
-17,A,N.A,0.00011560970648607727,0.05574550877309472,1.3950684141143507e-05
-16,C,0.0019216814569151847,0.09208643414636355,0.0009912734265383296,N.A.
-15,A,N.A,0.0001666896121976657,0.12175455417528482,9.20268535045647e-06
...
````

3\. **\[output\]\_spectrum.pdf**

A pdf file visualizing total frequencies of the three possible base transition patterns in every position across the target sequence.






### File format for base editor model

A base editor model needs to be prepared as a CSV file in the following format

````
#Model Name : Target-ACE
#TP         : Target transition probability
#CTP        : Conditional transition probability
Data type,Conditional transition,Target transition,Probability
TP,,-30:A>A,0
TP,,-30:A>T,1.74142232733626e-05
TP,,-30:A>G,1.1143076691010203e-05
TP,,-30:A>C,0
TP,,-30:T>A,0
TP,,-30:T>T,0
TP,,-30:T>G,1.8732512352396277e-05
TP,,-30:T>C,0.00013427188151391074
TP,,-30:G>A,6.913704645176088e-05
TP,,-30:G>T,9.668251416054742e-05
...
CTP,-30:A>A,-30:A>A,1.0
CTP,-30:A>A,-30:A>T,0
CTP,-30:A>A,-30:A>G,0
CTP,-30:A>A,-30:A>C,0
CTP,-30:A>T,-30:A>A,0
CTP,-30:A>T,-30:A>T,1.0
CTP,-30:A>T,-30:A>G,0
CTP,-30:A>T,-30:A>C,0
CTP,-30:A>G,-30:A>A,0
...
````

**Data type:** TP (transition probability) or CTP (Conditional transition probability)

**Conditional base transition:** {Relative position from the PAM}:{nucleotide transition pattern}. This should be left empty when the data type is TP, or ignored

**Target base transition:** {Relative position from the PAM}:{nucleotide transition pattern}

**Probability:** Probability of the target base transition (TP) or probability of the target base transition given the conditional base transition (CTP)

### Sample models 
Base editing models for xx different base editing methods used in Sakata, Ishiguro, Mori et al. (20xx) are provided in sample\_models/. All of the base editing models were created for a target sequence region from -30 bp to +10 bp to the PAM.

Cytosine base editors (CBEs):
-  ``sample_models/TargetAID.csv``
-   ``sample_models/xxxx.csv``
-   ``sample_models/xxxx.csv``
-   ``sample_models/xxxx.csv``

Adenine base editors (ABEs):
-   ``sample_models/ABE.csv``
-   ``sample_models/xxxx.csv``

Base editor mixes:
-   ``sample_models/TargetAID_plus_ABE.csv``
-   ``sample_models/xxxx.csv``
-   ``sample_models/xxxx.csv``
-   ``sample_models/xxxx.csv``

Dual function base editors
-  ``sample_models/TargetACE.csv``
-  ``sample_models/xxxx.csv``
-  ``sample_models/xxxx.csv``
