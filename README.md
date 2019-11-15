# BPHParkingAnalysis
A repository for analyzing the nanoAOD output of the [BParkingNANO framework](https://github.com/CMSBParking/BParkingNANO) with [RDataFrame](https://root.cern/doc/master/classROOT_1_1RDataFrame.html).

## Getting started
Recipe for setting up the analysis
```
cmsrel CMSSW_10_2_15
cd CMSSW_10_2_15/src
git clone https://github.com/VincenzoCa/BPHParkingAnalysis.git
```

## Run the analyzer
To run the analyzer, inside the folder `RDFanalysis`, do
```
source settings.sh

python analyzer_RDF.py --inList INPUT_LIST.txt --JOBid 0 --outFile /PATH/TO/STORAGE/FOLDER/output.root
```
where `INPUT_LIST.txt` is the txt file containing the paths of the nanoAOD files to be analyzed.
By default, the output will be stored as histograms. Please add `--tree` if you want a skimmed tree as output.

The output consists of both KEE and KMuMu quantities.

There is also a `C++` version of the analyzer, `analyzer_RDF.cpp`. To run `analyzer_RDF.cpp`, inside the folder `RDFanalysis`, do
```
source settings.sh

g++ -Wall -o analyzer_RDF `root-config --cflags --glibs ` analyzer_RDF.cpp
```
before doing
```
./analyzer_RDF --inList INPUT_LIST.txt --JOBid 0 --outFile /PATH/TO/STORAGE/FOLDER/output.root
```
Please add `--skTree 1` if you want a skimmed tree as output.

## Launch batch jobs
To launch batch jobs with `analyzer_RDF.py`, inside the folder `scripts`, do
```
python cmsSplit.py --anType scriptAndJOBID --cfg config_runAnalysis.sh --tag YOUR_TAG --listquery -i INPUT_LIST.txt --filesperjob (2, 10,...) --storeArea /PATH/TO/STORAGE/FOLDER

source launch_YOUR_TAG.sh
```
The output will be stored under the path `/PATH/TO/STORAGE/FOLDER/YOUR_TAG`. `config_runAnalysis.sh` contains the command to run the analyzer (add `--tree` to produce skimmed nanoAOD files). Set `--cfg config_runAnalysis_cpp.sh` if you want to launch jobs with `analyzer_RDF.cpp` (add `--skTree 1` to produce skimmed trees).

If you want to launch batch jobs with `macro/histoProducer_RDF.py`, to produce histograms from the skimmed nanoAOD files, set `--cfg config_histoProducer_RDF.sh` instead. In this case `INPUT_LIST.txt` will be the txt file containing the paths of the skimmed nanoAOD files.
