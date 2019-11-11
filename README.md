# BPHParkingAnalysis
A repository for analyzing the nanoAOD output of the [BParkingNANO framework](https://github.com/CMSBParking/BParkingNANO) with [RDataFrame](https://root.cern/doc/master/classROOT_1_1RDataFrame.html).

# Run the analyzer
To run the analyzer, inside the folder `RDFanalysis`, do
```
source settings.sh

python analyzer_RDF.py --inList INPUT_LIST.txt --JOBid 0 --outFile /PATH/TO/STORAGE/FOLDER/output.root
```
where `INPUT_LIST.txt` is the txt file containing the paths of the nanoAOD files to be analyzed.
By default, the output will be stored as histograms. Please add `--tree` if you want a (skimmed) tree as output.

The output consists of both KEE and KMuMu quantities.

# Launch batch jobs
To launch batch jobs, inside the folder `scripts`, do
```
python cmsSplit.py --anType scriptAndJOBID --cfg config_runAnalysis.sh --tag YOUR_TAG --listquery -i INPUT_LIST.txt --filesperjob (2, 10,...) --storeArea /PATH/TO/STORAGE/FOLDER

source launch_YOUR_TAG.sh
```
The output will be stored under the path `/PATH/TO/STORAGE/FOLDER/YOUR_TAG`. `config_runAnalysis.sh` contains the command to run the analyzer (add `--tree` to produce skimmed nanoAOD files). 

To produce histograms using the skimmed nanoAOD files, set `--cfg config_histoProducer.sh` after doing
```
g++ -Wall -o histoProducer `root-config --cflags --glibs` -lRooFitCore histoProducer.cpp
```
inside the folder `macro`.
In this case `INPUT_LIST.txt` will be the txt file containing the paths of the skimmed nanoAOD files.
