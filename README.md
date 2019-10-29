# BPHParkingAnalysis
A repository for the nanoAOD analysis with RDataFrame in Python

# Run the analyzer
To run the analyzer, inside the folder RDFanalysis, do
```
source settings.sh

python analyzer_RDF.py --isEleCh 1 --inList input_list.txt --JOBid 0 --outFile /path/output.root
```

# Launch batch jobs
To launch batch jobs, inside the folder scripts, do
```
python cmsSplit.py --anType scriptAndJOBID --cfg config_runAnalyis.sh --tag YOUR_TAG --listquery -i YOUR_INPUT_LIST.txt --filesperjob (2, 10,...) --storeArea /PATH_TO_YOUR_STORAGE_FOLDER

source launch_YOUR_TAG.sh
```
