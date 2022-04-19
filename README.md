# AutoESDv1

## Introduction

AutoESDv1 enables users to perform the precise, automated and high-throughput design of sequence manipulation tasks across species and technical variations, at any genomic locus for all manipulation types with adequate sequence length. The screening-marker-based homologous-recombination system and the overlap-based assembly method were chosen as the loading techniques. Homologous arms and primers required for sequence manipulation, vector DNA assembly and sequencing verification were provided as design results.

## Installation

### python packages

```shell
pip install -r requirements.txt

```

### blast+
```shell
wget https://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/LATEST/ncbi-blast-2.13.0+-x64-linux.tar.gz -O ~/ncbi-blast-2.13.0+-x64-linux.tar.gz

tar -zxvf ncbi-blast-2.13.0+-x64-linux.tar.gz

export PATH=~/ncbi-blast-2.13.0+/bin:$PATH

```

## Usage & Example

### Help message
```shell
python autoesd.py -h 
```

### PSS
This moudle is used for implementing the editing sequences design for plasmid-mediated single/single crossover HR technical variant.
"plasmid_single_single.py" is the design program.
"config.txt" is the parameter setting file.
The input files include the "Target genome" ("Cg_13032.fna"), the "Linear plasmid sequence" ("PSS_cgl_linear_plasmid.txt") and the "Target manipulations" ("PSS_cgl_input_example.csv").
The output files include the detailed information of the designed primers ("Design_results.xlsx"), the tasks without accessible primers and their failure reasons ("Failed_task.xlsx"), the list of primers provided to the primer synthesis company ("Primer_order.xlsx"), and the target sequences that can be mapped to multiple loci in the target genome and may lead to potential off-target events ("Evaluation_result.xlsx")
The series of ".json" files are used for visualization of the results.
Calculation for the example design tasks can be performed by running "cmd.sh" or run the following command in the root path:
```shell
python autoesd.py -t PSS -i ./PSS/test_input/PSS_cgl_input_example.csv  -p ./PSS/test_input/PSS_cgl_linear_plasmid.txt -r ./PSS/test_input/Cg_13032.fna  -c test_input/config.txt -d ./PSS/test_output
```

### FDS
This moudle is used for implementing the editing sequences design for fragment-mediated double/single crossover HR technical variant.
"fragment_double_single.py" is the design program.
"config.txt" is the parameter setting file.
The input files include the "Target genome" ("Ec_MG1655.fna"), the "Screening marker sequence" ("FDS_eco_screening_marker.txt") and the "Target manipulations" ("FDS_eco_input_example.csv").
The output files include the detailed information of the designed primers ("Design_results.xlsx"), the tasks without accessible primers and their failure reasons ("Failed_task.xlsx"), the list of primers provided to the primer synthesis company ("Primer_order.xlsx"), and the target sequences that can be mapped to multiple loci in the target genome and may lead to potential off-target events ("Evaluation_result.xlsx")
The series of ".json" files are used for visualization of the results.
Calculation for the example design tasks can be performed by running "cmd.sh" or run the following command in the root path:
```shell
python autoesd.py -t FDS -i ./FDS/test_input/FDS_eco_input_example.csv  -s ./FDS/test_input/FDS_eco_screening_marker.txt -r ./FDS/test_input/Ec_MG1655.fna  -c ./FDS/test_input/config.txt -d ./FDS/test_output -m "No"
```

### PDS
This moudle is used for implementing the editing sequences design for plasmid-mediated double/single crossover HR technical variant.
"plasmid_double_single.py" is the design program.
"config.txt" is the parameter setting file.
The input files include the "Target genome" ("bsu_168_upp_del.fna"), the "Linear plasmid sequence" ("PDS_bsu_linear_plasmid.txt"), the "Screening marker sequence" ("PDS_bsu_screening_marker.txt") and the "Target manipulations" ("PDS_bsu_input_example.csv").
The output files include the detailed information of the designed primers ("Design_results.xlsx"), the tasks without accessible primers and their failure reasons ("Failed_task.xlsx"), the list of primers provided to the primer synthesis company ("Primer_order.xlsx"), and the target sequences that can be mapped to multiple loci in the target genome and may lead to potential off-target events ("Evaluation_result.xlsx")
The series of ".json" files are used for visualization of the results.
Calculation for the example design tasks can be performed by running "cmd.sh" or run the following command in the root path:
```shell
python autoesd.py -t PDS -i ./PDS/test_input/PDS_bsu_input_example.csv -s ./PDS/test_input/PDS_bsu_screening_marker.txt -r ./PDS/test_input/bsu_168_upp_del.fna -p ./PDS/test_input/PDS_bsu_linear_plasmid.txt -c ./PDS/test_input/config.txt -d ./PDS/test_output -u "AGCTACACGCTGTCTTGCTTC" -m "No"
```

### FDD
This moudle is used for implementing the editing sequences design for fragment-mediated double/double crossover HR technical variant.
"fragment_double_double.py" is the design program.
"config.txt" is the parameter setting file.
The input files include the "Target genome" ("Ec_MG1655.fna"), the "Screening marker sequence" ("FDD_eco_screening_marker.txt") and the "Target manipulations" ("FDD_eco_input_example.csv").
The output files include the detailed information of the designed primers ("Design_results.xlsx"), the tasks without accessible primers and their failure reasons ("Failed_task.xlsx"), the list of primers provided to the primer synthesis company ("Primer_order.xlsx"), and the target sequences that can be mapped to multiple loci in the target genome and may lead to potential off-target events ("Evaluation_result.xlsx")
The series of ".json" files are used for visualization of the results.
Calculation for the example design tasks can be performed by running "cmd.sh" or run the following command in the root path:
```shell
python autoesd.py -t FDD -i ./FDD/test_input/FDD_eco_input_example.csv  -s ./FDD/test_input/FDD_eco_screening_marker.txt -r ./FDD/test_input/Ec_MG1655.fna  -c ./FDD/test_input/config.txt -d ./FDD/test_output  -u "CCGGAAACTCCGCTGGGCGA" -m "No"
```

### ODD
This moudle is used for implementing the editing sequences design for oligonucleotide-mediated double/double crossover HR technical variant.
"oligonucleotide_double_double.py" is the design program.
"config.txt" is the parameter setting file.
The input files include the "Target genome" ("Ec_MG1655.fna"), the "Screening marker sequence" ("ODD_eco_screening_marker.txt") and the "Target manipulations" ("ODD_eco_input_example.csv").
The output files include the detailed information of the designed primers ("Design_results.xlsx"), the tasks without accessible primers and their failure reasons ("Failed_task.xlsx"), the list of primers provided to the primer synthesis company ("Primer_order.xlsx"), and the target sequences that can be mapped to multiple loci in the target genome and may lead to potential off-target events ("Evaluation_result.xlsx")
The series of ".json" files are used for visualization of the results.
Calculation for the example design tasks can be performed by running "cmd.sh" or run the following command in the root path:
```shell
python autoesd.py -t ODD -i ./ODD/test_input/ODD_eco_input_example.csv  -s ./ODD/test_input/ODD_eco_screening_marker.txt -r ./ODD/test_input/Ec_MG1655.fna  -c ./ODD/test_input/config.txt -d ./ODD/test_output -u "TGGGCGAGCCGAAAAACAAATA" -m "No"
```

