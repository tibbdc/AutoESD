# AutoESDv1

## Introduction



## Installation

### python packages

```shell
pip install -r requirements.txt

```

### blast+
```shell
wget https://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/LATEST/ncbi-blast-2.13.0+-x64-linux.tar.gz -O ~/ncbi-blast-2.13.0+-x64-linux.tar.gz

tar zxvf ncbi-blast-2.13.0+-x64-linux.tar.gz

export PATH=~/ncbi-blast-2.13.0+/bin:$PATH

```

## Usage & Example

```shell
# help message
python autoesd.py -h 

# PSS
python autoesd.py -t PSS -i ./PSS/test_input/plasmid_single_single_input_test.csv  -p ./PSS/test_input/plasmid_seq.txt -r ./PSS/test_input/ASM1132v1.fna  -c test_input/config.txt -d ./PSS/test_output

# FDS
python autoesd.py -t FDS -i ./FDS/test_input/fragment_double_single_input_test.csv  -s ./FDS/test_input/screening_marker_seq.txt -r ./FDS/test_input/ASM1132v1.fna  -c ./FDS/test_input/config.txt -d ./FDS/test_output -m "No"

# PDS
python autoesd.py -t PDS -i ./PDS/test_input/plasmid_double_single_input_test.csv -s ./PDS/test_input/screening_marker_seq.txt -r ./PDS/test_input/ASM1132v1.fna -p ./PDS/test_input/plasmid_seq.txt -c ./PDS/test_input/config.txt -d ./PDS/test_output -u "AGTTGCTGTGGCGGAAAGCC" -m "No"

# FDD
python autoesd.py -t FDD -i ./FDD/test_input/fragment_double_double_input_test.csv  -s ./FDD/test_input/screening_marker_seq.txt -r ./FDD/test_input/ASM1132v1.fna  -c ./FDD/test_input/config.txt -d ./FDD/test_output  -u "AGTTGCTGTGGCGGAAAGCC" -m "No"

# ODD
python autoesd.py -t ODD -i ./ODD/test_input/oligonucleotide_double_double_input_test.csv  -s ./ODD/test_input/screening_marker_seq.txt -r ./ODD/test_input/ASM1132v1.fna  -c ./ODD/test_input/config.txt -d ./ODD/test_output -u "AGTTGCTGTGGCGGAAAGCC" -m "No"

```

