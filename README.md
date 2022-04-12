# LiGA-helper
LiGA helper scripts

# How to use
1. Run the MatLab scripts to generate the `AllUniqueFiles` directory
2. Run the python script `liga_helper` - Python3.x (7, 8 ...) would work for you.
```
 python3 liga_helper.py --all_uniq_dir AllUniqueFiles --barcode_file BarcodeSampleMap.csv --sdb_map_file PhageBarcodeModificationMap.csv --final_output_dir final_output_dir
```
3. Adjust the Rmarkdown file for sample names and meta data and run it.

## How to use liga_helper 
```
usage: liga_helper.py [-h] --all_uniq_dir ALL_UNIQ_DIR --barcode_file BARCODE_FILE --sdb_map_file SDB_MAP_FILE --dna_count_table_output_dir DNA_COUNT_TABLE_OUTPUT_DIR

Python script for generating DNA count tables.

optional arguments:
  -h, --help            show this help message and exit
  --all_uniq_dir ALL_UNIQ_DIR
                        AllUniqueFiles dir
  --barcode_file BARCODE_FILE
                        BarcodeSampleMap.csv
  --sdb_map_file SDB_MAP_FILE
                        PhageBarcodeModificationMap.csv
  --dna_count_table_output_dir DNA_COUNT_TABLE_OUTPUT_DIR
                        DNA count table output directory
```
