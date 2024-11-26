Visualization and statistics for brain connectivity data of 3xTgAD mice. 

We used permutation tests, t-tests and  Network Based Statisics (NBS; Zalesky et al., (2010)) to analyze differences in connectivity between different mouse strains :
- Wild type (C57BL/6J)
- Alzheimer's Model (3xTgAD)
- Wild type TSPO KO
- Alzheimer's Model TSPO KO

The goal is to observe the effect of the interaction of Alzheimer's and the deltion of an inflamation's pathway in brain connectivity, mainly between hippocampal regions and sensory cortex. 

The connectivity measures are based on blood flow in the brain and were aquired using functional ultrasound (fUS) imaging.

To run the analysis, one must have the raw data organized as such: 

main_dir|-- clone repo here
        |-- data
                |--code_animaux.xlsx
                |--3xTgAD
                    |--raw matrix 1
                    |--raw matrix n
                |--3xTgAD_TSPO_KO
                    |--raw matrix 1
                    |--raw matrix n
                |--TSPO_KO
                    |--raw matrix 1
                    |--raw matrix n
                |--WT
                    |--raw matrix 1
                    |--raw matrix n