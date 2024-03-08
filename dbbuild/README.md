`build_HMMs.sh` script automatically builds the HMM profile database from the following data resources.
This is a modified version of the database build [script](https://github.com/tseemann/barrnap/blob/0.9/build/build_HMMs.sh) provided in Barrnap v0.9.
pybarrnap uses only [Rfam](https://rfam.org/) to create the HMM profile database and does not use [SILVA](https://www.arb-silva.de/), which is used in Barrnap v0.9.

```diff
Bacteria (70S)  
        LSU 50S
                5S      RF00001
-               23S     SILVA-LSU-Bac (Barrnap v0.9)
+               23S     RF02541 (pybarrnap)
        SSU 30S
                16S     RF00177

Archaea (70S)   
        LSU 50S
                5S      RF00001
                5.8S    RF00002
-               23S     SILVA-LSU-Arc (Barrnap v0.9)
+               23S     RF02540 (pybarrnap)
        SSU 30S
                16S     RF01959

Eukaryote (80S)   
        LSU 60S
                5S      RF00001
                5.8S    RF00002
-               28S     SILVA-LSU-Euk (Barrnap v0.9)
+               28S     RF02543 (pybarrnap)
        SSU 40S
                18S     RF01960

Metazoan Mito
                12S     RefSeq (MT-RNR1, s-rRNA, rns)
                16S     RefSeq (MT-RNR2, l-rRNA, rnl)      
```
