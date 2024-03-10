- `build_HMMs.sh`  
  Rfam HMM profile database build script for pyhmmer.nhmmer rRNA search.

- `build_CMs.sh`  
  Rfam CM profile database build script for cmscan rRNA search.

**Rfam resource accession for rRNA search**:

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

Eukarya (80S)   
        LSU 60S
                5S      RF00001
                5.8S    RF00002
-               28S     SILVA-LSU-Euk (Barrnap v0.9)
+               28S     RF02543 (pybarrnap)
        SSU 40S
                18S     RF01960

-Metazoan Mito (Barrnap v0.9)
-               12S     RefSeq (MT-RNR1, s-rRNA, rns)
-               16S     RefSeq (MT-RNR2, l-rRNA, rnl)      
```

> [!NOTE]
>
> 1. pybarrnap no longer uses SILVA database.  
> 2. pybarrnap no longer includes metazoan mitochondria resources.  
