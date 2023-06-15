# PD MAMA Part 2: Meta-Analysis

`GP2 ❤️ Open Science 😍`

- **Module:** Multi-Ancestry Meta-Analysis Pipeline 
- **Authors:** Jonggeol Jeffrey Kim, Dan Vitale
- **Estimated Computation and Runtime:**
    - **Estimated Specifications:** 4 CPUs/28 GB
    - **Estimated Runtime:** < 5 minutes
- **Date Last Updated:** 14-MAR-2023
    - **Update Description:** Re-ran analysis to get estimated runtime

---

### Quick Description: 

Performing a multi-ancestry GWAS meta-analysis using the summary statistics of the largest Parkinson's disease GWASs in European, Asian, Latin American, and African ancestries followed by fine-mapping the loci using MR-MEGA and PESCA. To combine ancestrally diverse association studies to increase power to detect novel PD-associated loci and improve fine-mapping resolution.

### Background/Motivation:

Genome-wide association studies (GWASs) are a popular, cost-effective method to identify and test common genetic variants across the genomes of many individuals to identify genotype-phenotype associations. GWASs have been useful to begin explaining the heritability of common quantitative traits and better understanding the genetic architecture of complex diseases. Genetic variations found significantly in individuals with the disease compared to those without, those variants are then "associated" with disease.

Traditionally, due to lack of data collection in other ethnicities, GWASs have been conducted in European populations to focus on limiting population stratification effects on studies. Multi-ancestry studies aim to leverage the natural differences in linkage disequilibrium across diverse populations, the knowledge known about one population, and assess if those variants replicate across populations populations while increasing statistical power.

Additionally, multi-ancestry GWASs for the identification of SNPs associated with disease risk in mixed populations, allow for fine-mapping of functional variants, and prioritization of candidate genes.

Studies involved:

* European: **Identification of novel risk loci, causal insights, and heritable risk for Parkinson's disease: a meta-analysis of genome-wide association studies.** *Nalls et al. 2019* https://doi.org/10.1016/S1474-4422(19)30320-5
* Finnish European: **FinnGen Release 4, Parkinson's Disease (more controls excluded) G6_PARKINSON_EXMORE**
* Latin American: **Characterizing the Genetic Architecture of Parkinson's Disease in Latinos.** *Loesch et al. 2021* https://doi.org/10.1002/ana.26153
* Asian: **Identification of Risk Loci for Parkinson Disease in Asians and Comparison of Risk Between Asians and Europeans: A Genome-Wide Association Study.** *Foo et al. 2020* https://doi.org/10.1001/jamaneurol.2020.0428

23andMe has also contributed datasets from East Asian, Latin American, and African ancestries. In total, we have 7 different datasets from 4 different ancestries.

### Workflow Summary:
1. Meta-analysis
2. Fine-mapping `<- You are here!`
3. Functional mapping and annotation

The steps are as follows:
1. Read FUMA annotated data
2. subset to variants identified in the locus region defined by 1000 Genome reference panel
3. generate credible set until the posterior probability is greater than 95
---

### Table of Content:

#### [0. Getting started](#0)
#### [1. Read input data](#1)
#### [2. Fine-map using MR-MEGA Bayes Factor](#2)
#### [3. Add annotations to the credible sets](#3)

# 0. Getting Started

Here we will set up the workspace with python packages and downloading functionally annotated data.


```python
import pandas as pd
import os
import subprocess
import sys
import numpy as np
from shutil import copyfile
from scipy.stats import chi2
# import plotly.express as px
```


```bash
%%bash
# download FUMA results
mkdir -p data/FUMA/MRMEGA
gsutil cp gs://fc-12d0d22f-3861-4641-8236-1dbe3a3c4e9d/results/FUMA/PD_MAMA_MRMEGA_PC3_FUMA.zip data/FUMA/

unzip data/FUMA/PD_MAMA_MRMEGA_PC3_FUMA.zip -d data/FUMA/MRMEGA
# download MR-MEGA results
mkdir output/
GS_MRMEGA_RESULT='gs://fc-12d0d22f-3861-4641-8236-1dbe3a3c4e9d/results/meta_analysis/PD_TEMA.Joint_noSAS.all_sep.PC3.MR-MEGA.tsv.gz'
gsutil cp $GS_MRMEGA_RESULT output/
```

# 1. Read input data

We will read the MR-MEGA results and the FUMA annotated data.


```python
loci = pd.read_csv('data/FUMA/MRMEGA/GenomicRiskLoci.txt', delim_whitespace=True)
sumstat = pd.read_csv(
    'output/PD_MAMA.MR-MEGA.tsv.gz',
    delim_whitespace=True,
    engine='c'
)
loci.head()
```




<div>
<style scoped>
    .dataframe tbody tr th:only-of-type {
        vertical-align: middle;
    }

    .dataframe tbody tr th {
        vertical-align: top;
    }

    .dataframe thead th {
        text-align: right;
    }
</style>
<table border="1" class="dataframe">
  <thead>
    <tr style="text-align: right;">
      <th></th>
      <th>GenomicLocus</th>
      <th>uniqID</th>
      <th>rsID</th>
      <th>chr</th>
      <th>pos</th>
      <th>p</th>
      <th>start</th>
      <th>end</th>
      <th>nSNPs</th>
      <th>nGWASSNPs</th>
      <th>nIndSigSNPs</th>
      <th>IndSigSNPs</th>
      <th>nLeadSNPs</th>
      <th>LeadSNPs</th>
    </tr>
  </thead>
  <tbody>
    <tr>
      <th>0</th>
      <td>1</td>
      <td>1:93552187:C:G</td>
      <td>rs11164870</td>
      <td>1</td>
      <td>93552187</td>
      <td>2.636249e-09</td>
      <td>93539383</td>
      <td>94007922</td>
      <td>189</td>
      <td>150</td>
      <td>1</td>
      <td>rs11164870</td>
      <td>1</td>
      <td>rs11164870</td>
    </tr>
    <tr>
      <th>1</th>
      <td>2</td>
      <td>1:155033317:C:T</td>
      <td>rs11264304</td>
      <td>1</td>
      <td>155033317</td>
      <td>4.634658e-12</td>
      <td>155024309</td>
      <td>155049402</td>
      <td>10</td>
      <td>9</td>
      <td>3</td>
      <td>rs11264304;rs11264305;rs3765087</td>
      <td>1</td>
      <td>rs11264304</td>
    </tr>
    <tr>
      <th>2</th>
      <td>3</td>
      <td>1:155388851:A:T</td>
      <td>rs12734374</td>
      <td>1</td>
      <td>155388851</td>
      <td>3.852096e-68</td>
      <td>155388851</td>
      <td>155388851</td>
      <td>1</td>
      <td>1</td>
      <td>1</td>
      <td>rs12734374</td>
      <td>1</td>
      <td>rs12734374</td>
    </tr>
    <tr>
      <th>3</th>
      <td>4</td>
      <td>1:156063880:A:C</td>
      <td>rs10737170</td>
      <td>1</td>
      <td>156063880</td>
      <td>2.569624e-10</td>
      <td>156063880</td>
      <td>156063880</td>
      <td>1</td>
      <td>1</td>
      <td>1</td>
      <td>rs10737170</td>
      <td>1</td>
      <td>rs10737170</td>
    </tr>
    <tr>
      <th>4</th>
      <td>5</td>
      <td>1:161469054:C:G</td>
      <td>rs6658353</td>
      <td>1</td>
      <td>161469054</td>
      <td>2.302802e-11</td>
      <td>161444369</td>
      <td>161479745</td>
      <td>9</td>
      <td>9</td>
      <td>2</td>
      <td>rs6658353;rs12120358</td>
      <td>1</td>
      <td>rs6658353</td>
    </tr>
  </tbody>
</table>
</div>



# 2. Fine-map using MR-MEGA Bayes Factor

The folowing function does the following:

1. We use FUMA's definition of locus start and end
2. convert log Bayes Factor (lnBF) to Bayes Factor (BF)
3. get sum of BF to calculate Posterior Probability (PP)
4. Append each locus to finemapped df
5. grab all loci with PP > 0.95 or 0.99

The output produces a list of two pandas dataframes:
1. Finemap report that show the total number of SNPs, number of SNPs in the 95% credible set, and number of SNPs in the 99% credible set.
2. All variants in the finemap results.


```python
def finemapLocus(leadSNP, MRMEGA):
    MRMEGA.loc[:,'BF'] = np.exp(MRMEGA.loc[:,'lnBF'])
    metaChrList = [g for _,g in MRMEGA.groupby('Chromosome')]
    credibleSet = list()
    Finemapped = pd.DataFrame()
    for i, row in leadSNP.iterrows():
        row = row.copy()
        chrom = row.chr
        start = row.start
        stop = row.end
        listNum = chrom-1
        chrMeta = metaChrList[listNum]
        locus_selection = (chrMeta.Position >= start) & (chrMeta.Position <= stop)
        locus = chrMeta.loc[locus_selection]
        locus = locus.sort_values(
            by=['BF'],
            ascending=False
        ).reset_index(drop=True)
        denom = locus['BF'].sum()
        #print(locus)
        # 99% credible set
        numer = locus['BF'].iloc[0]
        PP = numer/denom
        count = 1
        while PP < 0.99:
            numer = numer + locus['BF'][count]
            PP = numer/denom
            count = count + 1
        credibleSetDF = locus[0:count]
        credibleSetDF.loc[:,'PP'] = credibleSetDF.loc[:,'BF']/denom
        credibleSetDF.loc[:,'Locus'] = i+1
        # PP > 0.95
        numer = locus['BF'].iloc[0]
        PP = numer/denom
        count_95 = 1
        while PP < 0.95:
            numer = numer + locus['BF'][count_95]
            PP = numer/denom
            count_95 = count_95 + 1
        credibleSet95DF = locus[0:count_95]
        credibleSet95DF.loc[:,'PP'] = credibleSet95DF.loc[:,'BF']/denom
        credibleSet95DF.loc[:,'Locus'] = i+1
        
        FinemapReport = pd.DataFrame(
            {
                'Locus':i+1,
                'NumSNPs':locus.shape[0],
                'NumSNPs in 99% Credible Set':count,
                'NumSNPs in 95% Credible Set':count_95
            },
            index=[i]
        )
        
        # make final DFs
        # for list of SNPs, keep 95% set since it is larger and should include 99%
        if i==0:
            OutputFinemapReport = FinemapReport
            #Finemapped8 = Finemapped8_temp
            OutputCredibleSet = credibleSet95DF
        else:
            OutputFinemapReport = pd.concat([OutputFinemapReport, FinemapReport])
            #Finemapped8 = pd.concat([Finemapped8, Finemapped8_temp])
            OutputCredibleSet = pd.concat([OutputCredibleSet, credibleSet95DF])
    Final = [OutputFinemapReport, OutputCredibleSet]
    return(Final)
```


```python
finemap_res = finemapLocus(loci, sumstat)
```

    /opt/conda/lib/python3.7/site-packages/pandas/core/indexing.py:1667: SettingWithCopyWarning: 
    A value is trying to be set on a copy of a slice from a DataFrame.
    Try using .loc[row_indexer,col_indexer] = value instead
    
    See the caveats in the documentation: https://pandas.pydata.org/pandas-docs/stable/user_guide/indexing.html#returning-a-view-versus-a-copy
      self.obj[key] = value


Let's look at the report.


```python
finemap_res[0]
```




<div>
<style scoped>
    .dataframe tbody tr th:only-of-type {
        vertical-align: middle;
    }

    .dataframe tbody tr th {
        vertical-align: top;
    }

    .dataframe thead th {
        text-align: right;
    }
</style>
<table border="1" class="dataframe">
  <thead>
    <tr style="text-align: right;">
      <th></th>
      <th>Locus</th>
      <th>NumSNPs</th>
      <th>NumSNPs in 99% Credible Set</th>
      <th>NumSNPs in 95% Credible Set</th>
    </tr>
  </thead>
  <tbody>
    <tr>
      <th>0</th>
      <td>1</td>
      <td>508</td>
      <td>116</td>
      <td>97</td>
    </tr>
    <tr>
      <th>1</th>
      <td>2</td>
      <td>28</td>
      <td>5</td>
      <td>3</td>
    </tr>
    <tr>
      <th>2</th>
      <td>3</td>
      <td>1</td>
      <td>1</td>
      <td>1</td>
    </tr>
    <tr>
      <th>3</th>
      <td>4</td>
      <td>1</td>
      <td>1</td>
      <td>1</td>
    </tr>
    <tr>
      <th>4</th>
      <td>5</td>
      <td>66</td>
      <td>5</td>
      <td>3</td>
    </tr>
    <tr>
      <th>...</th>
      <td>...</td>
      <td>...</td>
      <td>...</td>
      <td>...</td>
    </tr>
    <tr>
      <th>58</th>
      <td>59</td>
      <td>603</td>
      <td>2</td>
      <td>2</td>
    </tr>
    <tr>
      <th>59</th>
      <td>60</td>
      <td>1</td>
      <td>1</td>
      <td>1</td>
    </tr>
    <tr>
      <th>60</th>
      <td>61</td>
      <td>363</td>
      <td>84</td>
      <td>55</td>
    </tr>
    <tr>
      <th>61</th>
      <td>62</td>
      <td>25</td>
      <td>9</td>
      <td>4</td>
    </tr>
    <tr>
      <th>62</th>
      <td>63</td>
      <td>292</td>
      <td>27</td>
      <td>23</td>
    </tr>
  </tbody>
</table>
<p>63 rows × 4 columns</p>
</div>



Let's look at loci with a single SNP in the 95% credible set.


```python
finemap_res[0][finemap_res[0]['NumSNPs in 95% Credible Set']==1]
```




<div>
<style scoped>
    .dataframe tbody tr th:only-of-type {
        vertical-align: middle;
    }

    .dataframe tbody tr th {
        vertical-align: top;
    }

    .dataframe thead th {
        text-align: right;
    }
</style>
<table border="1" class="dataframe">
  <thead>
    <tr style="text-align: right;">
      <th></th>
      <th>Locus</th>
      <th>NumSNPs</th>
      <th>NumSNPs in 99% Credible Set</th>
      <th>NumSNPs in 95% Credible Set</th>
    </tr>
  </thead>
  <tbody>
    <tr>
      <th>2</th>
      <td>3</td>
      <td>1</td>
      <td>1</td>
      <td>1</td>
    </tr>
    <tr>
      <th>3</th>
      <td>4</td>
      <td>1</td>
      <td>1</td>
      <td>1</td>
    </tr>
    <tr>
      <th>10</th>
      <td>11</td>
      <td>6</td>
      <td>1</td>
      <td>1</td>
    </tr>
    <tr>
      <th>18</th>
      <td>19</td>
      <td>926</td>
      <td>1</td>
      <td>1</td>
    </tr>
    <tr>
      <th>22</th>
      <td>23</td>
      <td>1483</td>
      <td>1</td>
      <td>1</td>
    </tr>
    <tr>
      <th>23</th>
      <td>24</td>
      <td>121</td>
      <td>2</td>
      <td>1</td>
    </tr>
    <tr>
      <th>44</th>
      <td>45</td>
      <td>1371</td>
      <td>1</td>
      <td>1</td>
    </tr>
    <tr>
      <th>55</th>
      <td>56</td>
      <td>1060</td>
      <td>1</td>
      <td>1</td>
    </tr>
    <tr>
      <th>59</th>
      <td>60</td>
      <td>1</td>
      <td>1</td>
      <td>1</td>
    </tr>
  </tbody>
</table>
</div>



Locus 2, 3, and 4 are all part of the _GBA_ locus from a previous meta-analysis. Let's look at locus 11 which is the _TMEM163_ locus.


```python
finemap_res[1].loc[finemap_res[1]['Locus']==11]
```




<div>
<style scoped>
    .dataframe tbody tr th:only-of-type {
        vertical-align: middle;
    }

    .dataframe tbody tr th {
        vertical-align: top;
    }

    .dataframe thead th {
        text-align: right;
    }
</style>
<table border="1" class="dataframe">
  <thead>
    <tr style="text-align: right;">
      <th></th>
      <th>MarkerName</th>
      <th>Chromosome</th>
      <th>Position</th>
      <th>EA</th>
      <th>NEA</th>
      <th>EAF</th>
      <th>Nsample</th>
      <th>Ncohort</th>
      <th>Effects</th>
      <th>beta_0</th>
      <th>...</th>
      <th>ndf_ancestry_het</th>
      <th>P-value_ancestry_het</th>
      <th>chisq_residual_het</th>
      <th>ndf_residual_het</th>
      <th>P-value_residual_het</th>
      <th>lnBF</th>
      <th>Comments</th>
      <th>BF</th>
      <th>PP</th>
      <th>Locus</th>
    </tr>
  </thead>
  <tbody>
    <tr>
      <th>0</th>
      <td>chr2:135464616</td>
      <td>2</td>
      <td>135464616</td>
      <td>A</td>
      <td>G</td>
      <td>0.640883</td>
      <td>2613720.0</td>
      <td>7.0</td>
      <td>+-+++-+</td>
      <td>0.046178</td>
      <td>...</td>
      <td>3.0</td>
      <td>0.002509</td>
      <td>3.15769</td>
      <td>3.0</td>
      <td>0.367946</td>
      <td>25.4868</td>
      <td>NaN</td>
      <td>1.171593e+11</td>
      <td>1.0</td>
      <td>11</td>
    </tr>
  </tbody>
</table>
<p>1 rows × 31 columns</p>
</div>



We will be keeping loci with 95% credible sets with 5 or less variants.


```python
LociWithLessThan5 = finemap_res[0][finemap_res[0]['NumSNPs in 95% Credible Set'] < 5]
LociWithLessThan5
```




<div>
<style scoped>
    .dataframe tbody tr th:only-of-type {
        vertical-align: middle;
    }

    .dataframe tbody tr th {
        vertical-align: top;
    }

    .dataframe thead th {
        text-align: right;
    }
</style>
<table border="1" class="dataframe">
  <thead>
    <tr style="text-align: right;">
      <th></th>
      <th>Locus</th>
      <th>NumSNPs</th>
      <th>NumSNPs in 99% Credible Set</th>
      <th>NumSNPs in 95% Credible Set</th>
    </tr>
  </thead>
  <tbody>
    <tr>
      <th>1</th>
      <td>2</td>
      <td>28</td>
      <td>5</td>
      <td>3</td>
    </tr>
    <tr>
      <th>2</th>
      <td>3</td>
      <td>1</td>
      <td>1</td>
      <td>1</td>
    </tr>
    <tr>
      <th>3</th>
      <td>4</td>
      <td>1</td>
      <td>1</td>
      <td>1</td>
    </tr>
    <tr>
      <th>4</th>
      <td>5</td>
      <td>66</td>
      <td>5</td>
      <td>3</td>
    </tr>
    <tr>
      <th>6</th>
      <td>7</td>
      <td>250</td>
      <td>3</td>
      <td>2</td>
    </tr>
    <tr>
      <th>10</th>
      <td>11</td>
      <td>6</td>
      <td>1</td>
      <td>1</td>
    </tr>
    <tr>
      <th>11</th>
      <td>12</td>
      <td>559</td>
      <td>3</td>
      <td>2</td>
    </tr>
    <tr>
      <th>12</th>
      <td>13</td>
      <td>52</td>
      <td>2</td>
      <td>2</td>
    </tr>
    <tr>
      <th>18</th>
      <td>19</td>
      <td>926</td>
      <td>1</td>
      <td>1</td>
    </tr>
    <tr>
      <th>19</th>
      <td>20</td>
      <td>732</td>
      <td>4</td>
      <td>4</td>
    </tr>
    <tr>
      <th>21</th>
      <td>22</td>
      <td>1013</td>
      <td>4</td>
      <td>4</td>
    </tr>
    <tr>
      <th>22</th>
      <td>23</td>
      <td>1483</td>
      <td>1</td>
      <td>1</td>
    </tr>
    <tr>
      <th>23</th>
      <td>24</td>
      <td>121</td>
      <td>2</td>
      <td>1</td>
    </tr>
    <tr>
      <th>24</th>
      <td>25</td>
      <td>131</td>
      <td>6</td>
      <td>4</td>
    </tr>
    <tr>
      <th>28</th>
      <td>29</td>
      <td>132</td>
      <td>4</td>
      <td>2</td>
    </tr>
    <tr>
      <th>30</th>
      <td>31</td>
      <td>275</td>
      <td>2</td>
      <td>2</td>
    </tr>
    <tr>
      <th>41</th>
      <td>42</td>
      <td>352</td>
      <td>3</td>
      <td>2</td>
    </tr>
    <tr>
      <th>42</th>
      <td>43</td>
      <td>1508</td>
      <td>4</td>
      <td>2</td>
    </tr>
    <tr>
      <th>44</th>
      <td>45</td>
      <td>1371</td>
      <td>1</td>
      <td>1</td>
    </tr>
    <tr>
      <th>45</th>
      <td>46</td>
      <td>288</td>
      <td>12</td>
      <td>4</td>
    </tr>
    <tr>
      <th>55</th>
      <td>56</td>
      <td>1060</td>
      <td>1</td>
      <td>1</td>
    </tr>
    <tr>
      <th>58</th>
      <td>59</td>
      <td>603</td>
      <td>2</td>
      <td>2</td>
    </tr>
    <tr>
      <th>59</th>
      <td>60</td>
      <td>1</td>
      <td>1</td>
      <td>1</td>
    </tr>
    <tr>
      <th>61</th>
      <td>62</td>
      <td>25</td>
      <td>9</td>
      <td>4</td>
    </tr>
  </tbody>
</table>
</div>




```python
LociWithLessThan5CredibleSet = finemap_res[1][
    finemap_res[1]['Locus'].isin(LociWithLessThan5.Locus)
]
LociWithLessThan5CredibleSet = LociWithLessThan5CredibleSet.rename(columns={
    "Chromosome":"CHR",
    "Position":"BP"
})
```

# 3. Add annotations to the credible sets

FUMA already provides us with annotation data. Let's add it to our results.


```python
FUMAsnps = pd.read_csv("data/FUMA/MRMEGA/snps.txt", delim_whitespace=True)
FUMAsnps.head()
```




<div>
<style scoped>
    .dataframe tbody tr th:only-of-type {
        vertical-align: middle;
    }

    .dataframe tbody tr th {
        vertical-align: top;
    }

    .dataframe thead th {
        text-align: right;
    }
</style>
<table border="1" class="dataframe">
  <thead>
    <tr style="text-align: right;">
      <th></th>
      <th>uniqID</th>
      <th>rsID</th>
      <th>chr</th>
      <th>pos</th>
      <th>non_effect_allele</th>
      <th>effect_allele</th>
      <th>MAF</th>
      <th>gwasP</th>
      <th>r2</th>
      <th>IndSigSNP</th>
      <th>...</th>
      <th>nearestGene</th>
      <th>dist</th>
      <th>func</th>
      <th>CADD</th>
      <th>RDB</th>
      <th>minChrState</th>
      <th>commonChrState</th>
      <th>posMapFilt</th>
      <th>eqtlMapFilt</th>
      <th>ciMapFilt</th>
    </tr>
  </thead>
  <tbody>
    <tr>
      <th>0</th>
      <td>1:93539383:A:G</td>
      <td>rs12385720</td>
      <td>1</td>
      <td>93539383</td>
      <td>G</td>
      <td>A</td>
      <td>0.4908</td>
      <td>1.597183e-08</td>
      <td>0.943336</td>
      <td>rs11164870</td>
      <td>...</td>
      <td>MTF2</td>
      <td>5408</td>
      <td>intergenic</td>
      <td>3.249</td>
      <td>5</td>
      <td>5</td>
      <td>15</td>
      <td>1</td>
      <td>1</td>
      <td>0</td>
    </tr>
    <tr>
      <th>1</th>
      <td>1:93539656:G:T</td>
      <td>rs7532549</td>
      <td>1</td>
      <td>93539656</td>
      <td>G</td>
      <td>T</td>
      <td>0.4844</td>
      <td>5.638348e-08</td>
      <td>0.883138</td>
      <td>rs11164870</td>
      <td>...</td>
      <td>MTF2</td>
      <td>5135</td>
      <td>intergenic</td>
      <td>4.759</td>
      <td>7</td>
      <td>5</td>
      <td>15</td>
      <td>1</td>
      <td>1</td>
      <td>0</td>
    </tr>
    <tr>
      <th>2</th>
      <td>1:93543483:C:G</td>
      <td>rs3754184</td>
      <td>1</td>
      <td>93543483</td>
      <td>G</td>
      <td>C</td>
      <td>0.4910</td>
      <td>7.975444e-09</td>
      <td>0.965467</td>
      <td>rs11164870</td>
      <td>...</td>
      <td>MTF2</td>
      <td>1308</td>
      <td>intergenic</td>
      <td>11.150</td>
      <td>6</td>
      <td>5</td>
      <td>15</td>
      <td>1</td>
      <td>1</td>
      <td>0</td>
    </tr>
    <tr>
      <th>3</th>
      <td>1:93544131:C:T</td>
      <td>rs1856027</td>
      <td>1</td>
      <td>93544131</td>
      <td>C</td>
      <td>T</td>
      <td>0.4930</td>
      <td>1.106524e-08</td>
      <td>0.959348</td>
      <td>rs11164870</td>
      <td>...</td>
      <td>MTF2</td>
      <td>660</td>
      <td>upstream</td>
      <td>6.132</td>
      <td>NaN</td>
      <td>1</td>
      <td>2</td>
      <td>1</td>
      <td>1</td>
      <td>0</td>
    </tr>
    <tr>
      <th>4</th>
      <td>1:93546299:G:T</td>
      <td>rs12128131</td>
      <td>1</td>
      <td>93546299</td>
      <td>G</td>
      <td>T</td>
      <td>0.4896</td>
      <td>4.563418e-09</td>
      <td>0.963671</td>
      <td>rs11164870</td>
      <td>...</td>
      <td>MTF2</td>
      <td>0</td>
      <td>intronic</td>
      <td>14.900</td>
      <td>2b</td>
      <td>1</td>
      <td>1</td>
      <td>1</td>
      <td>1</td>
      <td>0</td>
    </tr>
  </tbody>
</table>
<p>5 rows × 21 columns</p>
</div>




```python
# rename before merge
FUMAsnps = FUMAsnps.rename(columns={"chr":"CHR","pos":"BP"})
```


```python
FUMAsnps = FUMAsnps.drop(columns=[
    "uniqID",
    "non_effect_allele",
    "effect_allele",
    "MAF",
    "gwasP",
    "r2"
    ])
FUMAsnps.head()
```




<div>
<style scoped>
    .dataframe tbody tr th:only-of-type {
        vertical-align: middle;
    }

    .dataframe tbody tr th {
        vertical-align: top;
    }

    .dataframe thead th {
        text-align: right;
    }
</style>
<table border="1" class="dataframe">
  <thead>
    <tr style="text-align: right;">
      <th></th>
      <th>rsID</th>
      <th>CHR</th>
      <th>BP</th>
      <th>IndSigSNP</th>
      <th>GenomicLocus</th>
      <th>nearestGene</th>
      <th>dist</th>
      <th>func</th>
      <th>CADD</th>
      <th>RDB</th>
      <th>minChrState</th>
      <th>commonChrState</th>
      <th>posMapFilt</th>
      <th>eqtlMapFilt</th>
      <th>ciMapFilt</th>
    </tr>
  </thead>
  <tbody>
    <tr>
      <th>0</th>
      <td>rs12385720</td>
      <td>1</td>
      <td>93539383</td>
      <td>rs11164870</td>
      <td>1</td>
      <td>MTF2</td>
      <td>5408</td>
      <td>intergenic</td>
      <td>3.249</td>
      <td>5</td>
      <td>5</td>
      <td>15</td>
      <td>1</td>
      <td>1</td>
      <td>0</td>
    </tr>
    <tr>
      <th>1</th>
      <td>rs7532549</td>
      <td>1</td>
      <td>93539656</td>
      <td>rs11164870</td>
      <td>1</td>
      <td>MTF2</td>
      <td>5135</td>
      <td>intergenic</td>
      <td>4.759</td>
      <td>7</td>
      <td>5</td>
      <td>15</td>
      <td>1</td>
      <td>1</td>
      <td>0</td>
    </tr>
    <tr>
      <th>2</th>
      <td>rs3754184</td>
      <td>1</td>
      <td>93543483</td>
      <td>rs11164870</td>
      <td>1</td>
      <td>MTF2</td>
      <td>1308</td>
      <td>intergenic</td>
      <td>11.150</td>
      <td>6</td>
      <td>5</td>
      <td>15</td>
      <td>1</td>
      <td>1</td>
      <td>0</td>
    </tr>
    <tr>
      <th>3</th>
      <td>rs1856027</td>
      <td>1</td>
      <td>93544131</td>
      <td>rs11164870</td>
      <td>1</td>
      <td>MTF2</td>
      <td>660</td>
      <td>upstream</td>
      <td>6.132</td>
      <td>NaN</td>
      <td>1</td>
      <td>2</td>
      <td>1</td>
      <td>1</td>
      <td>0</td>
    </tr>
    <tr>
      <th>4</th>
      <td>rs12128131</td>
      <td>1</td>
      <td>93546299</td>
      <td>rs11164870</td>
      <td>1</td>
      <td>MTF2</td>
      <td>0</td>
      <td>intronic</td>
      <td>14.900</td>
      <td>2b</td>
      <td>1</td>
      <td>1</td>
      <td>1</td>
      <td>1</td>
      <td>0</td>
    </tr>
  </tbody>
</table>
</div>



## Merge our results with FUMA annotations


```python
LociWithLessThan5CredibleSet = LociWithLessThan5CredibleSet.rename(columns={
    "Chromosome":"CHR",
    "Position":"BP"
})
```


```python
Finemap95ResultAnnotated = pd.merge(FUMAsnps, LociWithLessThan5CredibleSet, on=["CHR","BP"])
LociWithLessThan5CredibleSet.shape[0],Finemap95ResultAnnotated.shape[0]
```




    (51, 51)




```python
Finemap95ResultAnnotated
```




<div>
<style scoped>
    .dataframe tbody tr th:only-of-type {
        vertical-align: middle;
    }

    .dataframe tbody tr th {
        vertical-align: top;
    }

    .dataframe thead th {
        text-align: right;
    }
</style>
<table border="1" class="dataframe">
  <thead>
    <tr style="text-align: right;">
      <th></th>
      <th>rsID</th>
      <th>CHR</th>
      <th>BP</th>
      <th>IndSigSNP</th>
      <th>GenomicLocus</th>
      <th>nearestGene</th>
      <th>dist</th>
      <th>func</th>
      <th>CADD</th>
      <th>RDB</th>
      <th>...</th>
      <th>ndf_ancestry_het</th>
      <th>P-value_ancestry_het</th>
      <th>chisq_residual_het</th>
      <th>ndf_residual_het</th>
      <th>P-value_residual_het</th>
      <th>lnBF</th>
      <th>Comments</th>
      <th>BF</th>
      <th>PP</th>
      <th>Locus</th>
    </tr>
  </thead>
  <tbody>
    <tr>
      <th>0</th>
      <td>rs11264304</td>
      <td>1</td>
      <td>155033317</td>
      <td>rs11264304</td>
      <td>2</td>
      <td>ADAM15</td>
      <td>0</td>
      <td>intronic</td>
      <td>6.405</td>
      <td>5</td>
      <td>...</td>
      <td>3.0</td>
      <td>0.985364</td>
      <td>4.321360</td>
      <td>3.0</td>
      <td>0.228789</td>
      <td>25.6239</td>
      <td>NaN</td>
      <td>1.343750e+11</td>
      <td>0.836328</td>
      <td>2</td>
    </tr>
    <tr>
      <th>1</th>
      <td>rs12743272</td>
      <td>1</td>
      <td>155033334</td>
      <td>rs3765087</td>
      <td>2</td>
      <td>ADAM15</td>
      <td>0</td>
      <td>intronic</td>
      <td>6.694</td>
      <td>5</td>
      <td>...</td>
      <td>3.0</td>
      <td>0.238156</td>
      <td>1.548020</td>
      <td>2.0</td>
      <td>0.461160</td>
      <td>22.0701</td>
      <td>NaN</td>
      <td>3.845233e+09</td>
      <td>0.023932</td>
      <td>2</td>
    </tr>
    <tr>
      <th>2</th>
      <td>rs11264305</td>
      <td>1</td>
      <td>155033572</td>
      <td>rs11264305</td>
      <td>2</td>
      <td>ADAM15</td>
      <td>0</td>
      <td>intronic</td>
      <td>12.260</td>
      <td>5</td>
      <td>...</td>
      <td>3.0</td>
      <td>0.180145</td>
      <td>2.055370</td>
      <td>2.0</td>
      <td>0.357834</td>
      <td>23.4966</td>
      <td>NaN</td>
      <td>1.601193e+10</td>
      <td>0.099656</td>
      <td>2</td>
    </tr>
    <tr>
      <th>3</th>
      <td>rs12734374</td>
      <td>1</td>
      <td>155388851</td>
      <td>rs12734374</td>
      <td>3</td>
      <td>ASH1L</td>
      <td>0</td>
      <td>intronic</td>
      <td>0.037</td>
      <td>6</td>
      <td>...</td>
      <td>3.0</td>
      <td>0.000090</td>
      <td>1.454400</td>
      <td>2.0</td>
      <td>0.483260</td>
      <td>156.7270</td>
      <td>NaN</td>
      <td>1.163245e+68</td>
      <td>1.000000</td>
      <td>3</td>
    </tr>
    <tr>
      <th>4</th>
      <td>rs10737170</td>
      <td>1</td>
      <td>156063880</td>
      <td>rs10737170</td>
      <td>4</td>
      <td>LMNA</td>
      <td>0</td>
      <td>intronic</td>
      <td>0.528</td>
      <td>1f</td>
      <td>...</td>
      <td>3.0</td>
      <td>0.507340</td>
      <td>0.883665</td>
      <td>3.0</td>
      <td>0.829368</td>
      <td>21.4619</td>
      <td>NaN</td>
      <td>2.093075e+09</td>
      <td>1.000000</td>
      <td>4</td>
    </tr>
    <tr>
      <th>5</th>
      <td>rs6658353</td>
      <td>1</td>
      <td>161469054</td>
      <td>rs6658353</td>
      <td>5</td>
      <td>FCGR2A</td>
      <td>6165</td>
      <td>intergenic</td>
      <td>0.756</td>
      <td>5</td>
      <td>...</td>
      <td>3.0</td>
      <td>0.293572</td>
      <td>1.754040</td>
      <td>3.0</td>
      <td>0.624987</td>
      <td>23.9648</td>
      <td>NaN</td>
      <td>2.557292e+10</td>
      <td>0.614565</td>
      <td>5</td>
    </tr>
    <tr>
      <th>6</th>
      <td>rs4657041</td>
      <td>1</td>
      <td>161478859</td>
      <td>rs6658353</td>
      <td>5</td>
      <td>FCGR2A</td>
      <td>0</td>
      <td>intronic</td>
      <td>3.890</td>
      <td>6</td>
      <td>...</td>
      <td>3.0</td>
      <td>0.224335</td>
      <td>1.276000</td>
      <td>3.0</td>
      <td>0.734840</td>
      <td>22.8998</td>
      <td>NaN</td>
      <td>8.815699e+09</td>
      <td>0.211858</td>
      <td>5</td>
    </tr>
    <tr>
      <th>7</th>
      <td>rs1801274</td>
      <td>1</td>
      <td>161479745</td>
      <td>rs6658353</td>
      <td>5</td>
      <td>FCGR2A</td>
      <td>0</td>
      <td>exonic</td>
      <td>0.979</td>
      <td>NaN</td>
      <td>...</td>
      <td>3.0</td>
      <td>0.961453</td>
      <td>0.649296</td>
      <td>2.0</td>
      <td>0.722782</td>
      <td>22.5217</td>
      <td>NaN</td>
      <td>6.040182e+09</td>
      <td>0.145157</td>
      <td>5</td>
    </tr>
    <tr>
      <th>8</th>
      <td>rs12132270</td>
      <td>1</td>
      <td>205660940</td>
      <td>rs12118655</td>
      <td>7</td>
      <td>SLC45A3</td>
      <td>11352</td>
      <td>intergenic</td>
      <td>6.298</td>
      <td>4</td>
      <td>...</td>
      <td>3.0</td>
      <td>0.113914</td>
      <td>7.954620</td>
      <td>3.0</td>
      <td>0.046959</td>
      <td>72.0997</td>
      <td>NaN</td>
      <td>2.053534e+31</td>
      <td>0.198025</td>
      <td>7</td>
    </tr>
    <tr>
      <th>9</th>
      <td>rs12118655</td>
      <td>1</td>
      <td>205663478</td>
      <td>rs12118655</td>
      <td>7</td>
      <td>SLC45A3</td>
      <td>13890</td>
      <td>intergenic</td>
      <td>1.395</td>
      <td>5</td>
      <td>...</td>
      <td>3.0</td>
      <td>0.112337</td>
      <td>8.624070</td>
      <td>3.0</td>
      <td>0.034730</td>
      <td>73.4705</td>
      <td>NaN</td>
      <td>8.087843e+31</td>
      <td>0.779922</td>
      <td>7</td>
    </tr>
    <tr>
      <th>10</th>
      <td>rs57891859</td>
      <td>2</td>
      <td>135464616</td>
      <td>rs57891859</td>
      <td>11</td>
      <td>TMEM163</td>
      <td>0</td>
      <td>intronic</td>
      <td>6.746</td>
      <td>4</td>
      <td>...</td>
      <td>3.0</td>
      <td>0.002509</td>
      <td>3.157690</td>
      <td>3.0</td>
      <td>0.367946</td>
      <td>25.4868</td>
      <td>NaN</td>
      <td>1.171593e+11</td>
      <td>1.000000</td>
      <td>11</td>
    </tr>
    <tr>
      <th>11</th>
      <td>rs4613239</td>
      <td>2</td>
      <td>169119609</td>
      <td>rs12987123</td>
      <td>12</td>
      <td>STK39</td>
      <td>14957</td>
      <td>intergenic</td>
      <td>1.304</td>
      <td>4</td>
      <td>...</td>
      <td>3.0</td>
      <td>0.000720</td>
      <td>0.634202</td>
      <td>3.0</td>
      <td>0.888560</td>
      <td>92.8813</td>
      <td>NaN</td>
      <td>2.176888e+40</td>
      <td>0.034742</td>
      <td>12</td>
    </tr>
    <tr>
      <th>12</th>
      <td>rs12987123</td>
      <td>2</td>
      <td>169120541</td>
      <td>rs12987123</td>
      <td>12</td>
      <td>STK39</td>
      <td>15889</td>
      <td>intergenic</td>
      <td>0.754</td>
      <td>7</td>
      <td>...</td>
      <td>3.0</td>
      <td>0.001040</td>
      <td>0.777887</td>
      <td>3.0</td>
      <td>0.854749</td>
      <td>96.1788</td>
      <td>NaN</td>
      <td>5.887381e+41</td>
      <td>0.939583</td>
      <td>12</td>
    </tr>
    <tr>
      <th>13</th>
      <td>rs1461806</td>
      <td>3</td>
      <td>28700178</td>
      <td>rs6808178</td>
      <td>13</td>
      <td>LINC00693</td>
      <td>0</td>
      <td>ncRNA_intronic</td>
      <td>14.470</td>
      <td>6</td>
      <td>...</td>
      <td>3.0</td>
      <td>0.291694</td>
      <td>5.348370</td>
      <td>3.0</td>
      <td>0.147994</td>
      <td>25.0212</td>
      <td>NaN</td>
      <td>7.354770e+10</td>
      <td>0.428409</td>
      <td>13</td>
    </tr>
    <tr>
      <th>14</th>
      <td>rs6808178</td>
      <td>3</td>
      <td>28705690</td>
      <td>rs6808178</td>
      <td>13</td>
      <td>LINC00693</td>
      <td>0</td>
      <td>ncRNA_intronic</td>
      <td>2.285</td>
      <td>7</td>
      <td>...</td>
      <td>3.0</td>
      <td>0.439954</td>
      <td>5.472520</td>
      <td>3.0</td>
      <td>0.140291</td>
      <td>25.2980</td>
      <td>NaN</td>
      <td>9.700225e+10</td>
      <td>0.565030</td>
      <td>13</td>
    </tr>
    <tr>
      <th>15</th>
      <td>rs34311866</td>
      <td>4</td>
      <td>951947</td>
      <td>rs34311866</td>
      <td>19</td>
      <td>TMEM175</td>
      <td>0</td>
      <td>exonic</td>
      <td>11.090</td>
      <td>NaN</td>
      <td>...</td>
      <td>3.0</td>
      <td>0.221375</td>
      <td>2.513020</td>
      <td>3.0</td>
      <td>0.472943</td>
      <td>172.2960</td>
      <td>NaN</td>
      <td>6.717413e+74</td>
      <td>1.000000</td>
      <td>19</td>
    </tr>
    <tr>
      <th>16</th>
      <td>rs34559912</td>
      <td>4</td>
      <td>15730146</td>
      <td>rs4698412</td>
      <td>20</td>
      <td>BST1</td>
      <td>0</td>
      <td>intronic</td>
      <td>2.609</td>
      <td>7</td>
      <td>...</td>
      <td>3.0</td>
      <td>0.783888</td>
      <td>3.383080</td>
      <td>3.0</td>
      <td>0.336246</td>
      <td>73.2081</td>
      <td>NaN</td>
      <td>6.221196e+31</td>
      <td>0.077991</td>
      <td>20</td>
    </tr>
    <tr>
      <th>17</th>
      <td>rs11724635</td>
      <td>4</td>
      <td>15737101</td>
      <td>rs4698412</td>
      <td>20</td>
      <td>BST1</td>
      <td>0</td>
      <td>intronic</td>
      <td>2.859</td>
      <td>7</td>
      <td>...</td>
      <td>3.0</td>
      <td>0.865342</td>
      <td>3.172530</td>
      <td>3.0</td>
      <td>0.365782</td>
      <td>74.5933</td>
      <td>NaN</td>
      <td>2.485756e+32</td>
      <td>0.311623</td>
      <td>20</td>
    </tr>
    <tr>
      <th>18</th>
      <td>rs4698412</td>
      <td>4</td>
      <td>15737348</td>
      <td>rs4698412</td>
      <td>20</td>
      <td>BST1</td>
      <td>0</td>
      <td>intronic</td>
      <td>0.378</td>
      <td>6</td>
      <td>...</td>
      <td>3.0</td>
      <td>0.860964</td>
      <td>3.469910</td>
      <td>3.0</td>
      <td>0.324686</td>
      <td>74.6579</td>
      <td>NaN</td>
      <td>2.651637e+32</td>
      <td>0.332418</td>
      <td>20</td>
    </tr>
    <tr>
      <th>19</th>
      <td>rs4698413</td>
      <td>4</td>
      <td>15737882</td>
      <td>rs4698412</td>
      <td>20</td>
      <td>BST1</td>
      <td>0</td>
      <td>intronic</td>
      <td>0.310</td>
      <td>6</td>
      <td>...</td>
      <td>3.0</td>
      <td>0.853307</td>
      <td>3.510120</td>
      <td>3.0</td>
      <td>0.319452</td>
      <td>74.4775</td>
      <td>NaN</td>
      <td>2.213947e+32</td>
      <td>0.277548</td>
      <td>20</td>
    </tr>
    <tr>
      <th>20</th>
      <td>rs1946959</td>
      <td>4</td>
      <td>77196677</td>
      <td>rs6854006</td>
      <td>22</td>
      <td>FAM47E:FAM47E-STBD1:FAM47E-STBD1</td>
      <td>0:0:0</td>
      <td>intronic</td>
      <td>1.602</td>
      <td>4</td>
      <td>...</td>
      <td>3.0</td>
      <td>0.000815</td>
      <td>8.655390</td>
      <td>3.0</td>
      <td>0.034242</td>
      <td>42.3692</td>
      <td>NaN</td>
      <td>2.515995e+18</td>
      <td>0.140867</td>
      <td>22</td>
    </tr>
    <tr>
      <th>21</th>
      <td>rs6854006</td>
      <td>4</td>
      <td>77198054</td>
      <td>rs6854006</td>
      <td>22</td>
      <td>FAM47E:FAM47E-STBD1:FAM47E-STBD1</td>
      <td>0:0:0</td>
      <td>intronic</td>
      <td>3.577</td>
      <td>7</td>
      <td>...</td>
      <td>3.0</td>
      <td>0.000336</td>
      <td>8.817200</td>
      <td>3.0</td>
      <td>0.031823</td>
      <td>43.2435</td>
      <td>NaN</td>
      <td>6.031335e+18</td>
      <td>0.337686</td>
      <td>22</td>
    </tr>
    <tr>
      <th>22</th>
      <td>rs6858114</td>
      <td>4</td>
      <td>77198465</td>
      <td>rs6854006</td>
      <td>22</td>
      <td>FAM47E:FAM47E-STBD1:FAM47E-STBD1</td>
      <td>0:0:0</td>
      <td>intronic</td>
      <td>0.313</td>
      <td>7</td>
      <td>...</td>
      <td>3.0</td>
      <td>0.000328</td>
      <td>8.676510</td>
      <td>3.0</td>
      <td>0.033916</td>
      <td>42.6225</td>
      <td>NaN</td>
      <td>3.241280e+18</td>
      <td>0.181475</td>
      <td>22</td>
    </tr>
    <tr>
      <th>23</th>
      <td>rs6812193</td>
      <td>4</td>
      <td>77198986</td>
      <td>rs6854006</td>
      <td>22</td>
      <td>FAM47E:FAM47E-STBD1:FAM47E-STBD1</td>
      <td>0:0:0</td>
      <td>intronic</td>
      <td>3.887</td>
      <td>5</td>
      <td>...</td>
      <td>3.0</td>
      <td>0.000362</td>
      <td>8.636250</td>
      <td>3.0</td>
      <td>0.034539</td>
      <td>43.2227</td>
      <td>NaN</td>
      <td>5.907179e+18</td>
      <td>0.330735</td>
      <td>22</td>
    </tr>
    <tr>
      <th>24</th>
      <td>rs356182</td>
      <td>4</td>
      <td>90626111</td>
      <td>rs356182</td>
      <td>23</td>
      <td>RP11-115D19.1</td>
      <td>0</td>
      <td>ncRNA_intronic</td>
      <td>8.962</td>
      <td>NaN</td>
      <td>...</td>
      <td>3.0</td>
      <td>0.011963</td>
      <td>26.809400</td>
      <td>2.0</td>
      <td>0.000002</td>
      <td>392.4750</td>
      <td>NaN</td>
      <td>2.816610e+170</td>
      <td>1.000000</td>
      <td>23</td>
    </tr>
    <tr>
      <th>25</th>
      <td>rs13117519</td>
      <td>4</td>
      <td>114369065</td>
      <td>rs13117519</td>
      <td>24</td>
      <td>CAMK2D</td>
      <td>3122</td>
      <td>intergenic</td>
      <td>1.216</td>
      <td>3a</td>
      <td>...</td>
      <td>3.0</td>
      <td>0.450316</td>
      <td>0.579619</td>
      <td>3.0</td>
      <td>0.901082</td>
      <td>23.5522</td>
      <td>NaN</td>
      <td>1.692741e+10</td>
      <td>0.957090</td>
      <td>24</td>
    </tr>
    <tr>
      <th>26</th>
      <td>rs7434295</td>
      <td>4</td>
      <td>170563946</td>
      <td>rs7434295</td>
      <td>25</td>
      <td>CLCN3</td>
      <td>0</td>
      <td>intronic</td>
      <td>0.397</td>
      <td>7</td>
      <td>...</td>
      <td>3.0</td>
      <td>0.167195</td>
      <td>0.303847</td>
      <td>3.0</td>
      <td>0.959303</td>
      <td>24.3998</td>
      <td>NaN</td>
      <td>3.950922e+10</td>
      <td>0.522003</td>
      <td>25</td>
    </tr>
    <tr>
      <th>27</th>
      <td>rs13136720</td>
      <td>4</td>
      <td>170565154</td>
      <td>rs7434295</td>
      <td>25</td>
      <td>CLCN3</td>
      <td>0</td>
      <td>intronic</td>
      <td>1.071</td>
      <td>7</td>
      <td>...</td>
      <td>3.0</td>
      <td>0.150207</td>
      <td>0.315964</td>
      <td>3.0</td>
      <td>0.956999</td>
      <td>23.7489</td>
      <td>NaN</td>
      <td>2.060707e+10</td>
      <td>0.272264</td>
      <td>25</td>
    </tr>
    <tr>
      <th>28</th>
      <td>rs10866341</td>
      <td>4</td>
      <td>170568789</td>
      <td>rs7434295</td>
      <td>25</td>
      <td>CLCN3</td>
      <td>0</td>
      <td>intronic</td>
      <td>0.390</td>
      <td>7</td>
      <td>...</td>
      <td>3.0</td>
      <td>0.195838</td>
      <td>0.318706</td>
      <td>3.0</td>
      <td>0.956473</td>
      <td>22.0247</td>
      <td>NaN</td>
      <td>3.674563e+09</td>
      <td>0.048549</td>
      <td>25</td>
    </tr>
    <tr>
      <th>29</th>
      <td>rs6839362</td>
      <td>4</td>
      <td>170569481</td>
      <td>rs7434295</td>
      <td>25</td>
      <td>CLCN3</td>
      <td>0</td>
      <td>intronic</td>
      <td>4.124</td>
      <td>7</td>
      <td>...</td>
      <td>3.0</td>
      <td>0.200495</td>
      <td>0.297938</td>
      <td>3.0</td>
      <td>0.960416</td>
      <td>22.9138</td>
      <td>NaN</td>
      <td>8.939987e+09</td>
      <td>0.118117</td>
      <td>25</td>
    </tr>
    <tr>
      <th>30</th>
      <td>rs12528068</td>
      <td>6</td>
      <td>72487762</td>
      <td>rs12528068</td>
      <td>29</td>
      <td>RIMS1</td>
      <td>108643</td>
      <td>intergenic</td>
      <td>1.447</td>
      <td>7</td>
      <td>...</td>
      <td>3.0</td>
      <td>0.118501</td>
      <td>2.746140</td>
      <td>3.0</td>
      <td>0.432443</td>
      <td>25.1815</td>
      <td>NaN</td>
      <td>8.633492e+10</td>
      <td>0.520105</td>
      <td>29</td>
    </tr>
    <tr>
      <th>31</th>
      <td>rs66472504</td>
      <td>6</td>
      <td>72489955</td>
      <td>rs12528068</td>
      <td>29</td>
      <td>RIMS1</td>
      <td>106450</td>
      <td>intergenic</td>
      <td>3.381</td>
      <td>6</td>
      <td>...</td>
      <td>3.0</td>
      <td>0.119626</td>
      <td>2.733950</td>
      <td>3.0</td>
      <td>0.434489</td>
      <td>25.0262</td>
      <td>NaN</td>
      <td>7.391636e+10</td>
      <td>0.445292</td>
      <td>29</td>
    </tr>
    <tr>
      <th>32</th>
      <td>rs41286192</td>
      <td>6</td>
      <td>133118216</td>
      <td>rs41286192</td>
      <td>31</td>
      <td>SLC18B1</td>
      <td>0</td>
      <td>exonic</td>
      <td>17.920</td>
      <td>5</td>
      <td>...</td>
      <td>3.0</td>
      <td>0.726686</td>
      <td>6.371500</td>
      <td>3.0</td>
      <td>0.094870</td>
      <td>23.6716</td>
      <td>NaN</td>
      <td>1.907415e+10</td>
      <td>0.779946</td>
      <td>31</td>
    </tr>
    <tr>
      <th>33</th>
      <td>rs75859381</td>
      <td>6</td>
      <td>133210361</td>
      <td>rs41286192</td>
      <td>31</td>
      <td>HMGB1P13</td>
      <td>20362</td>
      <td>intergenic</td>
      <td>8.679</td>
      <td>7</td>
      <td>...</td>
      <td>3.0</td>
      <td>0.832996</td>
      <td>5.760440</td>
      <td>3.0</td>
      <td>0.123865</td>
      <td>22.3998</td>
      <td>NaN</td>
      <td>5.346992e+09</td>
      <td>0.218640</td>
      <td>31</td>
    </tr>
    <tr>
      <th>34</th>
      <td>rs3802921</td>
      <td>11</td>
      <td>133786993</td>
      <td>rs3802920</td>
      <td>42</td>
      <td>IGSF9B</td>
      <td>0</td>
      <td>UTR3</td>
      <td>4.303</td>
      <td>NaN</td>
      <td>...</td>
      <td>3.0</td>
      <td>0.080109</td>
      <td>7.571620</td>
      <td>3.0</td>
      <td>0.055747</td>
      <td>41.0956</td>
      <td>NaN</td>
      <td>7.040319e+17</td>
      <td>0.055693</td>
      <td>42</td>
    </tr>
    <tr>
      <th>35</th>
      <td>rs3802920</td>
      <td>11</td>
      <td>133787001</td>
      <td>rs3802920</td>
      <td>42</td>
      <td>IGSF9B</td>
      <td>0</td>
      <td>UTR3</td>
      <td>4.382</td>
      <td>NaN</td>
      <td>...</td>
      <td>3.0</td>
      <td>0.047665</td>
      <td>6.403810</td>
      <td>3.0</td>
      <td>0.093534</td>
      <td>43.8738</td>
      <td>NaN</td>
      <td>1.132790e+19</td>
      <td>0.896110</td>
      <td>42</td>
    </tr>
    <tr>
      <th>36</th>
      <td>rs28370649</td>
      <td>12</td>
      <td>40399149</td>
      <td>rs28370650</td>
      <td>43</td>
      <td>SLC2A13</td>
      <td>0</td>
      <td>intronic</td>
      <td>1.257</td>
      <td>7</td>
      <td>...</td>
      <td>3.0</td>
      <td>0.610961</td>
      <td>0.421330</td>
      <td>2.0</td>
      <td>0.810045</td>
      <td>75.6803</td>
      <td>NaN</td>
      <td>7.371174e+32</td>
      <td>0.433335</td>
      <td>43</td>
    </tr>
    <tr>
      <th>37</th>
      <td>rs28370650</td>
      <td>12</td>
      <td>40399948</td>
      <td>rs28370650</td>
      <td>43</td>
      <td>SLC2A13</td>
      <td>0</td>
      <td>intronic</td>
      <td>0.575</td>
      <td>6</td>
      <td>...</td>
      <td>3.0</td>
      <td>0.609166</td>
      <td>0.314997</td>
      <td>2.0</td>
      <td>0.854278</td>
      <td>75.8620</td>
      <td>NaN</td>
      <td>8.839913e+32</td>
      <td>0.519679</td>
      <td>43</td>
    </tr>
    <tr>
      <th>38</th>
      <td>rs10847864</td>
      <td>12</td>
      <td>123326598</td>
      <td>rs10847864</td>
      <td>45</td>
      <td>HIP1R</td>
      <td>0</td>
      <td>intronic</td>
      <td>2.403</td>
      <td>2b</td>
      <td>...</td>
      <td>3.0</td>
      <td>0.041826</td>
      <td>2.688250</td>
      <td>2.0</td>
      <td>0.260768</td>
      <td>94.7720</td>
      <td>NaN</td>
      <td>1.441971e+41</td>
      <td>0.999994</td>
      <td>45</td>
    </tr>
    <tr>
      <th>39</th>
      <td>rs11614702</td>
      <td>12</td>
      <td>133058157</td>
      <td>rs11610045</td>
      <td>46</td>
      <td>MUC8</td>
      <td>7430</td>
      <td>intergenic</td>
      <td>0.384</td>
      <td>4</td>
      <td>...</td>
      <td>3.0</td>
      <td>0.873849</td>
      <td>0.061361</td>
      <td>2.0</td>
      <td>0.969785</td>
      <td>19.6856</td>
      <td>NaN</td>
      <td>3.542807e+08</td>
      <td>0.019039</td>
      <td>46</td>
    </tr>
    <tr>
      <th>40</th>
      <td>rs11610045</td>
      <td>12</td>
      <td>133063768</td>
      <td>rs11610045</td>
      <td>46</td>
      <td>FBRSL1</td>
      <td>2368</td>
      <td>intergenic</td>
      <td>1.226</td>
      <td>5</td>
      <td>...</td>
      <td>3.0</td>
      <td>0.384693</td>
      <td>0.483400</td>
      <td>2.0</td>
      <td>0.785292</td>
      <td>23.5522</td>
      <td>NaN</td>
      <td>1.692741e+10</td>
      <td>0.909698</td>
      <td>46</td>
    </tr>
    <tr>
      <th>41</th>
      <td>rs35318451</td>
      <td>12</td>
      <td>133068484</td>
      <td>rs36098511</td>
      <td>46</td>
      <td>FBRSL1</td>
      <td>0</td>
      <td>intronic</td>
      <td>7.476</td>
      <td>2b</td>
      <td>...</td>
      <td>3.0</td>
      <td>0.182425</td>
      <td>1.352020</td>
      <td>3.0</td>
      <td>0.716819</td>
      <td>18.9510</td>
      <td>NaN</td>
      <td>1.699475e+08</td>
      <td>0.009133</td>
      <td>46</td>
    </tr>
    <tr>
      <th>42</th>
      <td>rs36098511</td>
      <td>12</td>
      <td>133080449</td>
      <td>rs36098511</td>
      <td>46</td>
      <td>FBRSL1</td>
      <td>0</td>
      <td>intronic</td>
      <td>3.374</td>
      <td>5</td>
      <td>...</td>
      <td>3.0</td>
      <td>0.095038</td>
      <td>1.119290</td>
      <td>3.0</td>
      <td>0.772420</td>
      <td>19.3137</td>
      <td>NaN</td>
      <td>2.442493e+08</td>
      <td>0.013126</td>
      <td>46</td>
    </tr>
    <tr>
      <th>43</th>
      <td>rs62053943</td>
      <td>17</td>
      <td>43744203</td>
      <td>rs62053943</td>
      <td>56</td>
      <td>CRHR1:RP11-105N13.4</td>
      <td>0:0</td>
      <td>ncRNA_intronic</td>
      <td>0.183</td>
      <td>5</td>
      <td>...</td>
      <td>3.0</td>
      <td>0.682861</td>
      <td>0.326851</td>
      <td>2.0</td>
      <td>0.849230</td>
      <td>160.1550</td>
      <td>NaN</td>
      <td>3.584534e+69</td>
      <td>1.000000</td>
      <td>56</td>
    </tr>
    <tr>
      <th>44</th>
      <td>rs4588066</td>
      <td>18</td>
      <td>40672964</td>
      <td>rs4588066</td>
      <td>59</td>
      <td>RIT2</td>
      <td>0</td>
      <td>intronic</td>
      <td>1.520</td>
      <td>7</td>
      <td>...</td>
      <td>3.0</td>
      <td>0.587721</td>
      <td>1.472170</td>
      <td>3.0</td>
      <td>0.688708</td>
      <td>64.0531</td>
      <td>NaN</td>
      <td>6.575184e+27</td>
      <td>0.594257</td>
      <td>59</td>
    </tr>
    <tr>
      <th>45</th>
      <td>rs4130047</td>
      <td>18</td>
      <td>40678235</td>
      <td>rs4588066</td>
      <td>59</td>
      <td>RIT2</td>
      <td>0</td>
      <td>intronic</td>
      <td>0.414</td>
      <td>6</td>
      <td>...</td>
      <td>3.0</td>
      <td>0.563879</td>
      <td>1.803620</td>
      <td>3.0</td>
      <td>0.614148</td>
      <td>63.6715</td>
      <td>NaN</td>
      <td>4.489326e+27</td>
      <td>0.405740</td>
      <td>59</td>
    </tr>
    <tr>
      <th>46</th>
      <td>rs55818311</td>
      <td>19</td>
      <td>2341047</td>
      <td>rs55818311</td>
      <td>60</td>
      <td>SPPL2B</td>
      <td>0</td>
      <td>ncRNA_exonic</td>
      <td>1.096</td>
      <td>5</td>
      <td>...</td>
      <td>3.0</td>
      <td>0.662466</td>
      <td>3.862420</td>
      <td>2.0</td>
      <td>0.144973</td>
      <td>20.0879</td>
      <td>NaN</td>
      <td>5.297417e+08</td>
      <td>1.000000</td>
      <td>60</td>
    </tr>
    <tr>
      <th>47</th>
      <td>rs1297255</td>
      <td>21</td>
      <td>16804367</td>
      <td>rs1736020</td>
      <td>62</td>
      <td>AJ006998.2</td>
      <td>29742</td>
      <td>intergenic</td>
      <td>10.050</td>
      <td>7</td>
      <td>...</td>
      <td>3.0</td>
      <td>0.000775</td>
      <td>1.297430</td>
      <td>3.0</td>
      <td>0.729743</td>
      <td>16.0256</td>
      <td>NaN</td>
      <td>9.116532e+06</td>
      <td>0.011867</td>
      <td>62</td>
    </tr>
    <tr>
      <th>48</th>
      <td>rs1736135</td>
      <td>21</td>
      <td>16805220</td>
      <td>rs1736020</td>
      <td>62</td>
      <td>AJ006998.2</td>
      <td>30595</td>
      <td>intergenic</td>
      <td>6.535</td>
      <td>7</td>
      <td>...</td>
      <td>3.0</td>
      <td>0.000198</td>
      <td>1.544240</td>
      <td>3.0</td>
      <td>0.672098</td>
      <td>17.4270</td>
      <td>NaN</td>
      <td>3.702115e+07</td>
      <td>0.048190</td>
      <td>62</td>
    </tr>
    <tr>
      <th>49</th>
      <td>rs1297257</td>
      <td>21</td>
      <td>16806607</td>
      <td>rs1736020</td>
      <td>62</td>
      <td>AJ006998.2</td>
      <td>31982</td>
      <td>intergenic</td>
      <td>0.799</td>
      <td>5</td>
      <td>...</td>
      <td>3.0</td>
      <td>0.000111</td>
      <td>1.443100</td>
      <td>3.0</td>
      <td>0.695464</td>
      <td>18.1021</td>
      <td>NaN</td>
      <td>7.271804e+07</td>
      <td>0.094656</td>
      <td>62</td>
    </tr>
    <tr>
      <th>50</th>
      <td>rs1736020</td>
      <td>21</td>
      <td>16812552</td>
      <td>rs1736020</td>
      <td>62</td>
      <td>AJ006998.2</td>
      <td>37927</td>
      <td>intergenic</td>
      <td>4.723</td>
      <td>5</td>
      <td>...</td>
      <td>3.0</td>
      <td>0.000047</td>
      <td>0.897962</td>
      <td>2.0</td>
      <td>0.638278</td>
      <td>20.2374</td>
      <td>NaN</td>
      <td>6.151643e+08</td>
      <td>0.800748</td>
      <td>62</td>
    </tr>
  </tbody>
</table>
<p>51 rows × 44 columns</p>
</div>




```python
Finemap95ResultAnnotated['func'].unique()
```




    array(['intronic', 'intergenic', 'exonic', 'ncRNA_intronic', 'UTR3',
           'ncRNA_exonic'], dtype=object)




```python
Finemap95ResultAnnotated[Finemap95ResultAnnotated['func']=='exonic']
```




<div>
<style scoped>
    .dataframe tbody tr th:only-of-type {
        vertical-align: middle;
    }

    .dataframe tbody tr th {
        vertical-align: top;
    }

    .dataframe thead th {
        text-align: right;
    }
</style>
<table border="1" class="dataframe">
  <thead>
    <tr style="text-align: right;">
      <th></th>
      <th>rsID</th>
      <th>CHR</th>
      <th>BP</th>
      <th>IndSigSNP</th>
      <th>GenomicLocus</th>
      <th>nearestGene</th>
      <th>dist</th>
      <th>func</th>
      <th>CADD</th>
      <th>RDB</th>
      <th>...</th>
      <th>ndf_ancestry_het</th>
      <th>P-value_ancestry_het</th>
      <th>chisq_residual_het</th>
      <th>ndf_residual_het</th>
      <th>P-value_residual_het</th>
      <th>lnBF</th>
      <th>Comments</th>
      <th>BF</th>
      <th>PP</th>
      <th>Locus</th>
    </tr>
  </thead>
  <tbody>
    <tr>
      <th>7</th>
      <td>rs1801274</td>
      <td>1</td>
      <td>161479745</td>
      <td>rs6658353</td>
      <td>5</td>
      <td>FCGR2A</td>
      <td>0</td>
      <td>exonic</td>
      <td>0.979</td>
      <td>NaN</td>
      <td>...</td>
      <td>3.0</td>
      <td>0.961453</td>
      <td>0.649296</td>
      <td>2.0</td>
      <td>0.722782</td>
      <td>22.5217</td>
      <td>NaN</td>
      <td>6.040182e+09</td>
      <td>0.145157</td>
      <td>5</td>
    </tr>
    <tr>
      <th>15</th>
      <td>rs34311866</td>
      <td>4</td>
      <td>951947</td>
      <td>rs34311866</td>
      <td>19</td>
      <td>TMEM175</td>
      <td>0</td>
      <td>exonic</td>
      <td>11.090</td>
      <td>NaN</td>
      <td>...</td>
      <td>3.0</td>
      <td>0.221375</td>
      <td>2.513020</td>
      <td>3.0</td>
      <td>0.472943</td>
      <td>172.2960</td>
      <td>NaN</td>
      <td>6.717413e+74</td>
      <td>1.000000</td>
      <td>19</td>
    </tr>
    <tr>
      <th>32</th>
      <td>rs41286192</td>
      <td>6</td>
      <td>133118216</td>
      <td>rs41286192</td>
      <td>31</td>
      <td>SLC18B1</td>
      <td>0</td>
      <td>exonic</td>
      <td>17.920</td>
      <td>5</td>
      <td>...</td>
      <td>3.0</td>
      <td>0.726686</td>
      <td>6.371500</td>
      <td>3.0</td>
      <td>0.094870</td>
      <td>23.6716</td>
      <td>NaN</td>
      <td>1.907415e+10</td>
      <td>0.779946</td>
      <td>31</td>
    </tr>
  </tbody>
</table>
<p>3 rows × 44 columns</p>
</div>




```python
Finemap95ResultAnnotated[Finemap95ResultAnnotated['GenomicLocus'].isin([5,8,31,47])]
```




<div>
<style scoped>
    .dataframe tbody tr th:only-of-type {
        vertical-align: middle;
    }

    .dataframe tbody tr th {
        vertical-align: top;
    }

    .dataframe thead th {
        text-align: right;
    }
</style>
<table border="1" class="dataframe">
  <thead>
    <tr style="text-align: right;">
      <th></th>
      <th>rsID</th>
      <th>CHR</th>
      <th>BP</th>
      <th>IndSigSNP</th>
      <th>GenomicLocus</th>
      <th>nearestGene</th>
      <th>dist</th>
      <th>func</th>
      <th>CADD</th>
      <th>RDB</th>
      <th>...</th>
      <th>ndf_ancestry_het</th>
      <th>P-value_ancestry_het</th>
      <th>chisq_residual_het</th>
      <th>ndf_residual_het</th>
      <th>P-value_residual_het</th>
      <th>lnBF</th>
      <th>Comments</th>
      <th>BF</th>
      <th>PP</th>
      <th>Locus</th>
    </tr>
  </thead>
  <tbody>
    <tr>
      <th>5</th>
      <td>rs6658353</td>
      <td>1</td>
      <td>161469054</td>
      <td>rs6658353</td>
      <td>5</td>
      <td>FCGR2A</td>
      <td>6165</td>
      <td>intergenic</td>
      <td>0.756</td>
      <td>5</td>
      <td>...</td>
      <td>3.0</td>
      <td>0.293572</td>
      <td>1.754040</td>
      <td>3.0</td>
      <td>0.624987</td>
      <td>23.9648</td>
      <td>NaN</td>
      <td>2.557292e+10</td>
      <td>0.614565</td>
      <td>5</td>
    </tr>
    <tr>
      <th>6</th>
      <td>rs4657041</td>
      <td>1</td>
      <td>161478859</td>
      <td>rs6658353</td>
      <td>5</td>
      <td>FCGR2A</td>
      <td>0</td>
      <td>intronic</td>
      <td>3.890</td>
      <td>6</td>
      <td>...</td>
      <td>3.0</td>
      <td>0.224335</td>
      <td>1.276000</td>
      <td>3.0</td>
      <td>0.734840</td>
      <td>22.8998</td>
      <td>NaN</td>
      <td>8.815699e+09</td>
      <td>0.211858</td>
      <td>5</td>
    </tr>
    <tr>
      <th>7</th>
      <td>rs1801274</td>
      <td>1</td>
      <td>161479745</td>
      <td>rs6658353</td>
      <td>5</td>
      <td>FCGR2A</td>
      <td>0</td>
      <td>exonic</td>
      <td>0.979</td>
      <td>NaN</td>
      <td>...</td>
      <td>3.0</td>
      <td>0.961453</td>
      <td>0.649296</td>
      <td>2.0</td>
      <td>0.722782</td>
      <td>22.5217</td>
      <td>NaN</td>
      <td>6.040182e+09</td>
      <td>0.145157</td>
      <td>5</td>
    </tr>
    <tr>
      <th>32</th>
      <td>rs41286192</td>
      <td>6</td>
      <td>133118216</td>
      <td>rs41286192</td>
      <td>31</td>
      <td>SLC18B1</td>
      <td>0</td>
      <td>exonic</td>
      <td>17.920</td>
      <td>5</td>
      <td>...</td>
      <td>3.0</td>
      <td>0.726686</td>
      <td>6.371500</td>
      <td>3.0</td>
      <td>0.094870</td>
      <td>23.6716</td>
      <td>NaN</td>
      <td>1.907415e+10</td>
      <td>0.779946</td>
      <td>31</td>
    </tr>
    <tr>
      <th>33</th>
      <td>rs75859381</td>
      <td>6</td>
      <td>133210361</td>
      <td>rs41286192</td>
      <td>31</td>
      <td>HMGB1P13</td>
      <td>20362</td>
      <td>intergenic</td>
      <td>8.679</td>
      <td>7</td>
      <td>...</td>
      <td>3.0</td>
      <td>0.832996</td>
      <td>5.760440</td>
      <td>3.0</td>
      <td>0.123865</td>
      <td>22.3998</td>
      <td>NaN</td>
      <td>5.346992e+09</td>
      <td>0.218640</td>
      <td>31</td>
    </tr>
  </tbody>
</table>
<p>5 rows × 44 columns</p>
</div>




```python
Finemap95ResultAnnotated['RDB'].unique()
```




    array(['5', '6', '1f', nan, '4', '7', '3a', '2b'], dtype=object)



## Organize before output


```python
# drop unnecessary columns
Finemap95ResultAnnotated = Finemap95ResultAnnotated.drop(columns=[
    "IndSigSNP",
    "GenomicLocus",
    "dist",
    "chisq_ancestry_het",
    "chisq_residual_het",
    "ndf_ancestry_het",
    "ndf_residual_het",
    "lnBF",
    "BF",
    "Comments",
    "MarkerName",
    "beta_0",
    "beta_1",
    "beta_2",
    "beta_3",
    "se_0",
    "se_1",
    "se_2",
    "se_3",
    "Ncohort",
    "Effects",
    "EAF",
    "chisq_association",
    "ndf_association",
    "Nsample",
    "posMapFilt",
    "eqtlMapFilt",
    "ciMapFilt"
])
```


```python
# sort by Locus + PP
Finemap95ResultAnnotated = Finemap95ResultAnnotated.sort_values(
    by=["Locus", "PP"],
    ascending=[True,False]
).reset_index(drop=True)
Finemap95ResultAnnotated.head()
```




<div>
<style scoped>
    .dataframe tbody tr th:only-of-type {
        vertical-align: middle;
    }

    .dataframe tbody tr th {
        vertical-align: top;
    }

    .dataframe thead th {
        text-align: right;
    }
</style>
<table border="1" class="dataframe">
  <thead>
    <tr style="text-align: right;">
      <th></th>
      <th>rsID</th>
      <th>CHR</th>
      <th>BP</th>
      <th>nearestGene</th>
      <th>func</th>
      <th>CADD</th>
      <th>RDB</th>
      <th>minChrState</th>
      <th>commonChrState</th>
      <th>EA</th>
      <th>NEA</th>
      <th>P-value_association</th>
      <th>P-value_ancestry_het</th>
      <th>P-value_residual_het</th>
      <th>PP</th>
      <th>Locus</th>
    </tr>
  </thead>
  <tbody>
    <tr>
      <th>0</th>
      <td>rs11264304</td>
      <td>1</td>
      <td>155033317</td>
      <td>ADAM15</td>
      <td>intronic</td>
      <td>6.405</td>
      <td>5</td>
      <td>3</td>
      <td>5</td>
      <td>T</td>
      <td>C</td>
      <td>4.634658e-12</td>
      <td>0.985364</td>
      <td>0.228789</td>
      <td>0.836328</td>
      <td>2</td>
    </tr>
    <tr>
      <th>1</th>
      <td>rs11264305</td>
      <td>1</td>
      <td>155033572</td>
      <td>ADAM15</td>
      <td>intronic</td>
      <td>12.260</td>
      <td>5</td>
      <td>3</td>
      <td>5</td>
      <td>A</td>
      <td>G</td>
      <td>4.871242e-11</td>
      <td>0.180145</td>
      <td>0.357834</td>
      <td>0.099656</td>
      <td>2</td>
    </tr>
    <tr>
      <th>2</th>
      <td>rs12743272</td>
      <td>1</td>
      <td>155033334</td>
      <td>ADAM15</td>
      <td>intronic</td>
      <td>6.694</td>
      <td>5</td>
      <td>3</td>
      <td>5</td>
      <td>A</td>
      <td>G</td>
      <td>1.925480e-10</td>
      <td>0.238156</td>
      <td>0.461160</td>
      <td>0.023932</td>
      <td>2</td>
    </tr>
    <tr>
      <th>3</th>
      <td>rs12734374</td>
      <td>1</td>
      <td>155388851</td>
      <td>ASH1L</td>
      <td>intronic</td>
      <td>0.037</td>
      <td>6</td>
      <td>4</td>
      <td>5</td>
      <td>A</td>
      <td>T</td>
      <td>3.852096e-68</td>
      <td>0.000090</td>
      <td>0.483260</td>
      <td>1.000000</td>
      <td>3</td>
    </tr>
    <tr>
      <th>4</th>
      <td>rs10737170</td>
      <td>1</td>
      <td>156063880</td>
      <td>LMNA</td>
      <td>intronic</td>
      <td>0.528</td>
      <td>1f</td>
      <td>2</td>
      <td>5</td>
      <td>A</td>
      <td>C</td>
      <td>2.569624e-10</td>
      <td>0.507340</td>
      <td>0.829368</td>
      <td>1.000000</td>
      <td>4</td>
    </tr>
  </tbody>
</table>
</div>




```python
Finemap95ResultAnnotated.columns
```




    Index(['rsID', 'CHR', 'BP', 'nearestGene', 'func', 'CADD', 'RDB',
           'minChrState', 'commonChrState', 'EA', 'NEA', 'P-value_association',
           'P-value_ancestry_het', 'P-value_residual_het', 'PP', 'Locus'],
          dtype='object')




```python
# rename columns
Finemap95ResultAnnotated = Finemap95ResultAnnotated.rename(columns={
    "nearestGene":"Nearest Gene/Feature",
    "func":"Functional Consequence",
    "minChrState":"Minimum Chromatin State",
    "commonChrState":"Most Common Chromatin State",
    "EA":"A1",
    "NEA":"A2",
    "P-value_association":"P(MR-MEGA)",
    "P-value_ancestry_het":"P(ANC-HET)",
    "P-value_residual_het":"P(RES-HET)",
})
```


```python
# reorder columns
Finemap95ResultAnnotated = Finemap95ResultAnnotated[[
    'Locus','rsID','CHR','BP','A1','A2',
    'P(MR-MEGA)',"P(ANC-HET)","P(RES-HET)",
    'Nearest Gene/Feature',
    'Functional Consequence',"CADD","RDB",
    "Minimum Chromatin State","Most Common Chromatin State",
    "PP"
]]
Finemap95ResultAnnotated.head()
```




<div>
<style scoped>
    .dataframe tbody tr th:only-of-type {
        vertical-align: middle;
    }

    .dataframe tbody tr th {
        vertical-align: top;
    }

    .dataframe thead th {
        text-align: right;
    }
</style>
<table border="1" class="dataframe">
  <thead>
    <tr style="text-align: right;">
      <th></th>
      <th>Locus</th>
      <th>rsID</th>
      <th>CHR</th>
      <th>BP</th>
      <th>A1</th>
      <th>A2</th>
      <th>P(MR-MEGA)</th>
      <th>P(ANC-HET)</th>
      <th>P(RES-HET)</th>
      <th>Nearest Gene/Feature</th>
      <th>Functional Consequence</th>
      <th>CADD</th>
      <th>RDB</th>
      <th>Minimum Chromatin State</th>
      <th>Most Common Chromatin State</th>
      <th>PP</th>
    </tr>
  </thead>
  <tbody>
    <tr>
      <th>0</th>
      <td>2</td>
      <td>rs11264304</td>
      <td>1</td>
      <td>155033317</td>
      <td>T</td>
      <td>C</td>
      <td>4.634658e-12</td>
      <td>0.985364</td>
      <td>0.228789</td>
      <td>ADAM15</td>
      <td>intronic</td>
      <td>6.405</td>
      <td>5</td>
      <td>3</td>
      <td>5</td>
      <td>0.836328</td>
    </tr>
    <tr>
      <th>1</th>
      <td>2</td>
      <td>rs11264305</td>
      <td>1</td>
      <td>155033572</td>
      <td>A</td>
      <td>G</td>
      <td>4.871242e-11</td>
      <td>0.180145</td>
      <td>0.357834</td>
      <td>ADAM15</td>
      <td>intronic</td>
      <td>12.260</td>
      <td>5</td>
      <td>3</td>
      <td>5</td>
      <td>0.099656</td>
    </tr>
    <tr>
      <th>2</th>
      <td>2</td>
      <td>rs12743272</td>
      <td>1</td>
      <td>155033334</td>
      <td>A</td>
      <td>G</td>
      <td>1.925480e-10</td>
      <td>0.238156</td>
      <td>0.461160</td>
      <td>ADAM15</td>
      <td>intronic</td>
      <td>6.694</td>
      <td>5</td>
      <td>3</td>
      <td>5</td>
      <td>0.023932</td>
    </tr>
    <tr>
      <th>3</th>
      <td>3</td>
      <td>rs12734374</td>
      <td>1</td>
      <td>155388851</td>
      <td>A</td>
      <td>T</td>
      <td>3.852096e-68</td>
      <td>0.000090</td>
      <td>0.483260</td>
      <td>ASH1L</td>
      <td>intronic</td>
      <td>0.037</td>
      <td>6</td>
      <td>4</td>
      <td>5</td>
      <td>1.000000</td>
    </tr>
    <tr>
      <th>4</th>
      <td>4</td>
      <td>rs10737170</td>
      <td>1</td>
      <td>156063880</td>
      <td>A</td>
      <td>C</td>
      <td>2.569624e-10</td>
      <td>0.507340</td>
      <td>0.829368</td>
      <td>LMNA</td>
      <td>intronic</td>
      <td>0.528</td>
      <td>1f</td>
      <td>2</td>
      <td>5</td>
      <td>1.000000</td>
    </tr>
  </tbody>
</table>
</div>




```python
finemap_res[0].to_csv("output/FineMap.report.txt", sep = "\t", index=False)
finemap_res[1].to_csv("output/FineMap.95CredibleSets.txt", sep = "\t", index=False)
LociWithLessThan5.to_csv("output/FineMap.report.Lessthan5.txt", sep = "\t", index=False)
Finemap95ResultAnnotated.to_csv("output/FineMap.95CredibleSets.Lessthan5.Annotated.txt", sep = "\t", index=False)
```


```python

```
