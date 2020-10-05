## Code and Data for: Crowding and the shape of COVID-19 epidemics

### Citation
Crowding and the shape of COVID-19 epidemics

Rader B, Scarpino SV, Nande A, Hill A, Adlam B, Reiner RC, Pigott DM, Gutierrez B, Zarebski A, Shrestha M,
open COVID-19 data working group, Brownstein JS, Castro MC, Tian H, Pybus OG, Kraemer MUG. In press.
Crowding and the shape of COVID-19 epidemics. _in press_ Nature Medicine. 

doi: https://doi.org/10.1038/s41591-020-1104-0

manuscript version of the code and data can be accessed via Zenodo:  [![DOI](https://zenodo.org/badge/299311245.svg)](https://zenodo.org/badge/latestdoi/299311245)

### Acknowledgements
We want to thank all the individuals and organizations across the world who have been willing and able to report data in as open and timely manner as possible. To see individuals involved in the often painstaking data curation process, please see [beoutbreakprepared/nCoV2019](https://github.com/beoutbreakprepared/nCoV2019) and our correspondence in The Lancet Infectious Diseases, ["Open access epidemiological data from the COVID-19 outbreak"](https://www.thelancet.com/journals/laninf/article/PIIS1473-3099(20)30119-5/fulltext).

*For a complete list of the COVID-19 data working group see [Specific Contributors](https://github.com/beoutbreakprepared/nCoV2019).

### Abstract
The Coronavirus Disease 2019 (COVID-19) pandemic is straining public health systems worldwide, and major non-pharmaceutical interventions have been implemented to slow its spread[1, 2, 3, 4]. During the initial phase of the outbreak, dissemination of severe acute respiratory syndrome coronavirus 2 (SARS-CoV-2) was primarily determined by human mobility from Wuhan, China[5, 6]. Yet empirical evidence on the effect of key geographic factors on local epidemic transmission is lacking[7]. In this study, we analyzed highly resolved spatial variables in cities, together with case count data, to investigate the role of climate, urbanization and variation in interventions. We show that the degree to which cases of COVID-19 are compressed into a short period of time (peakedness of the epidemic) is strongly shaped by population aggregation and heterogeneity, such that epidemics in crowded cities are more spread over time, and crowded cities have larger total attack rates than less populated cities. Observed differences in the peakedness of epidemics are consistent with a meta-population model of COVID-19 that explicitly accounts for spatial hierarchies. We paired our estimates with globally comprehensive data on human mobility and predict that crowded cities worldwide could experience more prolonged epidemics.

### Notes on the code
1. In order to run the analysis, you will need to obtain permission from Google to access the Google COVID-19 Aggregated Mobility Research Dataset.

### Data
1. The Google COVID-19 Aggregated Mobility Research Dataset used for this study is available with permission from Google. Once you have obtained written permission from Google, we can share the following data sets, which are necessary to run the code: "Google COVID-19 Aggregated Mobility Research Dataset.csv", italy_admin2_mobility_reductions_5_6_20.RDS, 

2. Mobility data were made public by [Baidu](https://www.baidu.com/) and manually transcribed from [Qianxi Baidu](https://qianxi.baidu.com/) on February 10th, 2020 and are stored [here](https://docs.google.com/spreadsheets/d/1ov7Z2IjEPRB41rmRe4oTF1nI6L3EB9_Bn64AJIkOX6g/edit#gid=0). We are publishing these data for research purposes under [Article 22 of the Copyright Law of the People's Republic of China](https://wipolex.wipo.int/en/text/466268) and US 9th Circuit Court of Appeals [No. 17-16783 D.C. No. 3:17-cv-03301-EMC](http://cdn.ca9.uscourts.gov/datastore/opinions/2019/09/09/17-16783.pdf).  Baidu's Intellectual Property statement may be found [here](https://www.baidu.com/duty/copyright.html). A version translated using Google translate is below.  Please see the Warranty section below as we make no representation as to the suitability of these data for any purpose.  We thank Baidu for making these data public and encourage them to continue to do so.  For more information please see: [Baidu service](https://qianxi.baidu.com/), this [gnews post](https://gnews.org/91700/), and this [blog post](https://zhuanlan.zhihu.com/p/104119625). 

Intellectual Property Statement (translated by Google translate via the [University of Virginia Biocomplexity Institute](https://dataverse.lib.virginia.edu/dataset.xhtml?persistentId=doi:10.18130/V3/YQLJ5W)) 

Baidu owns the copyright of all the materials in this website. If there are special provisions in the statement of rights of each sub-channel, those provisions shall prevail. Any authorized viewing, copying, printing and dissemination of materials belonging to this website must meet the following conditions:

* All materials and images are for informational purposes ;
* All materials and images must not be used for commercial purposes ;
* All materials, images and any part thereof must include this copyright notice ;

All products, technologies, and all programs on this website (www.baidu.com) belong to Baidu's intellectual property and are not authorized here. "Baidu", "Baidu" and related graphics are registered trademarks of Baidu.

Without Baidu's permission, no one may use (including but not limited to: illegally copying, disseminating, displaying, mirroring, uploading, downloading) or using unconventional means (such as maliciously interfering with Baidu data) to affect Baidu ’s normal Service, no one may obtain Baidu data automatically by software program without authorization. Otherwise, Baidu will pursue legal responsibility according to law.

3.   To see individuals involved in the often painstaking data curation process, please see [beoutbreakprepared/nCoV2019](https://github.com/beoutbreakprepared/nCoV2019) and our correspondence in The Lancet Infectious Diseases, ["Open access epidemiological data from the COVID-19 outbreak"](https://www.thelancet.com/journals/laninf/article/PIIS1473-3099(20)30119-5/fulltext).

4. The population data should be cited using, "Doxsey-Whitfield, E. et al. Taking advantage of the improved availability of census data: a first look at the gridded population of the world, version 4. Pap. Appl. Geogr. 1, 226–234 (2015)." Please see https://darksky.net/forecast/40.7127,-74.0059/us12/en for more information.

### License
(see LICENSE)

## Additional license, warranty, and copyright information
We provide a license for our code (see LICENSE) and do not claim ownership, nor the right to license, the data we have obtained nor any third-party software tools/code used in our analyses.  Please cite the appropriate agency, paper, and/or individual in publications and/or derivatives using these data, contact them regarding the legal use of these data, and remember to pass-forward any existing license/warranty/copyright information.  THE DATA AND SOFTWARE ARE PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NON-INFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE DATA AND/OR SOFTWARE OR THE USE OR OTHER DEALINGS IN THE DATA AND/OR SOFTWARE.