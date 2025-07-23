## Space-time data-driven modeling of precipitation-induced shallow landslides in South Tyrol, Italy

This is the R code to reproduce the analysis in "Space-time data-driven modeling of wildfire initiation in the mountainous region of Trentino--South Tyrol, Italy."
by Mateo Moreno <sup>[![](https://info.orcid.org/wp-content/uploads/2020/12/orcid_16x16.gif)](https://orcid.org/0000-0002-9530-3076)</sup>,
Luigi Lombardo <sup>[![](https://info.orcid.org/wp-content/uploads/2020/12/orcid_16x16.gif)](https://orcid.org/0000-0003-4348-7288)</sup>.
Alice Crespi <sup>[![](https://info.orcid.org/wp-content/uploads/2020/12/orcid_16x16.gif)](https://orcid.org/0000-0003-4186-8474)</sup>,
Peter Zellner <sup>[![](https://info.orcid.org/wp-content/uploads/2020/12/orcid_16x16.gif)](https://orcid.org/0000-0002-3394-9664)</sup>,
Volkmair Mair <sup>[![](https://info.orcid.org/wp-content/uploads/2020/12/orcid_16x16.gif)]</sup>,
Massimiliano Pittore <sup>[![](https://info.orcid.org/wp-content/uploads/2020/12/orcid_16x16.gif)](https://orcid.org/0000-0003-4940-3444)</sup>,
Cees van Westen <sup>[![](https://info.orcid.org/wp-content/uploads/2020/12/orcid_16x16.gif)](https://orcid.org/0000-0002-2992-902X)</sup>
and Stefan Steger <sup>[![](https://info.orcid.org/wp-content/uploads/2020/12/orcid_16x16.gif)](https://orcid.org/0000-0003-0886-5191)</sup>
[![DOI](https://zenodo.org/badge/887347458.svg)](https://doi.org/10.5281/zenodo.15033256)

> Moreno, M., Lombardo, L., Crespi, A., Zellner, P. J., Mair, V., Pittore, M., van Westen, C. J., & Steger, S. (2024). Space-time data-driven modeling of precipitation-induced shallow landslides in South Tyrol, Italy. Science of The Total Environment, 912(169166). <https://doi.org/10.1016/j.scitotenv.2023.169166>

## Abstract

Shallow landslides represent potentially damaging processes in mountain areas worldwide. These geomorphic processes are usually caused by an interplay of predisposing, preparatory, and triggering environmental factors. At regional scales, data-driven methods have been used to model shallow landslides by addressing the spatial and temporal components separately. So far, few studies have explored the integration of space and time for landslide prediction. This research leverages generalized additive mixed models to develop an integrated approach to model shallow landslides in space and time. We built upon data on precipitation-induced landslide records from 2000 to 2020 in South Tyrol, Italy (7400 km2). The slope unit-based model predicts landslide occurrence as a function of static and dynamic factors while seasonal effects are incorporated. The model also accounts for spatial and temporal biases inherent in the underlying landslide data. We validated the resulting predictions through a suite of cross-validation techniques, obtaining consistent performance scores above 0.85. The analyses revealed that the best-performing model combines static ground conditions and two precipitation time windows: a short-term cumulative precipitation of 2 days before the landslide event and a medium-term cumulative precipitation of 14 days. We demonstrated the model's predictive capabilities by predicting the dynamic landslide probabilities over historical data associated with a heavy precipitation event on August 4th and August 5th, 2016, and hypothetical non-spatially explicit precipitation (what-if) scenarios. The novel approach shows the potential to integrate static and dynamic landslide factors for large areas, accounting for the underlying data structure and data limitations.

## Acknowledgements

We thank the Faculty of Geo-information Science and Earth Observation (ITC), University of Twente, for covering the open-access publication costs. The authors also thank the Provincial Office for Geology and Building Materials Testing of the Autonomous Province of Bozen/Bolzano for assisting in preparing landslide data. The research leading to these results is framed within the PROSLIDE project integration of static and dynamic landslide controls at multiple scales using data-driven and physically-based methods – exploring new opportunities for the PRediction Of shallow landSLIDEs (<https://www.mountainresearch.at/proslide/>) that received funding from the research program Research Südtirol/Alto Adige 2019 of the Autonomous Province of Bozen/Bolzano – Südtirol/Alto Adige (grant no. 9/34). Lastly, we thank the nine reviewers who gave constructive feedback and helped improve the quality of this work.

## Further dissemination

Previous work was presented at the EGU General Assembly 2022. The abstract is available at:

> Moreno, M., Steger, S., Lombardo, L., Crespi, A., Zellner, P. J., Pittore, M., Mair, V., & van Westen, C. (2022). Space-time modeling of rainfall-induced shallow landslides in South Tyrol, Italy. EGU General Assembly 2022, EGU22-9175. <https://doi.org/10.5194/egusphere-egu22-9175>

## Repo structure

The general structure is as follows:
- `dat`: data sets
- `dev`: development (scripts)
- `plt`: plots
