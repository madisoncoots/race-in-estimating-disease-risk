# Race and Ethnicity in Estimating Disease Risk

Data and replication materials for Coots et al. (2024) "A Framework for Considering the Value of Race and Ethnicity in Estimating Disease Risk".

To reproduce our figures, run the `code/make_paper_figures.R` script. R data objects for the analyses are located in the `data/` repository, so it is not necessary to run any data preprocessing scripts to generate the figures. 

To reproduce the processed data files used in our analysis, run the `code/make_data_all_diseases.R` script followed by the `code/make_risk_predictions_all_diseases.R` script to make risk predictions on the processed data. 

The data used for our lung cancer analysis are not publicly available, so we are not able to provide these data for replication. Consequently, figures for lung cancer cannot be regenerated directly. However, all code used to produce our results for lung cancer is located in `code/lc`. Data from the National Lung Screening Trial is availably by request from the National Cancer Institute.

All figures that appear in our paper may be found in `figures`.

___

For further detail, please see the structure of this repository outlined below:

```bash
code/                     # Contains all code used for our analyses.
├── make_data_all_diseases.R              # Master script that runs data preprocessing scripts for each disease.
├── make_risk_predictions_all_diseases.R  # Master script that runs prediction scripts for each disease.
├── make_paper_figures.R                  # Master script that runs analysis scripts for each disease and constructs aggregate figures.
│
├── cvd/                   # Contains all code used for our cardiovascular disease analysis.
│   ├── make_data_cvd.R                 # Data preprocessing script for CVD. Pulls NHANES tables directly from the CDC website.
│   ├── 2013_pce.R                      # Script that generates CVD risk predictions using the 2013 pooled cohort equations.
│   ├── 2013_pce_constants.R            # Script containing weights for the 2013 pooled cohort equations.
│   ├── figures_analysis.R              # Script that produces all CVD subfigures and accompanying statistics.
│   └── colors.R                        # Script that specifies group colors for figures.
│
├── bc/                    # Contains all code used for our breast cancer analysis.
│   ├── make_data_bc.R                  # Data preprocessing script for breast cancer. Pulls NHANES tables directly from the CDC website.
│   ├── gail_equation.R                 # Script that generates breast cancer risk predictions using the Gail equations.
│   ├── figures_analysis.R              # Script that produces all breast cancer subfigures and accompanying statistics.
│   └── colors.R                        # Script that specifies group colors for figures.
│
├── lc/                    # Contains all code used for our lung cancer analysis.
│   ├── make_data_lc.R                  # Data preprocessing script for lung cancer using NLST data for generating risk predictions.
│   ├── lc_risk_predictions.R           # Script that generates lung cancer risk predictions using the LCRAT model.
│   ├── make_census_weights.R           # Script that generates survey weights for the NLST data using census data.
│   ├── figures_analysis.R              # Script that produces all lung cancer subfigures and accompanying statistics.
│   └── colors.R                        # Script that specifies group colors for figures.
│
data/                   # Contains raw and processed data used in our analyses.
figures/                # Contains PDF files of all figures in our paper.

```



