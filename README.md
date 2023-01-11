# AutoLC-BO

[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)
[![Open In Colab](https://colab.research.google.com/assets/colab-badge.svg)](https://colab.research.google.com/drive/1ebAcH8z0IRuR3kD7Kus-ZLIVqaUKVljo?usp=sharing)

Code related to publication titled: [Closed-loop automatic gradient design for liquid chromatography using Bayesian optimization](https://chemrxiv.org/engage/chemrxiv/article-details/62e2a383e7fc8f9e388caabc)

<img src="/figures/example.png" width="1000"/>

Please take a look at the Jupyter Notebook or checkout the [Google Collab](https://colab.research.google.com/drive/1ebAcH8z0IRuR3kD7Kus-ZLIVqaUKVljo?usp=sharing) which can be run without any prior installation. This covers similar code as used in Section 4.5 of the paper.

In the folder example interface code, we show C++ and Python code that can be used to do closed-loop Bayesian Optimization with Agilent software and hardware. This code is less straightforward to run as it requires some changes to the specific setup at hand (Paths, Chemstation instance, Slack support, etc.). Please do not hesitate to email us if you need any clarifications.

When using this code, please cite: [![DOI](https://data.caltech.edu/badge/110025475.svg)](https://doi.org/10.1016/j.aca.2023.340789) 
