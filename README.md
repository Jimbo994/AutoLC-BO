# AutoLC-BO

[![CC BY-NC-SA 4.0][cc-by-nc-sa-shield]][cc-by-nc-sa]
[![Open In Colab](https://colab.research.google.com/assets/colab-badge.svg)](https://colab.research.google.com/drive/1ebAcH8z0IRuR3kD7Kus-ZLIVqaUKVljo?usp=sharing)

Code related to publication titled: [Closed-loop automatic gradient design for liquid chromatography using Bayesian optimization](https://doi.org/10.1016/j.aca.2023.340789)

<img src="/figures/example.png" width="1000"/>

Please take a look at the Jupyter Notebook or checkout the [Google Collab](https://colab.research.google.com/drive/1ebAcH8z0IRuR3kD7Kus-ZLIVqaUKVljo?usp=sharing) which can be run without any prior installation. This covers similar code as used in Section 4.5 of the paper.

In the folder example interface code, we show C++ and Python code that can be used to do closed-loop Bayesian Optimization with Agilent software and hardware. This code is less straightforward to run as it requires some changes to the specific setup at hand (Paths, Chemstation instance, Slack support, etc.). Please do not hesitate to email us if you need any clarifications.

When using this code, please cite: [![DOI:10.1016/j.aca.2023.340789](http://img.shields.io/badge/DOI-10.1016/j.aca.2023.340789-B31B1B.svg)](https://doi.org/10.1016/j.aca.2023.340789)

[cc-by-nc-sa]: http://creativecommons.org/licenses/by-nc-sa/4.0/
[cc-by-nc-sa-image]: https://licensebuttons.net/l/by-nc-sa/4.0/88x31.png
[cc-by-nc-sa-shield]: https://img.shields.io/badge/License-CC%20BY--NC--SA%204.0-lightgrey.svg
