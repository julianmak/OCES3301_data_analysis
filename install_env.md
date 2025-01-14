# (Brief) Installation note

Below are the notes for making a virtual environment and installing all the packages for the course OCES 3301 Data Analysis in Ocean Science.

## Dependencies for ocesda (not limited to)

- ipython
- jupyter
- notebook
- numpy
- scipy
- pandas
- matplotlib
- scikit-learn
- statsmodels
- netCDF4
- bottleneck
- dask
- xarray
- cartopy

## Google Colab

Most of the packages should be available in the temporary virtual machine (VM) that Google provides, except a few less commonly used packages (e.g. cartopy). They can be installed every time you start the kernel of the notebook, like this (using the magic IPython command `!`):
```
!pip install netCDF4 cartopy
```

Also, it is suggested to mount your google drive to the temporary VM. The mounting can be done like this:
```
from google.colab import drive
drive.mount('/content/drive')
```

## Installing environment locally

If you prefer to run the notebook locally, it is best to make an environment to 'package' and 'manage' your libraries. This is to ensure that various libraries of various versions do not crush each other between projects.

Below are two ways to make an environement.

### Using `venv`

Make sure you have python (version that is not too outdated, e.g. >= 3.9) and virtualenv installed. You can install virtualenv like this:
```
pip install virtualenv
```
or
```
pip3 install virtualenv
```
(Most of the latest python3 version should have virtualenv installed already.)

Navigate to the project repository and make the environmenet inside the directory.

Create the environment (named `ocesda`):
```
python3.10 -m venv ocesda
```

Activate the environment and install all the dependencies:
```
source ocesda/bin/activate
pip install -r requirements.txt
```


### Using `conda`

If you have Anaconda3/Miniconda installed, you can create the environment, as well as install all the required packages using the following:
```
conda env create -f environment.yml
```

Activate the environment:
```
conda activate ocesda
```

Then verify your installation:
```
conda env list
```

## Note

- All the above installation is performed on a Linux machine (Ubuntu 22.04.5 LTS) with python3.10.
- Both ways of creating environment should (ideally) work, though I am more confident with the `venv` approach.
- It is not advisable to install the dependencies in base.