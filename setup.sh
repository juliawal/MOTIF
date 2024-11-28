#  Set up conda
source $(conda info --base)/etc/profile.d/conda.sh
conda create --name motif python=3.11.8
conda activate motif

#  Set the necessary environment variable for scikit-learn
export SKLEARN_ALLOW_DEPRECATED_SKLEARN_PACKAGE_INSTALL=True

#  Install the packages listed in requirements.txt
pip install -r requirements.txt