# print message
echo "Installing all required Python packages..."

# source research config file
source ~/.research_config

# location of Python directory inside Katsevich group software directory
py_dir=$KATS_SOFT_DIR/py

# create a Python virtual environment for this project if it does not already exist
if [ ! -d $py_dir/lowmoi-venv ] 
then
    python -m venv $py_dir/lowmoi-venv
fi

# activate the virtual environment
source $py_dir/lowmoi-venv/bin/activate

# upgrade the base installer packages to latest
python -m pip install -U pip
python -m pip install -U setuptools wheel

# install required Python packages
python -m pip install numpy==1.21.5
python -m pip install scipy
python -m pip install statsmodels
python -m pip install scikit-learn
python -m pip install pandas
python -m pip install six
python -m pip install tqdm
python -m pip install joblib
python -m pip install matplotlib
python -m pip install seaborn
python -m pip install numexpr

