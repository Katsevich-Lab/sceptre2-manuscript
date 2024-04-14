# print message
echo "Installing all required Python packages..."

# source research config file
source ~/.research_config

# location of Python directory inside Katsevich group software directory
py_dir=$KATS_SOFT_DIR/py

# create a Python virtual environment for this project if it does not already exist
if [ ! -d $py_dir/lowmoi-venv ]
then
    python3 -m venv $py_dir/lowmoi-venv
fi

# activate the virtual environment
source $py_dir/lowmoi-venv/bin/activate

# upgrade the base installer packages to latest
python3 -m pip install -U pip
python3 -m pip install -U setuptools wheel

# install required python3 packages
python3 -m pip install numpy==1.21.5
python3 -m pip install scipy
python3 -m pip install statsmodels
python3 -m pip install scikit-learn
python3 -m pip install pandas
python3 -m pip install six
python3 -m pip install tqdm
python3 -m pip install joblib
python3 -m pip install matplotlib
python3 -m pip install seaborn
python3 -m pip install numexpr
python3 -m pip install requests
