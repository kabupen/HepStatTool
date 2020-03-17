# HighEnergyPhysics Statistics interface
## How to compute the cross sction limit
1. Tune the hisgorms binning
```
root -l -b -q src/rebin.cxx
```

2. Create RooWorkspace
```
hist2workspace -standard_form /path/to/driver.xml
```

# Other tips
## Install pyhf
```
python3 -m venv --without-pip pyhf
source pyhf/bin/activate
curl https://bootstrap.pypa.io/get-pip.py | python
deactivate
source pyhf/bin/activate
```

conda create -n py27 python=2.7 anaconda
source ~/.pyenv/versions/anaconda3-4.0.0/bin/activate py27
pip install pyhf
pip install uproot
source ~/.pyenv/versions/anaconda3-4.0.0/bin/deactivate 

