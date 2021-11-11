# TiterScape - Determine COVID-19 Endpoint Titer Automatically - Stowell Lab
#### This program is written by Alex Ho. Contact him at alexho@g.harvard.edu if there are questions or bugs.
---
### Requirements:
1. Make sure Python 3 is installed; download link: https://www.python.org/downloads/release/python-383/
2. Install pip: https://www.activestate.com/resources/quick-reads/how-to-install-pip-on-windows/
3. Install related Python packages: pip install -r requirements.txt
4. TiterScape.py must appear in the same folder as the raw files (.xlsx)
---
### Execution:
1. Open the terminal in your computer by typing "cmd", "terminal", or "command prompt" in the search bar. A black screen will be prompted.
2. To execute, type "python titerscape.py <output.xlsx> <IgG_Cutoff> <IgA_Cutoff> <IgM_Cutoff>"
    * i.e. python titerscape.py Gal8.csv 0.3 0.2 0.3

---
### File Input Requirements
* Excel files only (no input number limit)
* There is some level of hard coding involved (for data cleaning purposes), which means the file input structure must match with TiterScape's desired format. Checkout "sample1.xlsx" for reference. 
    * Each excel file needs to have 24 samples and 3 antibodies in the order of IgG, IgA, and IgM
    * ID column (column C in the sample): sample names
    * Follow the input requirements video for furthur details

---
### File Output Requirements
* csv file type only
* Expected outputs:
    * one csv file
    * IgG, IgM, IgA graphs
---
### Tutorial
* Follow the tutorial video in case you have questions

