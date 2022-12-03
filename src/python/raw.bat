set PYTHON=python

# use Butterworth?
set SMOOTH=False

set WIDTH=1920
set HEIGHT=1080
set HERTZ=60
set DIST=21.65
set SCREEN=17

set DFWIDTH=5
set DFDEGREE=3

set BASELINE_T=2.0
set END_T=180.0
REM what I had for testing but we reverted back to 100.0
REM VT=5.0  # more fixations
set VT=10.0
REM VT=10.0
REM VT=44.0
REM VT=80.0
REM VT=100.0
REM VT=200.0
REM VT=240.0 # fewer fixations -- good results here

SET XTILES=8
SET YTILES=4

set INDIR=../../exp/data/
set AOIDIR=../../aois/1920x1080/
set IMGDIR=../../exp/stimuli/1920x1080/

set PLTDIR=./plots/
set OUTDIR=./data/
set RAWDIR=./data/raw/

REM raw

%PYTHON% ./hdf52csv.py --indir=%INDIR% --outdir=%RAWDIR% --width=%WIDTH% --height=%HEIGHT% --dist=%DIST%
