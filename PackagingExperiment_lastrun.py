#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
This experiment was created using PsychoPy3 Experiment Builder (v2021.2.3),
    on November 20, 2022, at 16:18
If you publish work using this script the most relevant publication is:

    Peirce J, Gray JR, Simpson S, MacAskill M, Höchenberger R, Sogo H, Kastman E, Lindeløv JK. (2019) 
        PsychoPy2: Experiments in behavior made easy Behav Res 51: 195. 
        https://doi.org/10.3758/s13428-018-01193-y

"""

from __future__ import absolute_import, division

from psychopy import locale_setup
from psychopy import prefs
from psychopy import sound, gui, visual, core, data, event, logging, clock, colors, iohub, hardware
from psychopy.constants import (NOT_STARTED, STARTED, PLAYING, PAUSED,
                                STOPPED, FINISHED, PRESSED, RELEASED, FOREVER)

import numpy as np  # whole numpy lib is available, prepend 'np.'
from numpy import (sin, cos, tan, log, log10, pi, average,
                   sqrt, std, deg2rad, rad2deg, linspace, asarray)
from numpy.random import random, randint, normal, shuffle, choice as randchoice
import os  # handy system and path functions
import sys  # to get file system encoding

import psychopy.iohub as io
from psychopy.hardware import keyboard



# Ensure that relative paths start from the same directory as this script
_thisDir = os.path.dirname(os.path.abspath(__file__))
os.chdir(_thisDir)

# Store info about the experiment session
psychopyVersion = '2021.2.3'
expName = 'PackagingExperiment'  # from the Builder filename that created this script
expInfo = {'participant': '', 'task': 'Experiment'}
dlg = gui.DlgFromDict(dictionary=expInfo, sortKeys=False, title=expName)
if dlg.OK == False:
    core.quit()  # user pressed cancel
expInfo['date'] = data.getDateStr()  # add a simple timestamp
expInfo['expName'] = expName
expInfo['psychopyVersion'] = psychopyVersion

# Data file name stem = absolute path + name; later add .psyexp, .csv, .log, etc
filename = _thisDir + os.sep + u'data/%s_%s_%s' % (expInfo['participant'], expName, expInfo['date'])

# An ExperimentHandler isn't essential but helps with data saving
thisExp = data.ExperimentHandler(name=expName, version='',
    extraInfo=expInfo, runtimeInfo=None,
    originPath='\\\\home.clemson.edu\\psmathe\\Desktop\\CPSC 4120\\Experiment\\eyetracking-experiment\\PackagingExperiment_lastrun.py',
    savePickle=True, saveWideText=True,
    dataFileName=filename)
# save a log file for detail verbose info
logFile = logging.LogFile(filename+'.log', level=logging.EXP)
logging.console.setLevel(logging.WARNING)  # this outputs to the screen, not a file

endExpNow = False  # flag for 'escape' or other condition => quit the exp
frameTolerance = 0.001  # how close to onset before 'same' frame

# Start Code - component code to be run after the window creation

# Setup the Window
win = visual.Window(
    size=[1920, 1080], fullscr=True, screen=0, 
    winType='pyglet', allowGUI=False, allowStencil=False,
    monitor='testMonitor', color=[0,0,0], colorSpace='rgb',
    blendMode='avg', useFBO=True, 
    units='norm')
# store frame rate of monitor if we can measure it
expInfo['frameRate'] = win.getActualFrameRate()
if expInfo['frameRate'] != None:
    frameDur = 1.0 / round(expInfo['frameRate'])
else:
    frameDur = 1.0 / 60.0  # could not measure, so guess

# Setup eyetracking
ioDevice = 'eyetracker.hw.gazepoint.gp3.EyeTracker'
ioConfig = {
    ioDevice: {
        'name': 'tracker',
        'network_settings': {
            'ip_address': '127.0.0.1',
            'port': 4242.0
        }
    }
}
ioSession = '1'
if 'session' in expInfo:
    ioSession = str(expInfo['session'])
ioServer = io.launchHubServer(window=win, experiment_code='PackagingExperiment', session_code=ioSession, datastore_name=filename, **ioConfig)
eyetracker = ioServer.getDevice('tracker')

# create a default keyboard (e.g. to check for escape)
defaultKeyboard = keyboard.Keyboard()

# Initialize components for Routine "calib_instr"
calib_instrClock = core.Clock()
calib_instr_text = visual.TextStim(win=win, name='calib_instr_text',
    text="Before we begin the experiment, we must FIRST calibrate the eye tracker.\n\nKeep your eye on the calibration dot, but don't anticipate it's movement.\n\nWhen you are ready to calibrate, press the SPACEBAR.\n\n*NOTE: It may take up to 7 seconds to move onto the next screen.",
    font='Open Sans',
    pos=(0, 0), height=0.07, wrapWidth=None, ori=0.0, 
    color='white', colorSpace='rgb', opacity=None, 
    languageStyle='LTR',
    depth=0.0);
calib_instr_resp = keyboard.Keyboard()

# Initialize components for Routine "calib_cont"
calib_contClock = core.Clock()
timeout = 10
ave_error = 0.0
valid_calib_points = 0
calib_text = visual.TextStim(win=win, name='calib_text',
    text='',
    font='Open Sans',
    pos=(0, 0), height=0.1, wrapWidth=None, ori=0.0, 
    color='white', colorSpace='rgb', opacity=None, 
    languageStyle='LTR',
    depth=-1.0);
calib_resp = keyboard.Keyboard()

# Initialize components for Routine "trial_instr"
trial_instrClock = core.Clock()
text = visual.TextStim(win=win, name='text',
    text='Now we are ready to begin the actual experiment.\n\n1. On the following screens, there will be 3 different items displayed.\n2. Focus your gaze on the items that stand out the most to you. There is no right or wrong thing you should be looking at on each slide! Just simply examine whatever grabs your attention the most.\n3. Once the experiment has ended, follow the on-screen instructions.\n\nPress the SPACEBAR to continue.',
    font='Open Sans',
    pos=(0, 0), height=0.07, wrapWidth=None, ori=0.0, 
    color='white', colorSpace='rgb', opacity=None, 
    languageStyle='LTR',
    depth=0.0);
key_resp = keyboard.Keyboard()

# Initialize components for Routine "stimulus"
stimulusClock = core.Clock()
#import random
#Healthy_food = ['stimuli/banana_chips.png', 'stimuli/fruit_cup.png', 'stimuli/clementines.png', 'stimuli/health_bar.png']
#Unhealthy_food = ['stimuli/nacho_chips.png','stimuli/lays_chips.png', 'stimuli/chocolate_bar.png', 'stimuli/chocolate_cookies.png']
#Random_image = ['stimuli/distractions/watch.png', 'stimuli/distractions/cart.png', 'stimuli/distractions/scissors.png', 'stimuli/distractions/basketball.png']
#random.seed()

ioServer.sendMessageEvent(text='%s' % win.units, category='units')
trial_eyetracker = hardware.eyetracker.EyetrackerControl(
    server=ioServer,
    tracker=eyetracker
)
layout = visual.ImageStim(
    win=win,
    name='layout', units='pix', 
    image='sin', mask=None,
    ori=0.0, pos=(0, 0), size=(1920, 1080),
    color=[1,1,1], colorSpace='rgb', opacity=None,
    flipHoriz=False, flipVert=False,
    texRes=128.0, interpolate=True, depth=-2.0)

# Initialize components for Routine "closing"
closingClock = core.Clock()
closing_text = visual.TextStim(win=win, name='closing_text',
    text='Thank you for your participation in the Food Packaging Eyetracking Experiment!\n\nThe experiment is now over. Please alert one of the researchers, and be sure to fill out a POST-EXPERIMENT SURVEY at this time. \n\nWe thank you again for your humble participation to help further our eyetracking research!\n\nPress the SPACEBAR confirming you have read all of this in order to end the experiment.',
    font='Open Sans',
    pos=(0, 0), height=0.07, wrapWidth=None, ori=0.0, 
    color='white', colorSpace='rgb', opacity=None, 
    languageStyle='LTR',
    depth=0.0);
closing_key_resp = keyboard.Keyboard()

# Create some handy timers
globalClock = core.Clock()  # to track the time since experiment started
routineTimer = core.CountdownTimer()  # to track time remaining of each (non-slip) routine 

# ------Prepare to start Routine "calib_instr"-------
continueRoutine = True
# update component parameters for each repeat
calib_instr_resp.keys = []
calib_instr_resp.rt = []
_calib_instr_resp_allKeys = []
# keep track of which components have finished
calib_instrComponents = [calib_instr_text, calib_instr_resp]
for thisComponent in calib_instrComponents:
    thisComponent.tStart = None
    thisComponent.tStop = None
    thisComponent.tStartRefresh = None
    thisComponent.tStopRefresh = None
    if hasattr(thisComponent, 'status'):
        thisComponent.status = NOT_STARTED
# reset timers
t = 0
_timeToFirstFrame = win.getFutureFlipTime(clock="now")
calib_instrClock.reset(-_timeToFirstFrame)  # t0 is time of first possible flip
frameN = -1

# -------Run Routine "calib_instr"-------
while continueRoutine:
    # get current time
    t = calib_instrClock.getTime()
    tThisFlip = win.getFutureFlipTime(clock=calib_instrClock)
    tThisFlipGlobal = win.getFutureFlipTime(clock=None)
    frameN = frameN + 1  # number of completed frames (so 0 is the first frame)
    # update/draw components on each frame
    
    # *calib_instr_text* updates
    if calib_instr_text.status == NOT_STARTED and tThisFlip >= 0.0-frameTolerance:
        # keep track of start time/frame for later
        calib_instr_text.frameNStart = frameN  # exact frame index
        calib_instr_text.tStart = t  # local t and not account for scr refresh
        calib_instr_text.tStartRefresh = tThisFlipGlobal  # on global time
        win.timeOnFlip(calib_instr_text, 'tStartRefresh')  # time at next scr refresh
        calib_instr_text.setAutoDraw(True)
    
    # *calib_instr_resp* updates
    waitOnFlip = False
    if calib_instr_resp.status == NOT_STARTED and tThisFlip >= 0.0-frameTolerance:
        # keep track of start time/frame for later
        calib_instr_resp.frameNStart = frameN  # exact frame index
        calib_instr_resp.tStart = t  # local t and not account for scr refresh
        calib_instr_resp.tStartRefresh = tThisFlipGlobal  # on global time
        win.timeOnFlip(calib_instr_resp, 'tStartRefresh')  # time at next scr refresh
        calib_instr_resp.status = STARTED
        # keyboard checking is just starting
        waitOnFlip = True
        win.callOnFlip(calib_instr_resp.clock.reset)  # t=0 on next screen flip
        win.callOnFlip(calib_instr_resp.clearEvents, eventType='keyboard')  # clear events on next screen flip
    if calib_instr_resp.status == STARTED and not waitOnFlip:
        theseKeys = calib_instr_resp.getKeys(keyList=['space'], waitRelease=False)
        _calib_instr_resp_allKeys.extend(theseKeys)
        if len(_calib_instr_resp_allKeys):
            calib_instr_resp.keys = _calib_instr_resp_allKeys[-1].name  # just the last key pressed
            calib_instr_resp.rt = _calib_instr_resp_allKeys[-1].rt
            # a response ends the routine
            continueRoutine = False
    
    # check for quit (typically the Esc key)
    if endExpNow or defaultKeyboard.getKeys(keyList=["escape"]):
        core.quit()
    
    # check if all components have finished
    if not continueRoutine:  # a component has requested a forced-end of Routine
        break
    continueRoutine = False  # will revert to True if at least one component still running
    for thisComponent in calib_instrComponents:
        if hasattr(thisComponent, "status") and thisComponent.status != FINISHED:
            continueRoutine = True
            break  # at least one component has not yet finished
    
    # refresh the screen
    if continueRoutine:  # don't flip if this routine is over or we'll get a blank screen
        win.flip()

# -------Ending Routine "calib_instr"-------
for thisComponent in calib_instrComponents:
    if hasattr(thisComponent, "setAutoDraw"):
        thisComponent.setAutoDraw(False)
thisExp.addData('calib_instr_text.started', calib_instr_text.tStartRefresh)
thisExp.addData('calib_instr_text.stopped', calib_instr_text.tStopRefresh)
# the Routine "calib_instr" was not non-slip safe, so reset the non-slip timer
routineTimer.reset()

# set up handler to look after randomisation of conditions etc
calib_loop = data.TrialHandler(nReps=10.0, method='sequential', 
    extraInfo=expInfo, originPath=-1,
    trialList=[None],
    seed=None, name='calib_loop')
thisExp.addLoop(calib_loop)  # add the loop to the experiment
thisCalib_loop = calib_loop.trialList[0]  # so we can initialise stimuli with some values
# abbreviate parameter names if possible (e.g. rgb = thisCalib_loop.rgb)
if thisCalib_loop != None:
    for paramName in thisCalib_loop:
        exec('{} = thisCalib_loop[paramName]'.format(paramName))

for thisCalib_loop in calib_loop:
    currentLoop = calib_loop
    # abbreviate parameter names if possible (e.g. rgb = thisCalib_loop.rgb)
    if thisCalib_loop != None:
        for paramName in thisCalib_loop:
            exec('{} = thisCalib_loop[paramName]'.format(paramName))
    
    # -------Run Routine 'calibration'-------
    
    # define target for calibration
    calibrationTarget = visual.TargetStim(win, 
        name='calibrationTarget',
        radius=0.01, fillColor='', borderColor='black', lineWidth=2.0,
        innerRadius=0.0035, innerFillColor='green', innerBorderColor='black', innerLineWidth=2.0,
        colorSpace='rgb', units=None
    )
    # define parameters for calibration
    calibration = hardware.eyetracker.EyetrackerCalibration(win, 
        eyetracker, calibrationTarget,
        units=None, colorSpace='rgb',
        progressMode='time', targetDur=1.0, expandScale=1.5,
        targetLayout='FIVE_POINTS', randomisePos=True,
        movementAnimation=True, targetDelay=1.0
    )
    # run calibration
    calibration.run()
    # clear any keypresses from during calibration so they don't interfere with the experiment
    defaultKeyboard.clearEvents()
    # the Routine "calibration" was not non-slip safe, so reset the non-slip timer
    routineTimer.reset()
    
    # ------Prepare to start Routine "calib_cont"-------
    continueRoutine = True
    # update component parameters for each repeat
    if calibration is not None:
        print(calibration.last)
        ave_error = calibration.last['SUMMARY']['AVE_ERROR']
        valid_calib_points = calibration.last['SUMMARY']['VALID_POINTS']
    
    calib_query = "Average error is: %f\n\nValid points: %d\n\nRecalibrate (y/n)?" % (ave_error,valid_calib_points)
    calib_text.setText(calib_query)
    calib_resp.keys = []
    calib_resp.rt = []
    _calib_resp_allKeys = []
    # keep track of which components have finished
    calib_contComponents = [calib_text, calib_resp]
    for thisComponent in calib_contComponents:
        thisComponent.tStart = None
        thisComponent.tStop = None
        thisComponent.tStartRefresh = None
        thisComponent.tStopRefresh = None
        if hasattr(thisComponent, 'status'):
            thisComponent.status = NOT_STARTED
    # reset timers
    t = 0
    _timeToFirstFrame = win.getFutureFlipTime(clock="now")
    calib_contClock.reset(-_timeToFirstFrame)  # t0 is time of first possible flip
    frameN = -1
    
    # -------Run Routine "calib_cont"-------
    while continueRoutine:
        # get current time
        t = calib_contClock.getTime()
        tThisFlip = win.getFutureFlipTime(clock=calib_contClock)
        tThisFlipGlobal = win.getFutureFlipTime(clock=None)
        frameN = frameN + 1  # number of completed frames (so 0 is the first frame)
        # update/draw components on each frame
        
        # *calib_text* updates
        if calib_text.status == NOT_STARTED and tThisFlip >= 0.0-frameTolerance:
            # keep track of start time/frame for later
            calib_text.frameNStart = frameN  # exact frame index
            calib_text.tStart = t  # local t and not account for scr refresh
            calib_text.tStartRefresh = tThisFlipGlobal  # on global time
            win.timeOnFlip(calib_text, 'tStartRefresh')  # time at next scr refresh
            calib_text.setAutoDraw(True)
        
        # *calib_resp* updates
        waitOnFlip = False
        if calib_resp.status == NOT_STARTED and tThisFlip >= 0.0-frameTolerance:
            # keep track of start time/frame for later
            calib_resp.frameNStart = frameN  # exact frame index
            calib_resp.tStart = t  # local t and not account for scr refresh
            calib_resp.tStartRefresh = tThisFlipGlobal  # on global time
            win.timeOnFlip(calib_resp, 'tStartRefresh')  # time at next scr refresh
            calib_resp.status = STARTED
            # keyboard checking is just starting
            waitOnFlip = True
            win.callOnFlip(calib_resp.clock.reset)  # t=0 on next screen flip
            win.callOnFlip(calib_resp.clearEvents, eventType='keyboard')  # clear events on next screen flip
        if calib_resp.status == STARTED and not waitOnFlip:
            theseKeys = calib_resp.getKeys(keyList=['y', 'n'], waitRelease=False)
            _calib_resp_allKeys.extend(theseKeys)
            if len(_calib_resp_allKeys):
                calib_resp.keys = [key.name for key in _calib_resp_allKeys]  # storing all keys
                calib_resp.rt = [key.rt for key in _calib_resp_allKeys]
                # a response ends the routine
                continueRoutine = False
        
        # check for quit (typically the Esc key)
        if endExpNow or defaultKeyboard.getKeys(keyList=["escape"]):
            core.quit()
        
        # check if all components have finished
        if not continueRoutine:  # a component has requested a forced-end of Routine
            break
        continueRoutine = False  # will revert to True if at least one component still running
        for thisComponent in calib_contComponents:
            if hasattr(thisComponent, "status") and thisComponent.status != FINISHED:
                continueRoutine = True
                break  # at least one component has not yet finished
        
        # refresh the screen
        if continueRoutine:  # don't flip if this routine is over or we'll get a blank screen
            win.flip()
    
    # -------Ending Routine "calib_cont"-------
    for thisComponent in calib_contComponents:
        if hasattr(thisComponent, "setAutoDraw"):
            thisComponent.setAutoDraw(False)
    #calib_resp Keyboard component must have set Store: all keys
    if calib_resp is not None and 'n' in calib_resp.keys:
        calib_loop.finished = True
    calib_loop.addData('calib_text.started', calib_text.tStartRefresh)
    calib_loop.addData('calib_text.stopped', calib_text.tStopRefresh)
    # check responses
    if calib_resp.keys in ['', [], None]:  # No response was made
        calib_resp.keys = None
    calib_loop.addData('calib_resp.keys',calib_resp.keys)
    if calib_resp.keys != None:  # we had a response
        calib_loop.addData('calib_resp.rt', calib_resp.rt)
    calib_loop.addData('calib_resp.started', calib_resp.tStartRefresh)
    calib_loop.addData('calib_resp.stopped', calib_resp.tStopRefresh)
    # the Routine "calib_cont" was not non-slip safe, so reset the non-slip timer
    routineTimer.reset()
    thisExp.nextEntry()
    
# completed 10.0 repeats of 'calib_loop'


# ------Prepare to start Routine "trial_instr"-------
continueRoutine = True
# update component parameters for each repeat
key_resp.keys = []
key_resp.rt = []
_key_resp_allKeys = []
# keep track of which components have finished
trial_instrComponents = [text, key_resp]
for thisComponent in trial_instrComponents:
    thisComponent.tStart = None
    thisComponent.tStop = None
    thisComponent.tStartRefresh = None
    thisComponent.tStopRefresh = None
    if hasattr(thisComponent, 'status'):
        thisComponent.status = NOT_STARTED
# reset timers
t = 0
_timeToFirstFrame = win.getFutureFlipTime(clock="now")
trial_instrClock.reset(-_timeToFirstFrame)  # t0 is time of first possible flip
frameN = -1

# -------Run Routine "trial_instr"-------
while continueRoutine:
    # get current time
    t = trial_instrClock.getTime()
    tThisFlip = win.getFutureFlipTime(clock=trial_instrClock)
    tThisFlipGlobal = win.getFutureFlipTime(clock=None)
    frameN = frameN + 1  # number of completed frames (so 0 is the first frame)
    # update/draw components on each frame
    
    # *text* updates
    if text.status == NOT_STARTED and tThisFlip >= 0.0-frameTolerance:
        # keep track of start time/frame for later
        text.frameNStart = frameN  # exact frame index
        text.tStart = t  # local t and not account for scr refresh
        text.tStartRefresh = tThisFlipGlobal  # on global time
        win.timeOnFlip(text, 'tStartRefresh')  # time at next scr refresh
        text.setAutoDraw(True)
    
    # *key_resp* updates
    waitOnFlip = False
    if key_resp.status == NOT_STARTED and tThisFlip >= 0.0-frameTolerance:
        # keep track of start time/frame for later
        key_resp.frameNStart = frameN  # exact frame index
        key_resp.tStart = t  # local t and not account for scr refresh
        key_resp.tStartRefresh = tThisFlipGlobal  # on global time
        win.timeOnFlip(key_resp, 'tStartRefresh')  # time at next scr refresh
        key_resp.status = STARTED
        # keyboard checking is just starting
        waitOnFlip = True
        win.callOnFlip(key_resp.clock.reset)  # t=0 on next screen flip
        win.callOnFlip(key_resp.clearEvents, eventType='keyboard')  # clear events on next screen flip
    if key_resp.status == STARTED and not waitOnFlip:
        theseKeys = key_resp.getKeys(keyList=['space'], waitRelease=False)
        _key_resp_allKeys.extend(theseKeys)
        if len(_key_resp_allKeys):
            key_resp.keys = _key_resp_allKeys[-1].name  # just the last key pressed
            key_resp.rt = _key_resp_allKeys[-1].rt
            # a response ends the routine
            continueRoutine = False
    
    # check for quit (typically the Esc key)
    if endExpNow or defaultKeyboard.getKeys(keyList=["escape"]):
        core.quit()
    
    # check if all components have finished
    if not continueRoutine:  # a component has requested a forced-end of Routine
        break
    continueRoutine = False  # will revert to True if at least one component still running
    for thisComponent in trial_instrComponents:
        if hasattr(thisComponent, "status") and thisComponent.status != FINISHED:
            continueRoutine = True
            break  # at least one component has not yet finished
    
    # refresh the screen
    if continueRoutine:  # don't flip if this routine is over or we'll get a blank screen
        win.flip()

# -------Ending Routine "trial_instr"-------
for thisComponent in trial_instrComponents:
    if hasattr(thisComponent, "setAutoDraw"):
        thisComponent.setAutoDraw(False)
thisExp.addData('text.started', text.tStartRefresh)
thisExp.addData('text.stopped', text.tStopRefresh)
# check responses
if key_resp.keys in ['', [], None]:  # No response was made
    key_resp.keys = None
thisExp.addData('key_resp.keys',key_resp.keys)
if key_resp.keys != None:  # we had a response
    thisExp.addData('key_resp.rt', key_resp.rt)
thisExp.addData('key_resp.started', key_resp.tStartRefresh)
thisExp.addData('key_resp.stopped', key_resp.tStopRefresh)
thisExp.nextEntry()
# the Routine "trial_instr" was not non-slip safe, so reset the non-slip timer
routineTimer.reset()

# set up handler to look after randomisation of conditions etc
stimulus_loop = data.TrialHandler(nReps=1.0, method='random', 
    extraInfo=expInfo, originPath=-1,
    trialList=data.importConditions('stim.csv'),
    seed=None, name='stimulus_loop')
thisExp.addLoop(stimulus_loop)  # add the loop to the experiment
thisStimulus_loop = stimulus_loop.trialList[0]  # so we can initialise stimuli with some values
# abbreviate parameter names if possible (e.g. rgb = thisStimulus_loop.rgb)
if thisStimulus_loop != None:
    for paramName in thisStimulus_loop:
        exec('{} = thisStimulus_loop[paramName]'.format(paramName))

for thisStimulus_loop in stimulus_loop:
    currentLoop = stimulus_loop
    # abbreviate parameter names if possible (e.g. rgb = thisStimulus_loop.rgb)
    if thisStimulus_loop != None:
        for paramName in thisStimulus_loop:
            exec('{} = thisStimulus_loop[paramName]'.format(paramName))
    
    # ------Prepare to start Routine "stimulus"-------
    continueRoutine = True
    routineTimer.add(5.250000)
    # update component parameters for each repeat
    #Healthy_index = random.randrange(4)
    #Unhealthy_index = random.randrange(4)
    #Random_index = random.randrange(4)
    #
    #images = [Healthy_food[Healthy_index], Unhealthy_food[Healthy_index], Random_image[Random_index]]
    #random.shuffle(images)
    #
    ioServer.sendMessageEvent(text='%s' % ('beginStimuliFrame'), category='trial') 
    #ioServer.sendMessageEvent(text='top: %s, bottomLeft: %s, bottomRight: %s' % (images[2], images[0], images[1]), category='trial') 
    #
    layout.setImage(stim)
    # keep track of which components have finished
    stimulusComponents = [trial_eyetracker, layout]
    for thisComponent in stimulusComponents:
        thisComponent.tStart = None
        thisComponent.tStop = None
        thisComponent.tStartRefresh = None
        thisComponent.tStopRefresh = None
        if hasattr(thisComponent, 'status'):
            thisComponent.status = NOT_STARTED
    # reset timers
    t = 0
    _timeToFirstFrame = win.getFutureFlipTime(clock="now")
    stimulusClock.reset(-_timeToFirstFrame)  # t0 is time of first possible flip
    frameN = -1
    
    # -------Run Routine "stimulus"-------
    while continueRoutine and routineTimer.getTime() > 0:
        # get current time
        t = stimulusClock.getTime()
        tThisFlip = win.getFutureFlipTime(clock=stimulusClock)
        tThisFlipGlobal = win.getFutureFlipTime(clock=None)
        frameN = frameN + 1  # number of completed frames (so 0 is the first frame)
        # update/draw components on each frame
        # *trial_eyetracker* updates
        if trial_eyetracker.status == NOT_STARTED and t >= 0.0-frameTolerance:
            # keep track of start time/frame for later
            trial_eyetracker.frameNStart = frameN  # exact frame index
            trial_eyetracker.tStart = t  # local t and not account for scr refresh
            trial_eyetracker.tStartRefresh = tThisFlipGlobal  # on global time
            win.timeOnFlip(trial_eyetracker, 'tStartRefresh')  # time at next scr refresh
            trial_eyetracker.status = STARTED
        if trial_eyetracker.status == STARTED:
            # is it time to stop? (based on global clock, using actual start)
            if tThisFlipGlobal > trial_eyetracker.tStartRefresh + 5.25-frameTolerance:
                # keep track of stop time/frame for later
                trial_eyetracker.tStop = t  # not accounting for scr refresh
                trial_eyetracker.frameNStop = frameN  # exact frame index
                win.timeOnFlip(trial_eyetracker, 'tStopRefresh')  # time at next scr refresh
                trial_eyetracker.status = FINISHED
        
        # *layout* updates
        if layout.status == NOT_STARTED and tThisFlip >= 0.0-frameTolerance:
            # keep track of start time/frame for later
            layout.frameNStart = frameN  # exact frame index
            layout.tStart = t  # local t and not account for scr refresh
            layout.tStartRefresh = tThisFlipGlobal  # on global time
            win.timeOnFlip(layout, 'tStartRefresh')  # time at next scr refresh
            layout.setAutoDraw(True)
        if layout.status == STARTED:
            # is it time to stop? (based on global clock, using actual start)
            if tThisFlipGlobal > layout.tStartRefresh + 5-frameTolerance:
                # keep track of stop time/frame for later
                layout.tStop = t  # not accounting for scr refresh
                layout.frameNStop = frameN  # exact frame index
                win.timeOnFlip(layout, 'tStopRefresh')  # time at next scr refresh
                layout.setAutoDraw(False)
        
        # check for quit (typically the Esc key)
        if endExpNow or defaultKeyboard.getKeys(keyList=["escape"]):
            core.quit()
        
        # check if all components have finished
        if not continueRoutine:  # a component has requested a forced-end of Routine
            break
        continueRoutine = False  # will revert to True if at least one component still running
        for thisComponent in stimulusComponents:
            if hasattr(thisComponent, "status") and thisComponent.status != FINISHED:
                continueRoutine = True
                break  # at least one component has not yet finished
        
        # refresh the screen
        if continueRoutine:  # don't flip if this routine is over or we'll get a blank screen
            win.flip()
    
    # -------Ending Routine "stimulus"-------
    for thisComponent in stimulusComponents:
        if hasattr(thisComponent, "setAutoDraw"):
            thisComponent.setAutoDraw(False)
    ioServer.sendMessageEvent(text='%s' % ('endStimuliFrame'), category='trial') 
    # make sure the eyetracker recording stops
    if trial_eyetracker.status != FINISHED:
        trial_eyetracker.status = FINISHED
    stimulus_loop.addData('layout.started', layout.tStartRefresh)
    stimulus_loop.addData('layout.stopped', layout.tStopRefresh)
    thisExp.nextEntry()
    
# completed 1.0 repeats of 'stimulus_loop'


# ------Prepare to start Routine "closing"-------
continueRoutine = True
# update component parameters for each repeat
closing_key_resp.keys = []
closing_key_resp.rt = []
_closing_key_resp_allKeys = []
# keep track of which components have finished
closingComponents = [closing_text, closing_key_resp]
for thisComponent in closingComponents:
    thisComponent.tStart = None
    thisComponent.tStop = None
    thisComponent.tStartRefresh = None
    thisComponent.tStopRefresh = None
    if hasattr(thisComponent, 'status'):
        thisComponent.status = NOT_STARTED
# reset timers
t = 0
_timeToFirstFrame = win.getFutureFlipTime(clock="now")
closingClock.reset(-_timeToFirstFrame)  # t0 is time of first possible flip
frameN = -1

# -------Run Routine "closing"-------
while continueRoutine:
    # get current time
    t = closingClock.getTime()
    tThisFlip = win.getFutureFlipTime(clock=closingClock)
    tThisFlipGlobal = win.getFutureFlipTime(clock=None)
    frameN = frameN + 1  # number of completed frames (so 0 is the first frame)
    # update/draw components on each frame
    
    # *closing_text* updates
    if closing_text.status == NOT_STARTED and tThisFlip >= 0.0-frameTolerance:
        # keep track of start time/frame for later
        closing_text.frameNStart = frameN  # exact frame index
        closing_text.tStart = t  # local t and not account for scr refresh
        closing_text.tStartRefresh = tThisFlipGlobal  # on global time
        win.timeOnFlip(closing_text, 'tStartRefresh')  # time at next scr refresh
        closing_text.setAutoDraw(True)
    
    # *closing_key_resp* updates
    waitOnFlip = False
    if closing_key_resp.status == NOT_STARTED and tThisFlip >= 0.0-frameTolerance:
        # keep track of start time/frame for later
        closing_key_resp.frameNStart = frameN  # exact frame index
        closing_key_resp.tStart = t  # local t and not account for scr refresh
        closing_key_resp.tStartRefresh = tThisFlipGlobal  # on global time
        win.timeOnFlip(closing_key_resp, 'tStartRefresh')  # time at next scr refresh
        closing_key_resp.status = STARTED
        # keyboard checking is just starting
        waitOnFlip = True
        win.callOnFlip(closing_key_resp.clock.reset)  # t=0 on next screen flip
        win.callOnFlip(closing_key_resp.clearEvents, eventType='keyboard')  # clear events on next screen flip
    if closing_key_resp.status == STARTED and not waitOnFlip:
        theseKeys = closing_key_resp.getKeys(keyList=['space'], waitRelease=False)
        _closing_key_resp_allKeys.extend(theseKeys)
        if len(_closing_key_resp_allKeys):
            closing_key_resp.keys = _closing_key_resp_allKeys[-1].name  # just the last key pressed
            closing_key_resp.rt = _closing_key_resp_allKeys[-1].rt
            # a response ends the routine
            continueRoutine = False
    
    # check for quit (typically the Esc key)
    if endExpNow or defaultKeyboard.getKeys(keyList=["escape"]):
        core.quit()
    
    # check if all components have finished
    if not continueRoutine:  # a component has requested a forced-end of Routine
        break
    continueRoutine = False  # will revert to True if at least one component still running
    for thisComponent in closingComponents:
        if hasattr(thisComponent, "status") and thisComponent.status != FINISHED:
            continueRoutine = True
            break  # at least one component has not yet finished
    
    # refresh the screen
    if continueRoutine:  # don't flip if this routine is over or we'll get a blank screen
        win.flip()

# -------Ending Routine "closing"-------
for thisComponent in closingComponents:
    if hasattr(thisComponent, "setAutoDraw"):
        thisComponent.setAutoDraw(False)
thisExp.addData('closing_text.started', closing_text.tStartRefresh)
thisExp.addData('closing_text.stopped', closing_text.tStopRefresh)
# check responses
if closing_key_resp.keys in ['', [], None]:  # No response was made
    closing_key_resp.keys = None
thisExp.addData('closing_key_resp.keys',closing_key_resp.keys)
if closing_key_resp.keys != None:  # we had a response
    thisExp.addData('closing_key_resp.rt', closing_key_resp.rt)
thisExp.addData('closing_key_resp.started', closing_key_resp.tStartRefresh)
thisExp.addData('closing_key_resp.stopped', closing_key_resp.tStopRefresh)
thisExp.nextEntry()
# the Routine "closing" was not non-slip safe, so reset the non-slip timer
routineTimer.reset()

# Flip one final time so any remaining win.callOnFlip() 
# and win.timeOnFlip() tasks get executed before quitting
win.flip()

# these shouldn't be strictly necessary (should auto-save)
thisExp.saveAsWideText(filename+'.csv', delim='auto')
thisExp.saveAsPickle(filename)
logging.flush()
# make sure everything is closed down
thisExp.abort()  # or data files will save again on exit
win.close()
core.quit()
