#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
This experiment was created using PsychoPy3 Experiment Builder (v2021.2.3),
    on Sun Oct 16 21:05:06 2022
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
expInfo = {'participant': '', 'task': 'calib'}
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
    originPath='/Users/psmathe/Library/CloudStorage/GoogleDrive-psmathe@g.clemson.edu/My Drive/CPSC 4120/eyetracking-experiment/PackagingExperiment_lastrun.py',
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
    text="We must first calibrate the eye tracker.\n\nKeep your eye on the calibration dot but don't anticipate it's movement.\n\nPress space to continue.",
    font='Open Sans',
    pos=(0, 0), height=0.1, wrapWidth=None, ori=0.0, 
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

# Initialize components for Routine "stimulus"
stimulusClock = core.Clock()
Healthy_food = ['stimuli/banana_chips.png', 'stimuli/fruit_cup.png']
Unhealthy_food = ['stimuli/nacho_chips.png','stimuli/potato_chips.png']
Random_images = ['stimuli/distractions/watch.jpg', 'stimuli/distractions/shoppingcart.jpg']
TopLeft = visual.ImageStim(
    win=win,
    name='TopLeft', 
    image='sin', mask=None,
    ori=0.0, pos=(0, 0), size=(0.5, 0.5),
    color=[1,1,1], colorSpace='rgb', opacity=None,
    flipHoriz=False, flipVert=False,
    texRes=128.0, interpolate=True, depth=-1.0)
TopRight = visual.ImageStim(
    win=win,
    name='TopRight', 
    image='sin', mask=None,
    ori=0.0, pos=(0, 0), size=(0.5, 0.5),
    color=[1,1,1], colorSpace='rgb', opacity=None,
    flipHoriz=False, flipVert=False,
    texRes=128.0, interpolate=True, depth=-2.0)
Bottom = visual.ImageStim(
    win=win,
    name='Bottom', 
    image='sin', mask=None,
    ori=0.0, pos=(0, 0), size=(0.5, 0.5),
    color=[1,1,1], colorSpace='rgb', opacity=None,
    flipHoriz=False, flipVert=False,
    texRes=128.0, interpolate=True, depth=-3.0)
etRecord = hardware.eyetracker.EyetrackerControl(
    server=ioServer,
    tracker=eyetracker
)

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
calib_loop = data.TrialHandler(nReps=0.0, method='sequential', 
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
    
    # -------Run Routine 'validation'-------
    
    # define target for validation
    validationTarget = visual.TargetStim(win, 
        name='validationTarget',
        radius=0.01, fillColor='', borderColor='black', lineWidth=2.0,
        innerRadius=0.0035, innerFillColor='green', innerBorderColor='black', innerLineWidth=2.0,
        colorSpace='rgb', units=None
    )
    # define parameters for validation
    validation = iohub.ValidationProcedure(win,
        target=validationTarget,
        gaze_cursor='green', 
        positions='FIVE_POINTS', randomize_positions=True,
        expand_scale=1.0, target_duration=1.0,
        enable_position_animation=True, target_delay=1.0,
        progress_on_key=None,
        show_results_screen=True, save_results_screen=False,
        color_space='rgb', unit_type=None
    )
    # run validation
    validation.run()
    # clear any keypresses from during validation so they don't interfere with the experiment
    defaultKeyboard.clearEvents()
    # the Routine "validation" was not non-slip safe, so reset the non-slip timer
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
    
# completed 0.0 repeats of 'calib_loop'


# set up handler to look after randomisation of conditions etc
stimulus_loop = data.TrialHandler(nReps=2.0, method='random', 
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
    # update component parameters for each repeat
    random.seed()
    image_num = 3
    Image_index = random.randrange(3)
    
    images = [Healthy_food[Image_index], Unhealthy_food[Image_index], Random_image[Image_index]]
    random.shuffle(images)
    
    
    
    TopLeft.setImage(images[0])
    TopRight.setImage(images[1])
    Bottom.setImage(images[2])
    # keep track of which components have finished
    stimulusComponents = [TopLeft, TopRight, Bottom, etRecord]
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
    while continueRoutine:
        # get current time
        t = stimulusClock.getTime()
        tThisFlip = win.getFutureFlipTime(clock=stimulusClock)
        tThisFlipGlobal = win.getFutureFlipTime(clock=None)
        frameN = frameN + 1  # number of completed frames (so 0 is the first frame)
        # update/draw components on each frame
        
        # *TopLeft* updates
        if TopLeft.status == NOT_STARTED and tThisFlip >= 0.0-frameTolerance:
            # keep track of start time/frame for later
            TopLeft.frameNStart = frameN  # exact frame index
            TopLeft.tStart = t  # local t and not account for scr refresh
            TopLeft.tStartRefresh = tThisFlipGlobal  # on global time
            win.timeOnFlip(TopLeft, 'tStartRefresh')  # time at next scr refresh
            TopLeft.setAutoDraw(True)
        if TopLeft.status == STARTED:
            # is it time to stop? (based on global clock, using actual start)
            if tThisFlipGlobal > TopLeft.tStartRefresh + 10.0-frameTolerance:
                # keep track of stop time/frame for later
                TopLeft.tStop = t  # not accounting for scr refresh
                TopLeft.frameNStop = frameN  # exact frame index
                win.timeOnFlip(TopLeft, 'tStopRefresh')  # time at next scr refresh
                TopLeft.setAutoDraw(False)
        
        # *TopRight* updates
        if TopRight.status == NOT_STARTED and tThisFlip >= 0.0-frameTolerance:
            # keep track of start time/frame for later
            TopRight.frameNStart = frameN  # exact frame index
            TopRight.tStart = t  # local t and not account for scr refresh
            TopRight.tStartRefresh = tThisFlipGlobal  # on global time
            win.timeOnFlip(TopRight, 'tStartRefresh')  # time at next scr refresh
            TopRight.setAutoDraw(True)
        
        # *Bottom* updates
        if Bottom.status == NOT_STARTED and tThisFlip >= 0.0-frameTolerance:
            # keep track of start time/frame for later
            Bottom.frameNStart = frameN  # exact frame index
            Bottom.tStart = t  # local t and not account for scr refresh
            Bottom.tStartRefresh = tThisFlipGlobal  # on global time
            win.timeOnFlip(Bottom, 'tStartRefresh')  # time at next scr refresh
            Bottom.setAutoDraw(True)
        # *etRecord* updates
        if etRecord.status == NOT_STARTED and t >= 0.0-frameTolerance:
            # keep track of start time/frame for later
            etRecord.frameNStart = frameN  # exact frame index
            etRecord.tStart = t  # local t and not account for scr refresh
            etRecord.tStartRefresh = tThisFlipGlobal  # on global time
            win.timeOnFlip(etRecord, 'tStartRefresh')  # time at next scr refresh
            etRecord.status = STARTED
        if etRecord.status == STARTED:
            # is it time to stop? (based on global clock, using actual start)
            if tThisFlipGlobal > etRecord.tStartRefresh + 10-frameTolerance:
                # keep track of stop time/frame for later
                etRecord.tStop = t  # not accounting for scr refresh
                etRecord.frameNStop = frameN  # exact frame index
                win.timeOnFlip(etRecord, 'tStopRefresh')  # time at next scr refresh
                etRecord.status = FINISHED
        
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
    stimulus_loop.addData('TopLeft.started', TopLeft.tStartRefresh)
    stimulus_loop.addData('TopLeft.stopped', TopLeft.tStopRefresh)
    stimulus_loop.addData('TopRight.started', TopRight.tStartRefresh)
    stimulus_loop.addData('TopRight.stopped', TopRight.tStopRefresh)
    stimulus_loop.addData('Bottom.started', Bottom.tStartRefresh)
    stimulus_loop.addData('Bottom.stopped', Bottom.tStopRefresh)
    # make sure the eyetracker recording stops
    if etRecord.status != FINISHED:
        etRecord.status = FINISHED
    # the Routine "stimulus" was not non-slip safe, so reset the non-slip timer
    routineTimer.reset()
    thisExp.nextEntry()
    
# completed 2.0 repeats of 'stimulus_loop'


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
