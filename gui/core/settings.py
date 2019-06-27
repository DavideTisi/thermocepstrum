'''
--------------------------------------------
    Thermocepstrum graphic user interface

    Settings file
--------------------------------------------

This file contains the main settings and preferences of the Thermocepstrum GUI
'''

import os

# -------- FILE SETTINGS --------

BASE_PATH = '.'
DATA_PATH = ''
OUTPUT_PATH = ''
LOG_PATH = ''
ASSETS_PATH = os.path.join(BASE_PATH, 'assets')

# todo: Add/Remove extensions
FILE_EXTENSIONS = ['.dat', '.log', '.txt', '.bin']


# -------- GUI SETTINGS --------

STATUS = ['Loading', 'Configuring', 'Calculating', 'Idle', 'Sleep']
STATUS_NOW = STATUS[0]

# LANGUAGE = 'en-EN'

X_RESIZE = True
Y_RESIZE = True

X_SIZE = 1024
Y_SIZE = 450

X_SPACING = 300
Y_SPACING = 200


# -------- GRAPHIC SETTINGS --------

BG_COLOR = '#ffffff'
FONT = 'Arial'
FONT_SIZE = 18


