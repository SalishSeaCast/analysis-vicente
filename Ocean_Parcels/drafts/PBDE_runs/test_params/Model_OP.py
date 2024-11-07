############################## THIS IS THE FUNCTION FOR RUNNING THE MODEL WITH SPECIFIED PARAMETERS ################################
import sys
import random
import xarray as xr
import numpy as np
import os
import yaml
import math
from datetime import datetime, timedelta
from parcels import FieldSet, Field, VectorField, ParticleSet, JITParticle, ParcelsRandom, Variable, Kernel, AdvectionRK4
from cartopy import crs, feature
import zarr 
#
sys.path.append('/ocean/vvalenzuela/MOAD/Ocean_Parcels')
#
from OP_functions import *
#
import warnings
warnings.simplefilter(action='ignore', category=FutureWarning)


def model_simulation(inarg=None):
    # Model based in Ocean Parcels for simulations of PBDEs discharge from the Iona Outfall


    return