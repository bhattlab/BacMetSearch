from os.path import join, dirname, isfile

__version__ = '0.0.1_dev'

BACMET2_EXPERIMENTAL_DMND = join(dirname(__file__), 'data/BacMet2_EXP_database.dmnd')
BACMET2_PREDICTED_DMND = join(dirname(__file__), 'data/BacMet2_predicted_database.dmnd')

PRODIGAL_PATH = join(dirname(__file__), 'bin/prodigal.linux')
DIAMOND_PATH = join(dirname(__file__), 'bin/diamond')
