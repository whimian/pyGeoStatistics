from .variogram_model import spherical, exponential, gaussian, power, hole_effect
from .krige3d import Krige3d
from .sgsim import Sgsim
from .gamv import Gamv

from ._version import get_versions
__version__ = get_versions()['version']
del get_versions
