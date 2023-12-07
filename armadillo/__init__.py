"""
armadillo
pi stacking project
"""

# Add imports here
#from .armadillo import canvas ##edit HNT canvas does not exist in armadillo.armadillo

# Handle versioneer
from ._version import get_versions
versions = get_versions()
__version__ = versions['version']
__git_revision__ = versions['full-revisionid']
del get_versions, versions
