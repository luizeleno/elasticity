# read version from installed package
from importlib.metadata import version
__version__ = version("elasticityproject")
from elasticityproject.elasticityproject import ElasticityTheory,DirectionalProperties,VRH
from elasticityproject.plots3d import mayaviplot,LinearCompressibility3DPlot,YoungModulus3DPlot
