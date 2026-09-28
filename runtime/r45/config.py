"""The single supported R installation for this repository."""
from pathlib import Path

R_HOME = Path("D:/R/R-4.5.0")
R_VERSION = "4.5.0"
RSCRIPT = R_HOME / "bin/Rscript.exe"
R_LIBRARY = R_HOME / "library"
REPO = Path(__file__).resolve().parents[2]
PROFILE = Path(__file__).with_name("r45_profile.R")
RUNTIME_LOG = REPO / ".runtime/r45"
