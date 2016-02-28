#! usr/bin/python

"""
    This file contains paths to different executables 
    and template files, so that this software can be 
    configured easily on different machines.
"""


# Local Paths
Aromodel_Path = "/Users/Sam/Desktop/Research/Code/2016/Aromodel/"
Template_Path = Aromodel_Path + "/Templates/"


# Remote Paths
Comet_Login = "seroot@comet.sdsc.edu"
Comet_Path = "/oasis/scratch/comet/seroot/temp_project/Aromodel/%s" # % Directory_Name
Orca_Path = "/oasis/scratch/comet/seroot/temp_project/Orca/orca_3_0_3_linux_x86-64/orca"


c2c = "scp %s " + Comet_Login + ":" + Comet_Path # % (File, Directory_Name)
c2l = "scp " + Comet_Login + ":" + Comet_Path + "/%s ./" # % (Directory_Name, File)
SBATCH = "sbatch " + Comet_Path + "/%s" # %( Directory_Name, Submit_Script)

