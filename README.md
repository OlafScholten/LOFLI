# LOFLI 2.0

This update contains a new folder organization, as well as some bug fixws in the main program.
The (somewhat outdated, still for the old organization of folders) instruction can be found in 
"Nts_LofarImaging.pdf" residing in folder 'docu'
Before running anything you will need to make a subfolder  bin  in the mail LOFLI folder
You will need to edit the initial part of "ShortCuts.sh" to point to the right directories for libraries on your machine.
Please look in the instructions on how to install GLE.
It is advisible to add the following line to your  .bashrc  script
   export LL_BaseDir=/home/olaf/LOFLI  
Note that this may not yet have been included correctly in all running scripts.

To get started on analyzing a flash, use the script "NewFlas.sh" in the main directory, or copy all of folder "flash" to the new working directory (including subfolders!)
