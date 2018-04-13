# NTA-GUI
A Graphical User Interface designed to process Non-targeted analysis of chemical samples.
Input into the GUI is a csv or tsv as output by Agilent's MPP software. The backend of
of this GUI performs a data reduction of ducplicate chemicals, calculates a number of relevant
statistical variables for every sample while combining them in replicates, performs a Toxpi calculation,
hits the comptox dashboard automatically using selenium, and processes the positive and negative mode data
in multithreads for speed and efficiency. 
