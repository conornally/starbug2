****************
Trouble Shooting
****************

*StarbugII* crashed!
    - If the code crashes, check the output information from the routine. Warnings about incorrect inputs or bad completion of steps usually give a warning into the stderr stream with **Warning** at the start of the line. If there was no warnings given, try running the routine again with the verbose flag :code:`-v` enabled, this may show some earlier causes of problems that weren't considered bad, such as no sources passing the quality requirements.

    - Make sure the input products and parameters are formatted correctly and have been loaded into the program as intended. For example has a table been loaded into the place that an image was meant to go? Or is the parameter a string when it should be a float?

    - It may be that there was a known bug that has been fixed in an up to date version. Try updating *starbug2* with :code:`pip install -U starbug2`. 

    - If everything looks correct and it seems to be a programming error, I apologies for the inconvenience, please submit a full report on the `github <https://github.com/conornally/starbug2/issues>`_ page.

