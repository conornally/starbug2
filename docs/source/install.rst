************
Installation
************


Installing the latest release::

    $~ pip install starbug2

Installing the development version::

    $~ git clone https://github.com/conornally/starbug2.git
    $~ cd starbug2
    $~ python -m build
    $~ pip install .




Initialising
    After the package is installed, there are a few steps required to initialise *starbug2*:

    `WEBBPSF <https://github.com/spacetelescope/webbpsf>`_ is a dependency of *starbug* that has its own installation process which is not done automatically. This process is documented `here <https://webbpsf.readthedocs.io/en/latest/installation.html>`_ but requires two main steps. Download the data file on the website, named something like webbpsf-data-X.X.X.tar.gz and expand it into a directory, then append to your .bashrc (or equivalent)::

        $~ echo "export 'WEBBPSF_PATH=PATH/TO/DIRECTORY'" >> ~/.bashrc



    *Starbug2* needs to generate the WEBBPSFs, and collect some CRDS.
    It will generate these files by default into :code:`${HOME}/.local/share/starbug` however if you wish to use a different directory, set the environment variable :code:`"STARBUG_DATDIR"` to the desired destination::

        $~ echo "export 'STARBUG_DATDIR=PATH/TO/DESTINATION'" >> ~/.bashrc 
        $~ starbug2 --init 
    

Verification
    Verify the installation has been successful by running::
        
        $~ starbug2 --version

Tab Completion
    *Starbug2* has a bash completion script :code:`starbug2/extras/starbug2.completion`. This can be installed directly into :code:`/etc/bash_completion.d/` or :code:`"source /PATH/TO/COMPLETION/FILE"` can be place within your .bashrc. Unfortunately this completion script works only in bash shells.
    
Updating
    Update your installation of *StarbugII* with::

        $~ pip install -U starbug2
