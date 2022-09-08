#python3 setup.py sdist
#pip install -e .
from setuptools import setup, find_packages
VERSION="0.2.2"

setup(
        name="starbug2",
        version=VERSION,
        author="Conor Nally",
        author_email="conor.nally@ed.ac.uk",
        packages=find_packages(),
        include_package_data=True,
        scripts=["bin/starbug2", "bin/starbug2-match"],
        #package_data={ "":["default.param"]},
        #data_files=[ ("dat", ["starbug2/dat/default.param"]) ],
        data_files=[("/etc/bash_completion.d",["extras/starbug2.completion"]) ],
        #url="https://conornally.github.io/starbug2",
        license="LICENSE.txt",
        description="JWST PSF photometry in complex crowded fields",
        long_description=open("README.md").read(),
        long_description_content_type="text/markdown",
        install_requires=[
            "photutils",
            "astropy",
            "parse",
            "scipy",
            "webbpsf"
            ],
        )
