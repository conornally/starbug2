#python3 setup.py sdist
#pip install -e .
from setuptools import setup
if __name__=="__main__":
    setup(  scripts=["bin/starbug2","bin/starbug2-match"],
            data_files=[("extras",["extras/starbug2.completion"])],
            long_description_content_type = "text/markdown")
