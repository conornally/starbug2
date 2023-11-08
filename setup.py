#python3 setup.py sdist
#pip install -e .
from setuptools import setup
if __name__=="__main__":
    setup(  scripts=["bin/starbug2-match"],
            entry_points={
                "console_scripts":[
                    "starbug2 = starbug2.scripts.main:starbug_main"
                    ]
                },
            long_description_content_type = "text/markdown")
