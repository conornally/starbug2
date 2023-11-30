import os
from starbug2 import param

def test_parse_param():
    assert type(param.parse_param("A=1//."))==dict

    assert param.parse_param("A = 1 //.")=={'A':1}
    assert param.parse_param("A = B //.")=={'A':'B'}
    assert param.parse_param("A = B //.\n")=={'A':'B'}

    assert param.parse_param("A = //.\n")=={'A':''}
    assert param.parse_param(" = //.")=={}

    assert param.parse_param("A=B")=={"A":"B"}
    assert param.parse_param("A=B/")=={"A":"B/"}
    assert param.parse_param("A=B/.")=={"A":"B/."}
    assert param.parse_param("A=1/.")=={"A":"1/."}

    assert param.parse_param("A      =1")=={"A":1}
    assert param.parse_param("A=1      ")=={"A":1}
    assert param.parse_param("A=1     a")=={"A":"1     a"}

def test_load_default_params():
    assert param.load_default_params()!={}
    assert type(param.load_default_params()) == dict
    assert "PARAM" in param.load_default_params().keys()
    assert param.load_default_params().get("PARAM")=="STARBUGII PARAMETERS"


def test_load_params():
    assert param.load_default_params() == param.load_params(None)

    assert param.load_params("doesnotexist")=={}

    os.system("starbug2 --local-param")
    assert param.load_params("starbug.param")!={}
    assert "PARAM" in param.load_params("starbug.param").keys()
    assert param.load_params("starbug.param").get("PARAM")=="STARBUGII PARAMETERS"
    os.remove("starbug.param")

def test_update_params():
    os.system("starbug2 --local-param")
    os.system("sed -i s/PARAM/PARAM1/g starbug.param")

    assert "PARAM" not in param.load_params("starbug.param").keys()
    assert param.update_paramfile("starbug.param") is None
    assert param.update_paramfile("starbug.param") is None
    os.remove("starbug.param")

