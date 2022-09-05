
## HASHDEFS
MIRI=1
NIRCAM=2

## SOURCE FLAGS
SRC_GOOD=0
SRC_BAD=0x01
SRC_JMP=0x02
SRC_VAR=0x04 ##source frame mean >5% differnet than median


##DQ FLAGS
DQ_DO_NOT_USE=0x01
DQ_SATURATED =0x02
DQ_JUMP_DET  =0x04

## DEFAULT MATCHING COLS
match_cols=["RA","DEC","flag","flux","eflux"]

# ZERO POINT...
ZP={   "F070W"	:[3631,0],
       "F090W"	:[3631,0],
       "F115W"	:[3631,0],
       "F140M"	:[3631,0],
       "F150W"	:[3631,0],
       "F162M"	:[3631,0],
       "F164N"	:[3631,0],
       "F150W2" :[3631,0],
       "F182M"	:[3631,0],
       "F187N"	:[3631,0],
       "F200W"	:[3631,0],
       "F210M"	:[3631,0],
       "F212N"	:[3631,0],
       "F250M"	:[3631,0],
       "F277W"	:[3631,0],
       "F300M"	:[3631,0],
       "F322W2" :[3631,0],
       "F323N"	:[3631,0],
       "F335M"	:[3631,0],
       "F356W"	:[3631,0],
       "F360M"	:[3631,0],
       "F405N"	:[3631,0],
       "F410M"	:[3631,0],
       "F430M"	:[3631,0],
       "F444W"	:[3631,0],
       "F460M"	:[3631,0],
       "F466N"	:[3631,0],
       "F470N"	:[3631,0],
       "F480M"	:[3631,0],

       "F560W"  :[3631,0],
       "F770W"  :[3631,0],
       "F1000W" :[3631,0],
       "F1130W" :[3631,0],
       "F1280W" :[3631,0],
       "F1500W" :[3631,0],
       "F1800W" :[3631,0],
       "F2100W" :[3631,0],
       "F2550W" :[3631,0],
       }

#            name    wavelen FWHM(arcsec,pix)   zp zper  instrument
filters={   "F070W"	:[0.704,  0.030, 0.987,     0, 0, NIRCAM],
            "F090W"	:[0.902,  0.034, 1.103,     0, 0, NIRCAM],
            "F115W"	:[1.154,  0.040, 1.298,     0, 0, NIRCAM],
            "F140M"	:[1.405,  0.048, 1.553,     0, 0, NIRCAM],
            "F150W"	:[1.501,  0.050, 1.628,     0, 0, NIRCAM],
            "F162M"	:[1.627,  0.055, 1.770,     0, 0, NIRCAM],
            "F164N"	:[1.645,  0.056, 1.801,     0, 0, NIRCAM],
            "F150W2":[1.659,  0.046, 1.494,     0, 0, NIRCAM],
            "F182M"	:[1.845,  0.062, 1.990,     0, 0, NIRCAM],
            "F187N"	:[1.874,  0.064, 2.060,     0, 0, NIRCAM],
            "F200W"	:[1.989,  0.066, 2.141,     0, 0, NIRCAM],
            "F210M"	:[2.095,  0.071, 2.304,     0, 0, NIRCAM],
            "F212N"	:[2.121,  0.072, 2.341,     0, 0, NIRCAM],
            "F250M"	:[2.503,  0.084, 1.340,     0, 0, NIRCAM],
            "F277W"	:[2.762,  0.091, 1.444,     0, 0, NIRCAM],
            "F300M"	:[2.989,  0.100, 1.585,     0, 0, NIRCAM],
            "F322W2":[3.232,  0.097, 1.547,     0, 0, NIRCAM],
            "F323N"	:[3.237,  0.108, 1.711,     0, 0, NIRCAM],
            "F335M"	:[3.362,  0.111, 1.760,     0, 0, NIRCAM],
            "F356W"	:[3.568,  0.115, 1.830,     0, 0, NIRCAM],
            "F360M"	:[3.624,  0.120, 1.901,     0, 0, NIRCAM],
            "F405N"	:[4.052,  0.136, 2.165,     0, 0, NIRCAM],
            "F410M"	:[4.082,  0.137, 2.179,     0, 0, NIRCAM],
            "F430M"	:[4.281,  0.145, 2.300,     0, 0, NIRCAM],
            "F444W"	:[4.408,  0.145, 2.302,     0, 0, NIRCAM],
            "F460M"	:[4.630,  0.155, 2.459,     0, 0, NIRCAM],
            "F466N"	:[4.654,  0.158, 2.507,     0, 0, NIRCAM],
            "F470N"	:[4.708,  0.160, 2.535,     0, 0, NIRCAM],
            "F480M"	:[4.874,  0.162, 2.574,     0, 0, NIRCAM],

            "F560W" :[5.589,  0.182, 1.636,     0, 0, MIRI],
            "F770W" :[7.528,  0.243, 2.187,     0, 0, MIRI],
            "F1000W":[9.883,  0.321, 2.888,     0, 0, MIRI],
            "F1130W":[11.298, 0.368, 3.318,     0, 0, MIRI],
            "F1280W":[12.712, 0.412, 3.713,     0, 0, MIRI],
            "F1500W":[14.932, 0.483, 4.354,     0, 0, MIRI],
            "F1800W":[17.875, 0.580, 5.224,     0, 0, MIRI],
            "F2100W":[20.563, 0.665, 5.989,     0, 0, MIRI],
            "F2550W":[25.147, 0.812, 7.312,     0, 0, MIRI],
            }

