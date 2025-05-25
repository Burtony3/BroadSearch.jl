# =====================================================================
# === Basic Ephemeris Default Data

# Solar System Barycenter
ssb = BasicEphemeris(1.3211e11; name="Solar System Barycenter")

# Mercury
mercury = BasicEphemeris(
    SA[ 5.016821503445289E+07, -3.002552617224754E+07, -7.031821182995008E+06],
    SA[ 1.472812803662450E+01,  4.446008830051881E+01,  2.283740460202777E+00],
    Dates.julian2datetime(2460806.5); 
    μ=22_031.86855, 
    R=2_440.53,
    parent=ssb, 
    name="Mercury",
    id=199
)

# Venus
venus = BasicEphemeris(
    SA[-1.826584837488302E+07, -1.080324707474637E+08, -4.362302843822464E+05], 
    SA[ 3.433811751339912E+01, -5.798531552765095E+00, -2.060422993618770E+00], 
    Dates.julian2datetime(2460806.5); 
    μ=324_858.592, 
    R=6_051.893,
    parent=ssb, 
    name="Venus",
    id=299
)

# Earth
earth = BasicEphemeris(
    SA[-9.739065086101781E+07, -1.168835899939733E+08, 3.215712293333560E+04], 
    SA[ 2.241415755885093E+01, -1.916764497189185E+01, 1.647478203533836E-03], 
    Dates.julian2datetime(2460806.5); 
    μ=398_600.435436, 
    R=6_378.137,
    parent=ssb, 
    name="Earth",
    id=399
)

# Mars
mars = BasicEphemeris(
    SA[-2.427376587813935E+08,  5.648226215753347E+07,  7.159669810561325E+06], 
    SA[-4.661523228103734E+00, -2.151121188496233E+01, -3.363192259658989E-01], 
    Dates.julian2datetime(2460806.5);
    μ=42_828.375214, 
    R=3_396.19,
    parent=ssb, 
    name="Mars",
    id=499
)

# Jupiter
jupiter = BasicEphemeris(
    SA[ 9.808456526822060E+06, 7.661066057682021E+08, -3.396683274403751E+06], 
    SA[-1.321463753771015E+01, 7.881160546847092E-01,  2.924111662178345E-01], 
    Dates.julian2datetime(2460806.5);
    μ=126_686_531.900, 
    R=71_492.0,
    parent=ssb, 
    name="Jupiter",
    id=599
)

# Saturn
saturn = BasicEphemeris(
    SA[ 1.424436296067085E+09, -1.576537880843841E+08, -5.397292156412704E+07],
    SA[ 5.283134994599743E-01,  9.580258264540364E+00, -1.875537018150935E-01],
    Dates.julian2datetime(2460806.5);
    μ=37_931_206.234, 
    R=60_268.0,
    parent=ssb, 
    name="Saturn",
    id=699
)

# Uranus
uranus = BasicEphemeris(
    SA[ 1.596144490155079E+09,  2.446018815406066E+09, -1.159387277380586E+07],
    SA[-5.753214248792700E+00,  3.404148738627178E+00,  8.724873530889643E-02],
    Dates.julian2datetime(2460806.5);
    μ=5_793_950.6103, 
    R=25_559.0,
    parent=ssb, 
    name="Uranus",
    id=799
)

# Neptune
neptune = BasicEphemeris(
    SA[ 4.469599984022156E+09, -3.420369314639731E+07, -1.023021695072326E+08],
    SA[ 5.432314291365605E-03,  5.466592668941207E+00, -1.131408667304561E-01],
    Dates.julian2datetime(2460806.5);
    μ=6_835_099.97, 
    R=24_766.0,
    parent=ssb, 
    name="Neptune",
    id=899
)

# Final struct
basic = (; 
    ssb=ssb, 
    mercury=mercury, venus=venus, earth=earth, mars=mars, 
    jupiter=jupiter, saturn=saturn, uranus=uranus, neptune=neptune
)

# =====================================================================
# === Setup and exports

export DefaultEphemeris

DefaultEphemeris = (; basic=basic)

