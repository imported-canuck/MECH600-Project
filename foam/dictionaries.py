"""Write steady incompressible OpenFOAM dictionaries for the airfoil study."""

from math import cos, isfinite, radians, sin, sqrt
from pathlib import Path
from textwrap import indent, dedent


MODELS = {
    "SpalartAllmaras": ("nuTilda",),
    "kOmegaSST": ("k", "omega"),
    "kOmegaSSTLM": ("k", "omega", "gammaInt", "ReThetat"),
}
FLAVORS = ("openfoam-com", "foundation")


def _scalar(value):
    return f"{value:.12g}"


def _vector(*values):
    return "(" + " ".join(_scalar(value) for value in values) + ")"


def _write(path, body, cls="dictionary", object_name=None):
    header = f"""FoamFile
{{
    version     2.0;
    format      ascii;
    class       {cls};
    location    "{path.parent.name}";
    object      {object_name or path.name};
}}

"""
    with path.open("w", encoding="ascii", newline="\n") as file:
        file.write(header + dedent(body).strip() + "\n")


def _patch(name, entries):
    return f"{name}\n{{\n{indent(entries, '    ')}\n}}\n"


def _field(case_dir, name, dimensions, initial, wall, exterior, vector=False):
    body = f"""dimensions      [{dimensions}];
internalField   uniform {initial};

boundaryField
{{
{indent(_patch("airfoil", wall), '    ')}
{indent(_patch("farfield", exterior), '    ')}
{indent(_patch("frontAndBack", "type empty;").rstrip(), '    ')}
}}
"""
    _write(case_dir / "0" / name, body, "volVectorField" if vector else "volScalarField")


def _freestream(patch_type, value):
    return f"type {patch_type};\nfreestreamValue uniform {value};\nvalue uniform {value};"


def _with_value(patch_type, value):
    return f"type {patch_type};\nvalue uniform {value};"


def _length_scale(spec):
    length_scale = spec.get("length_scale")
    return 0.05 * spec["chord"] if length_scale is None else length_scale


def _initial_fields(case_dir, spec, direction):
    speed = spec["speed"]
    nu = spec["nu"]
    velocity = _vector(speed * direction[0], speed * direction[1], 0)
    _field(case_dir, "U", "0 1 -1 0 0 0 0", velocity, "type noSlip;",
           _freestream("freestreamVelocity", velocity), vector=True)
    _field(case_dir, "p", "0 2 -2 0 0 0 0", "0", "type zeroGradient;",
           _freestream("freestreamPressure", "0"))

    model = spec["model"]
    wall_function = spec["wall_treatment"] == "function"
    nut_type = "nutkWallFunction" if wall_function else "nutLowReWallFunction"
    if model == "SpalartAllmaras":
        nu_tilda = 3 * nu
        nut = nu_tilda * 3**3 / (3**3 + 7.1**3)
        _field(case_dir, "nuTilda", "0 2 -1 0 0 0 0", _scalar(nu_tilda),
               _with_value("fixedValue", "0"), _freestream("freestream", _scalar(nu_tilda)))
    else:
        tu = spec.get("tu")
        if tu is None:
            # TMR SSTm farfield: 0.052% intensity and nut/nu = 0.009.
            k = 1.5 * (speed * 0.00052)**2
            omega = k / (0.009 * nu)
        else:
            k = 1.5 * (speed * tu / 100)**2
            omega = sqrt(k) / (0.09**0.25 * _length_scale(spec))
        nut = k / omega
        k_wall_type = "kqRWallFunction" if wall_function else "fixedValue"
        _field(case_dir, "k", "0 2 -2 0 0 0 0", _scalar(k),
               _with_value(k_wall_type, "0"), _freestream("freestream", _scalar(k)))
        _field(case_dir, "omega", "0 0 -1 0 0 0 0", _scalar(omega),
               _with_value("omegaWallFunction", _scalar(omega)),
               _freestream("freestream", _scalar(omega)))
        if model == "kOmegaSSTLM":
            if tu <= 1.3:
                re_theta = 1173.51 - 589.428 * tu + 0.2196 / tu**2
            else:
                re_theta = 331.5 / (tu - 0.5658)**0.671
            for name, value in (("gammaInt", "1"), ("ReThetat", _scalar(re_theta))):
                _field(case_dir, name, "0 0 0 0 0 0 0", value, "type zeroGradient;",
                       _freestream("freestream", value))
    _field(case_dir, "nut", "0 2 -1 0 0 0 0", _scalar(nut),
           _with_value(nut_type, "0"), _with_value("calculated", _scalar(nut)))


def _properties(case_dir, spec, flavor):
    if flavor == "foundation":
        transport_name, model_name, ras_key = "physicalProperties", "momentumTransport", "model"
        transport = "viscosityModel constant;\nrho [1 -3 0 0 0 0 0] 1;"
        print_coeffs = ""
    else:
        transport_name, model_name, ras_key = "transportProperties", "turbulenceProperties", "RASModel"
        transport = "transportModel Newtonian;"
        print_coeffs = "    printCoeffs on;\n"
    _write(case_dir / "constant" / transport_name, f"{transport}\nnu [0 2 -1 0 0 0 0] {_scalar(spec['nu'])};")
    _write(case_dir / "constant" / model_name, f"""simulationType RAS;

RAS
{{
    {ras_key} {spec["model"]};
    turbulence on;
{print_coeffs}}}
""")


def _schemes(spec, startup):
    if startup:
        velocity = "bounded Gauss upwind"
        scalar = "bounded Gauss upwind"
    else:
        velocity = "bounded Gauss limitedLinearV 1"
        scalar = "bounded Gauss limitedLinear 1"
    advected = "\n".join(f"    div(phi,{name}) {scalar};" for name in MODELS[spec["model"]])
    return f"""ddtSchemes
{{
    default steadyState;
}}
gradSchemes
{{
    default Gauss linear;
    grad(U) cellLimited Gauss linear 1;
}}
divSchemes
{{
    default none;
    div(phi,U) {velocity};
{advected}
    div(div(phi,U)) Gauss linear;
    div((nuEff*dev2(T(grad(U))))) Gauss linear;
}}
laplacianSchemes
{{
    default Gauss linear limited 0.5;
}}
interpolationSchemes
{{
    default linear;
}}
snGradSchemes
{{
    default limited 0.5;
}}
wallDist
{{
    method meshWave;
}}
"""


def _solution(spec):
    fields = "|".join(("U", *MODELS[spec["model"]]))
    return f"""solvers
{{
    p
    {{
        solver GAMG;
        tolerance 1e-9;
        relTol 0.02;
        smoother GaussSeidel;
    }}

    Phi
    {{
        solver GAMG;
        tolerance 1e-8;
        relTol 0;
        smoother GaussSeidel;
        maxIter 400;
    }}

    "({fields})"
    {{
        solver smoothSolver;
        smoother symGaussSeidel;
        nSweeps 2;
        tolerance 1e-9;
        relTol 0.05;
        maxIter 100;
    }}
}}

SIMPLE
{{
    nNonOrthogonalCorrectors 1;
    pRefCell 0;
    pRefValue 0;
}}

potentialFlow
{{
    nNonOrthogonalCorrectors 2;
    PhiRefCell 0;
    PhiRefValue 0;
}}

relaxationFactors
{{
    fields
    {{
        p 0.15;
    }}
    equations
    {{
        "({fields})" 0.35;
    }}
}}
"""


def _control(spec, flavor, direction, end_time, write_interval):
    application = "application simpleFoam;"
    if flavor == "foundation":
        application = "application foamRun;\nsolver incompressibleFluid;"
    chord, span, speed = spec["chord"], spec["span"], spec["speed"]
    cosine, sine = direction
    point = _vector(-2 * chord * cosine, -2 * chord * sine, 0.5 * span)
    probe_fields = "U" if spec["model"] == "SpalartAllmaras" else "U k omega"
    surface_select = "patch farfield;" if flavor == "foundation" else "regionType patch;\n        name farfield;"
    field_lib = '("libfieldFunctionObjects.so")'
    return f"""{application}

startFrom latestTime;
startTime 0;
stopAt endTime;
endTime {end_time};
deltaT 1;

writeControl runTime;
writeInterval {write_interval};
purgeWrite 3;
writeFormat ascii;
writePrecision 12;
writeCompression off;
timeFormat general;
timePrecision 12;
runTimeModifiable true;

functions
{{
    forceCoeffs
    {{
        type forceCoeffs;
        libs ("libforces.so");
        patches (airfoil);
        rho rhoInf;
        rhoInf 1;
        pRef 0;
        CofR ({_scalar(chord / 4)} 0 0);
        dragDir {_vector(cosine, sine, 0)};
        liftDir {_vector(-sine, cosine, 0)};
        pitchAxis (0 0 1);
        magUInf {_scalar(speed)};
        lRef {_scalar(chord)};
        Aref {_scalar(chord * span)};
        writeControl timeStep;
        writeInterval 5;
        writePrecision 12;
        log false;
    }}

    forces
    {{
        type forces;
        libs ("libforces.so");
        patches (airfoil);
        rho rhoInf;
        rhoInf 1;
        pRef 0;
        CofR ({_scalar(chord / 4)} 0 0);
        writeControl timeStep;
        writeInterval 5;
        writePrecision 12;
        log false;
    }}

    yPlus
    {{
        type yPlus;
        libs {field_lib};
        executeControl writeTime;
        writeControl writeTime;
        writeFields true;
    }}

    wallShearStress
    {{
        type wallShearStress;
        libs {field_lib};
        patches (airfoil);
        executeControl writeTime;
        writeControl writeTime;
        writeFields true;
    }}

    upstream
    {{
        type probes;
        libs ("libsampling.so");
        fields ({probe_fields});
        probeLocations ({point});
        interpolationScheme cell;
        writeControl timeStep;
        writeInterval 10;
    }}

    netFlux
    {{
        type surfaceFieldValue;
        libs {field_lib};
        {surface_select}
        operation sum;
        fields (phi);
        writeFields false;
        writeControl timeStep;
        writeInterval 10;
        log false;
    }}

    absoluteFlux
    {{
        type surfaceFieldValue;
        libs {field_lib};
        {surface_select}
        operation sumMag;
        fields (phi);
        writeFields false;
        writeControl timeStep;
        writeInterval 10;
        log false;
    }}
}}
"""


def _validate(spec, flavor, nprocs, end_time, write_interval):
    if flavor not in FLAVORS:
        raise ValueError(f"Unsupported OpenFOAM flavor: {flavor!r}. Choose from {FLAVORS}.")
    if spec["model"] not in MODELS:
        raise ValueError(f"Unsupported turbulence model: {spec['model']!r}.")
    if spec["wall_treatment"] not in ("resolved", "function"):
        raise ValueError("wall_treatment must be 'resolved' or 'function'.")
    if spec["wall_treatment"] == "function" and spec["model"] != "kOmegaSST":
        raise ValueError("Wall-function cases use kOmegaSST.")
    for key in ("chord", "span", "speed", "nu"):
        if not isfinite(spec[key]) or spec[key] <= 0:
            raise ValueError(f"{key} must be positive and finite.")
    if not isfinite(spec["alpha"]):
        raise ValueError("alpha must be finite.")
    tu = spec.get("tu")
    if tu is not None and (not isfinite(tu) or tu <= 0):
        raise ValueError("tu must be positive and finite, in percent.")
    if spec["model"] == "kOmegaSSTLM" and tu is None:
        raise ValueError("A value of tu in percent is required for kOmegaSSTLM.")
    length_scale = _length_scale(spec)
    if not isfinite(length_scale) or length_scale <= 0:
        raise ValueError("length_scale must be positive and finite, in meters.")
    if min(nprocs, end_time, write_interval) < 1:
        raise ValueError("Processor count and time controls must be positive.")


def write_case(case_dir, spec, flavor="openfoam-com", nprocs=1, end_time=5000, write_interval=250):
    """Prepare field and solver dictionaries, leaving mesh generation to the runner."""
    _validate(spec, flavor, nprocs, end_time, write_interval)
    case_dir = Path(case_dir)
    for directory in ("0", "constant", "system"):
        (case_dir / directory).mkdir(parents=True, exist_ok=True)

    alpha = radians(spec["alpha"])
    direction = (cos(alpha), sin(alpha))
    _initial_fields(case_dir, spec, direction)
    _properties(case_dir, spec, flavor)
    system = case_dir / "system"
    _write(system / "controlDict", _control(spec, flavor, direction, end_time, write_interval))
    _write(system / "fvSolution", _solution(spec))
    final = _schemes(spec, startup=False)
    _write(system / "fvSchemes", final)
    _write(system / "fvSchemes.final", final, object_name="fvSchemes")
    _write(system / "fvSchemes.startup", _schemes(spec, startup=True), object_name="fvSchemes")
    _write(system / "decomposeParDict", f"numberOfSubdomains {nprocs};\nmethod scotch;")
